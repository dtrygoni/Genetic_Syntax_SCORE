# Genetic SYNTAX Score (GESS)

R code and fitted model objects accompanying the study **“Predicting coronary artery disease severity through genomic profiling and machine learning modelling.”**

The GEnetic SYNTAX Score (GESS) framework is an exploratory hierarchical machine-learning pipeline that combines clinical and laboratory predictors with candidate single-nucleotide polymorphisms (SNPs) to predict the presence and severity of coronary artery disease as represented by the SYNTAX score.

> **Status and intended use:** This repository documents internally validated model development in a single-center cohort. The models have not undergone external or prospective clinical validation and are not a clinical decision-support tool. They must not be used to select treatment or guide patient care.

## Analytical framework

The pipeline contains three sequential prediction components:

1. **Zero/non-zero component:** classifies SYNTAX score as 0 versus >0.
2. **Positive-score severity classifier:** classifies patients predicted to have a positive score as low (<23) versus intermediate/high (≥23) complexity.
3. **Stratum-specific regression:** applies the corresponding low- or intermediate/high-complexity regressor to generate the final SYNTAX-score prediction.

Two predictor configurations are compared using the same cross-validation partitions:

- **Model A:** clinical and laboratory predictors.
- **Model B:** the corresponding clinical and laboratory predictors plus SNP-derived predictors.

The primary evaluation uses pooled out-of-fold predictions from stratified 10-fold cross-validation. The repository also contains repeated cross-validation analyses and sensitivity analyses using fold-specific imputation and simpler logistic/linear models.

## Repository contents

| Path | Purpose |
|---|---|
| `Main_model/` | Primary complete-case analysis: SNP quality control, preprocessing, feature selection, model fitting, component-level evaluation, complete-pipeline evaluation, calibration/decision-curve analyses, missingness summaries, and subgroup analyses. |
| `Repeated_cross_validation/` | Ten repetitions of stratified 10-fold cross-validation, including saved fold indices and fitted fold-specific objects. |
| `Sensitivity_Analysis_Imputation/` | Sensitivity analysis using fold-specific median imputation for continuous predictors and mode imputation for categorical/discrete predictors and SNPs. |
| `Sensitivity_Analysis_SimpleModels/` | Sensitivity analysis using logistic-regression classifiers and linear-regression models as transparent comparators. |
| `GESS_model/` | Fitted full-data model objects, preprocessing objects, predictor schema, and a prediction wrapper for structured new data. |
| `Utils/` | Supporting SNP-level and performance-plotting functions. |
| `sessionInfo_complete.txt` | R, operating-system, and installed-package information recorded for the repository. |

Saved `.Rda` files contain fitted models, fold assignments, preprocessing objects, feature selections, or intermediate analysis results. They do not constitute a public release of the participant-level development dataset.

## Data availability and privacy

The participant-level clinical and genomic data are **not included** in this public repository. Unrestricted release is not permitted because the dataset contains sensitive health and genetic information and is subject to participant-consent, ethics, and data-protection restrictions.

Accordingly:

- the public repository permits inspection of the analytical logic, model objects, seeds, and fold-specific procedures;
- complete end-to-end regeneration of the reported results requires authorized access to the original data;
- access requests, if permitted under the study governance framework, must be reviewed by the responsible institution and ethics/data-access authorities.

The expected variables and coding for the fitted prediction implementation are described in `GESS_model/Dataset_schema.txt`. Researchers receiving controlled access should retain the original variable names and factor/genotype coding.

## Software environment

The project is implemented in R. The exact environment recorded for this repository is reported in `sessionInfo_complete.txt`. The principal packages used across the analysis include:

```r
c(
  "Boruta", "caret", "DescTools", "dplyr", "ggplot2",
  "kernelshap", "missForestPredict", "MLmetrics", "mltools",
  "pROC", "PRROC", "precrec", "randomForest", "recipes",
  "shapviz", "SNPassoc", "tidyr", "xgboost"
)
```

Package requirements differ across scripts. Install missing packages before running the relevant analysis. For example:

```r
packages <- c(
  "Boruta", "caret", "dplyr", "ggplot2", "MLmetrics",
  "pROC", "PRROC", "randomForest", "recipes", "SNPassoc", "tidyr"
)

to_install <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
```

Bioconductor packages, if required for the bioinformatic analyses, should be installed through `BiocManager`.

## Reproducing the development analyses

### Important portability note

The analysis scripts are preserved as research code and currently contain absolute Windows paths referring to the authors’ local analysis directories. They are therefore **not directly portable after cloning**. Before running them with an authorized copy of the data:

1. Replace absolute paths such as `C:/Users/.../GESS_experiments_code/` with a project-root variable or relative paths.
2. Store the restricted analytical dataset outside version control.
3. Confirm that loading the source file creates an R object named `dataset` with the variables and coding expected by the scripts.
4. Create separate output directories or modify `save()` calls as required.
5. Run the pipeline in a clean R session and retain the generated log and session information.

The primary scripts expect an R data file commonly referenced as `Whole_dataset.Rda`; this restricted file is not included.

### Primary analysis

The primary analysis is organized into the following stages:

1. **Outer-fold construction and zero/non-zero component**
   - `Main_model/Zero_part_classification/Cross_validation_construction.R`
   - `Main_model/Zero_part_classification/ML_zero_model_train.R`
   - `Main_model/Zero_part_classification/Preparation_before_predictions.R`
   - `Main_model/Zero_part_classification/CI_ROC.R`

2. **Positive-score severity component**
   - `Main_model/Count_part_regression/CV_count_classifier.R`
   - `Main_model/Count_part_regression/Classifiers_model_train.R`
   - the component-specific evaluation scripts in the same directory

3. **Complete hierarchical pipeline**
   - scripts in `Main_model/Whole_model_experiments/`

4. **Additional analyses**
   - `Main_model/Brier_Calibration_DCA_ConfusionMatrices_ZeroPart.R`
   - `Main_model/Subgroup_analysis_of_complete_pipeline.R`
   - `Main_model/Missing_values_predictors.R`
   - `Main_model/Comparison_simple_with_RF_fullmodel.R`

5. **Full-data refit for exploratory interpretation and fitted prediction objects**
   - scripts in `Main_model/Full_set_models/`

Full-data refits are used for fitted-model construction and exploratory interpretation; internally validated performance estimates are obtained from out-of-fold predictions rather than from the full-data refit.

### Repeated cross-validation

`Repeated_cross_validation/Created_repeated_crossvalidation_inds.R` uses `caret::createMultiFolds()` with `k = 10`, `times = 10`, and seed 123 to construct repeated stratified folds. The remaining top-level scripts in `Repeated_cross_validation/` fit and evaluate the component and complete-pipeline models across those partitions. The `Fold1/`–`Fold10/` directories contain repetition-specific indices and fitted objects.

### Missing-data sensitivity analysis

`Sensitivity_Analysis_Imputation/` mirrors the principal pipeline while replacing complete-case deletion with fold-specific imputation. Imputation parameters are estimated from each outer training fold and then applied unchanged to its corresponding test fold. Continuous predictors are imputed using the training-fold median; categorical/discrete predictors and SNPs are imputed using the training-fold mode.

### Simpler-model sensitivity analysis

`Sensitivity_Analysis_SimpleModels/` evaluates conventional logistic-regression classifiers and linear-regression models using the same analytical structure as the principal pipeline.

## Preprocessing and leakage control

For internally validated analyses, supervised preprocessing is performed within the outer training data. This includes genotype quality-control decisions, inheritance-mode selection, LD-based reduction with outcome-informed selection among correlated variants, univariate screening, Boruta feature selection, SNP-interaction construction, preprocessing, inner-cross-validation tuning, and classification-threshold selection. Learned transformations are then applied to the corresponding held-out fold. Model A and Model B use identical outer test folds, and paired comparisons use participants with predictions available from both models.

The primary analysis uses complete cases for the fold-specific selected predictors. The imputation sensitivity analysis instead estimates median/mode values exclusively from the training fold.

## Using the fitted GESS model

The repository contains a fitted hierarchical model and its preprocessing objects in `GESS_model/`. The wrapper currently resolves model-object paths relative to that directory, so run it with `GESS_model` as the working directory:

```r
setwd("GESS_model")
source("GESS_prediction_source_code.R")

# new_data must follow Dataset_schema.txt exactly.
# Do not include direct identifiers.
predicted_syntax <- GESS_prediction(new_data)
```

### Input requirements

- One row per participant.
- Column names and encodings must match `GESS_model/Dataset_schema.txt`.
- SNP genotypes must use the allele-pair strings and levels documented in the schema.
- Required predictors must be present.
- Rows with missing values in predictors required by either component are excluded by the current wrapper; impute or otherwise resolve missing values before prediction if retaining every row is necessary.

The function returns rounded predicted SYNTAX scores for eligible complete rows. Preserve an external, nonidentifying row key if predictions must be joined back to a source table.

The fitted model is supplied for research transparency only. Its availability does not establish external validity, calibration in new settings, clinical utility, or readiness for deployment.

## Validation and reported quantities

The analysis code includes procedures for:

- pooled out-of-fold AUC-ROC and regression-error estimates;
- fold-level performance summaries;
- participant-level bootstrap confidence intervals;
- paired permutation comparisons of classifier predictions;
- calibration, Brier score, confusion-matrix measures, and exploratory decision-curve analysis for the zero/non-zero classifier;
- complete-pipeline and clinically relevant subgroup analyses;
- feature-selection frequency and exploratory SHAP analyses.

Refer to the accompanying manuscript and supplementary material for the prespecified definitions, final numerical estimates, and interpretation. Results may vary across pipeline components and should not be interpreted as uniformly demonstrating incremental clinical value from SNP information.

## Reproducibility notes

- Random seed 123 is specified in the principal fold-generation and multiple analysis scripts.
- Saved fold-index and model objects are included where available.
- `sessionInfo_complete.txt` documents the environment in which the repository information was generated. If this differs from the historical environment used for a particular analysis, it should not be interpreted as a retrospective lockfile for that earlier run.
- Because the original data are restricted and several scripts require path configuration, this release provides **code and model transparency**, not unrestricted one-command reproduction of the study.

## Citation

If you use or adapt this repository, please cite the associated publication when available:

> *Predicting coronary artery disease severity through genomic profiling and machine learning modelling.* ClinicalTrials.gov identifier: NCT03150680. [Full citation and DOI to be added after publication.]

## License

This repository is distributed under the GNU General Public License version 2; see `LICENSE`.

## Contact and issues

For questions about the code, open a GitHub issue. Requests involving the restricted participant-level dataset must follow the study’s institutional data-access and ethics procedures and should not include participant information in a public issue.
