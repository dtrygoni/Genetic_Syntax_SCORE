# Genetic SYNTAX Score (GESS)

R code and fitted model objects accompanying the study **“Predicting coronary artery disease severity through genomic profiling and machine learning modelling.”**

The GEnetic SYNTAX Score (GESS) framework combines clinical and laboratory measurements with candidate single-nucleotide polymorphisms (SNPs) to predict coronary artery disease presence and severity, expressed as the SYNTAX score.

The repository contains model development, cross-validation, component and complete-pipeline evaluation, missing-data sensitivity analyses, and linear mixed-effects (LME) comparisons of prediction errors.

## Prediction framework

The hierarchical pipeline contains three stages:

1. **Zero/non-zero classification:** predicts a SYNTAX score of 0 versus >0.
2. **Positive-score severity classification:** distinguishes low complexity (<23) from intermediate/high complexity (≥23).
3. **Stratum-specific regression:** uses the low- or intermediate/high-complexity regressor to produce a positive-score estimate. The zero/non-zero classifier determines whether the final prediction is zero or the positive-score estimate.

Two predictor configurations are evaluated:

- **Model A — Clinical:** clinical and laboratory predictors.
- **Model B — Clinical + SNPs:** clinical and laboratory predictors together with SNP-derived predictors and pairwise SNP interactions.

The analyses use predefined outer folds to evaluate both configurations. Component-level performance and complete-pipeline performance are evaluated separately.

## Repository contents

| Path | Contents |
|---|---|
| `Zero_non_zero_component/` | Fold-specific data preparation, feature selection, classifier training, ROC/PR evaluation, and calibration. |
| `Positive_score_component/` | Positive-score data preparation, severity classifier and stratum-specific regressor training, and evaluation. |
| `Whole_model_experiments/` | Combination of the prediction components and complete-pipeline evaluation. |
| `Repeated_Cross_Validation/` | Repeated stratified cross-validation, repetition-specific fold indices and fitted objects, complete-pipeline evaluation, and feature-selection stability. |
| `Sensitivity_Analysis_Imputation/` | Alternative analysis using training-derived median/mode imputation. |
| `LME_experiments/` | Linear mixed-effects comparisons of prediction errors and diagnostic/estimated-marginal-mean plots. |
| `Additional_analysis/` | Decision-curve and ROC plots, missingness summaries, SNP summaries, and subgroup analyses. |
| `GESS_model/` | Fitted prediction objects, input schema, and the `GESS_prediction()` wrapper. |
| `Test_inds.Rda` | Saved primary outer-test fold assignments. |
| `Count_part.Rda` | Saved positive-score prediction results. |
| `Full_model_results.Rda` | Saved complete-pipeline prediction results. |
| `sessionInfo_complete.txt` | Recorded R environment and package information. |

Saved `.Rda` files include fitted models, fold assignments, preprocessing objects, and analysis outputs. The participant-level development dataset is not included.

## Data and input requirements

The development scripts load a restricted file named `Whole_dataset.Rda`, which must create an R object named `dataset`. It contains clinical predictors, SNP predictors, and observed SYNTAX outcomes. The scripts use the original variable names and coding, including `SYNTAX.SCORE`, `Binary.SS`, and `SYNTAX.Tertiles`.

The participant-level clinical and genetic data are subject to consent, ethics, and data-protection restrictions. Reproduction with the original data requires authorized access through the responsible institution.

The fitted-model input variables and genotype coding are described in [`GESS_model/Dataset_schema.txt`](GESS_model/Dataset_schema.txt).

## Software requirements

The project is implemented in R. Consult [`sessionInfo_complete.txt`](sessionInfo_complete.txt) for the recorded software environment. Dependencies vary by script; the main packages include:

```r
packages <- c(
  "Boruta", "caret", "dplyr", "tidyr", "tibble", "ggplot2",
  "randomForest", "xgboost", "PRROC", "precrec", "mltools",
  "pROC", "SNPassoc", "recipes", "MLmetrics", "missForestPredict",
  "reticulate", "rms", "CalibrationCurves", "tidyverse",
  "reshape2", "lme4", "lmerTest", "patchwork", "performance",
  "emmeans", "ggExtra"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) install.packages(missing_packages)
```

Install any additional packages requested by the particular scripts being run.

## Analysis workflow

Run the development analyses from the repository root so that paths beginning with `./` resolve correctly. Configure any script-specific paths and create the output directories referenced by `save()` calls. The saved fold assignments must correspond to the row order of the analytical dataset after the script's initial outcome filtering.

### 1. Zero/non-zero component

Run the following scripts in order:

1. [`1.Cross_validation_construction.R`](Zero_non_zero_component/1.Cross_validation_construction.R): constructs fold-specific training/test datasets and performs SNP processing, univariate screening, complete-case filtering, and Boruta selection.
2. [`2.ML_zero_model_train.R`](Zero_non_zero_component/2.ML_zero_model_train.R): constructs SNP interactions, applies preprocessing, tunes and fits classifiers, and saves fold-specific outputs.
3. [`3.CI_ROC.R`](Zero_non_zero_component/3.CI_ROC.R): evaluates classifier predictions, including pooled ROC/PR performance and a paired permutation comparison of AUCs.
4. [`4.Calibration_with_RMS.R`](Zero_non_zero_component/4.Calibration_with_RMS.R): evaluates calibration.

### 2. Positive-score component

Run:

1. [`1.CV_count_classifier.R`](Positive_score_component/1.CV_count_classifier.R): prepares positive-score training/test datasets and selects predictors.
2. [`2.ML_components_train.R`](Positive_score_component/2.ML_components_train.R): trains the severity classifier and the low- and intermediate/high-complexity regressors.
3. [`3.Risk_Classifier_Evaluation_CI_ROC.R`](Positive_score_component/3.Risk_Classifier_Evaluation_CI_ROC.R): evaluates severity classification.
4. [`4.Low_Risk_Regression_Evaluation.R`](Positive_score_component/4.Low_Risk_Regression_Evaluation.R) and [`5.High_Risk_regression_Evaluation.R`](Positive_score_component/5.High_Risk_regression_Evaluation.R): evaluate the stratum-specific regressors.
5. [`6.Count_part_prediction.R`](Positive_score_component/6.Count_part_prediction.R): combines severity classification and regression predictions and writes `Count_part.Rda`.

### 3. Complete hierarchical model

Run [`Whole_model_experiments/1.SYNTAX_prediction_final.R`](Whole_model_experiments/1.SYNTAX_prediction_final.R) after producing the component models and fold-specific datasets. This script combines the zero/non-zero classifier, severity classifier, and regressors, and writes `Full_model_results.Rda`.

### Fold-specific processing and threshold selection

The development scripts contain the following processing stages:

- HWE testing with FDR-adjusted p-values and minor-allele-frequency filtering.
- LD assessment using squared correlations and reduction of correlated SNP groups using training-set associations.
- Training-derived inheritance-model selection and genotype coding.
- Univariate screening at p < 0.10 and complete-case filtering before Boruta.
- Boruta selection for the zero/non-zero and positive-score components.
- Pairwise categorical SNP interactions among retained SNPs.
- One-hot encoding and min–max scaling estimated from training data.
- Inner four-fold cross-validation for Random Forest hyperparameter tuning.

The revised `Optimal_threshold()` function performs an additional four-fold procedure within the outer training data, holding the selected RF hyperparameters fixed. It pools the held-out probabilities, selects the Youden-index threshold, and returns the largest finite optimal threshold when there is a tie. This calculation runs in memory. The final outer-training classifier supplies the probabilities to which that threshold is applied during test prediction.

## Repeated cross-validation

[`Repeated_Cross_Validation/1.Created_repeated_crossvalidation_inds.R`](Repeated_Cross_Validation/1.Created_repeated_crossvalidation_inds.R) uses `caret::createMultiFolds()` with `k = 10`, `times = 10`, and seed 123.

The remaining numbered scripts prepare data, fit components, evaluate positive-score and complete-pipeline predictions, and summarize feature-selection stability. The directories `Fold1/` through `Fold10/` represent repetitions, each containing its own ten outer folds. Set the repetition variable used by each script before running it, and match any differently capitalized directory strings to the checkout on case-sensitive systems.

## Missing-data sensitivity analysis

`Sensitivity_Analysis_Imputation/` mirrors the component and complete-pipeline workflow. Its imputation helper estimates replacement values from the outer training data and applies them to the corresponding test data:

- median imputation for continuous numeric variables;
- mode imputation for categorical and discrete variables, including SNP predictors.

The helper returns imputed datasets, imputation values, and a summary of the methods used. Consult its arguments for excluded variables and the treatment of discrete numeric variables.

## Linear mixed-effects comparison of prediction errors

[`LME_experiments/LME Performance Evaluation.R`](LME_experiments/LME%20Performance%20Evaluation.R) compares the predictive models using square-root-transformed absolute errors.

### Input table

The script loads `Count_part.Rda` and expects a long-format object named `res_df`. Prepare that object from the prediction results before fitting the models if the loaded file contains only a prediction table named `res`.

| Column | Meaning |
|---|---|
| `ID` | Participant identifier shared by the paired Model A and Model B observations. |
| `Model` | Predictive-model configuration, represented as a factor. |
| `Risk` | Observed SYNTAX complexity stratum, represented as a factor. |
| `Fold` | Outer fold; represent it as a factor for categorical fold adjustment. |
| `AE` | Absolute error, calculated as `abs(observed_score - predicted_score)`. |

Each participant contributes one error for each predictive model being compared. Preserve the same participant identifiers and fold assignments across the paired observations.

### Model specification

The reduced and interaction models are fitted by maximum likelihood:

```r
model_reduced <- lmerTest::lmer(
  sqrt(AE) ~ Model + Risk + Fold + (1 | ID),
  data = res_df, REML = FALSE
)

model_full <- lmerTest::lmer(
  sqrt(AE) ~ Model * Risk + Fold + (1 | ID),
  data = res_df, REML = FALSE
)

anova(model_reduced, model_full)
```

`Model`, `Risk`, and `Fold` enter as fixed effects. The participant random intercept accounts for the paired errors. The Model × Risk interaction assesses whether the difference between predictive configurations varies across complexity strata.

The script compares the models using a likelihood-ratio test and then fits the interaction model using REML for final estimation. It reports coefficient summaries and fixed-effect tests, and produces residual-versus-fitted, residual Q–Q, scale–location, and random-intercept Q–Q plots.

Estimated marginal means are calculated using `emmeans(model_interaction, ~ Model | Risk)`, followed by within-stratum pairwise comparisons with the requested Tukey adjustment. These comparisons are on the square-root AE scale. The script also uses `regrid(..., transform = "response")` to display back-transformed estimated marginal means and confidence intervals; these should be distinguished from raw arithmetic mean absolute errors.

The supplied plotting labels cover the two positive-score strata. For complete-pipeline analysis, use the corresponding long-format complete-pipeline errors and adapt the risk labels to include zero, low, and intermediate/high SYNTAX scores.

## Additional analyses

`Additional_analysis/` contains scripts for decision-curve analysis, missing-predictor summaries, SNP QC summaries, SNP table construction, ROC plotting, and subgroup evaluation of the complete pipeline. Run these after generating the prediction outputs they load.

## Using the fitted GESS model

The prediction wrapper resolves fitted-object paths relative to `GESS_model/`:

```r
setwd("GESS_model")
source("GESS_prediction_source_code.R")

# new_data must follow Dataset_schema.txt.
predicted_syntax <- GESS_prediction(new_data)
```

Supply one row per participant, with the exact predictor names and genotype encodings described in the schema. The wrapper applies the saved component objects and returns predictions for eligible rows. Preserve an external row key when joining predictions back to a source dataset.

## Reproducibility and interpretation

The repository supplies source code, seeds, fold indices, and saved model objects. Regenerating results requires the authorized development dataset, the expected R environment, and the intermediate files produced by the preceding scripts. Saved outputs should be regenerated after changes to preprocessing or model development.

Performance quantities are defined in the individual evaluation scripts. These include classification discrimination, regression errors, pooled and fold-level summaries, resampling-based comparisons, calibration analyses, and LME comparisons of paired errors. Consult the manuscript and supplementary material for final numerical results and the study's analysis definitions.

GESS is a research model developed with internal validation. External and prospective validation are required before clinical application.

## Citation

**Predicting coronary artery disease severity through genomic profiling and machine learning modelling.** ClinicalTrials.gov identifier: NCT03150680. Full publication details and DOI will be added when available.

## License

Distributed under the GNU General Public License version 2. See [`LICENSE`](LICENSE).

## Contact

Use GitHub issues for code questions. Requests for restricted participant-level data should follow the responsible institution's data-access procedures.
