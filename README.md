# Genetic_Syntax_SCORE
Machine learning framework integrating clinical, laboratory, and genomic (SNP) data to predict the presence and severity of coronary artery disease (CAD) via the GEnetic SYNTAX Score (GESS) model. Enhances CAD risk stratification and supports precision cardiovascular medicine.

Predicting coronary artery disease severity through genomic profiling and machine learning modelling

This repository contains the implementation of the GEnetic SYNTAX Score (GESS) framework — a machine learning approach that integrates clinical, laboratory, and genomic (SNP) data to predict both the presence and severity of obstructive coronary artery disease (CAD) as expressed by the SYNTAX score.

The code supports the full modeling pipeline, including data preprocessing, feature selection, cross-validation experiments, and the final trained model for practical use.

# Repository Structure

Framework_modelling – End-to-end ML workflow integrating preprocessing, feature selection, model training, and performance evaluation.
Zero_part_modelling – Classification component predicting the presence of obstructive CAD (SYNTAX score = 0 vs >0).
Count_part_modelling – Regression component predicting the severity of CAD (numeric SYNTAX score) among obstructive cases.
GESS_model – Final fitted model(s) and wrapper functions for direct prediction on new datasets.
Utils – Utility scripts for data handling, feature encoding, performance metrics, and visualization (ROC, PR, SHAP).
results – Folder to store output figures, metrics, and validation results (optional).

# Requirements

The project is implemented in R (version 4.2 or later).

Required packages:
caret, pROC, Boruta, lme4, lmerTest, emmeans, SNPassoc, kernelshap, shapviz, tidyverse, data.table

Install them in R with:
pkgs <- c("caret","pROC","Boruta","lme4","lmerTest","emmeans","SNPassoc","kernelshap","shapviz","tidyverse","data.table")
inst <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst, dependencies = TRUE)

# Data

The workflow expects two data sources:

Clinical, laboratory, and demographic data

Genomic data (SNP genotypes) after quality control (HWE, MAF, LD filtering)

Create a local folder named data (not tracked in GitHub) and place your input files there.
Paths can be specified at the beginning of each script.

# How to Run the Experiments

Step 1. Preprocessing and setup
Run the main scripts in Framework_modelling to clean, normalize, and encode data for both clinical and genomic predictors.

Step 2. Train and cross-validate models
Run:

Zero_part_modelling for binary classification of CAD presence

Count_part_modelling for regression on SYNTAX score (severity)

Each pipeline:

Uses stratified 10-fold cross-validation

Computes AUC/ROC, PR curves, and MAE

Applies Boruta feature selection and permutation testing for model comparison

Step 3. Review results
All output figures, metrics, and SHAP summaries are saved to the results folder.

# Using the Final Model

You can use the trained model in the GESS_model folder to make predictions on new patients.

Example:

source("Utils/helper_functions.R")
source("GESS_model/GESS_predict.R")

Examine carefully the variable names and the modelling of the SNPs selected. The variables named of each model can be reviewed with the function predictors(zero_part) or predictors(count_part).

The prediction function:

Classifies patients for presence of obstructive CAD (zero-part model)

If CAD is likely, estimates severity (SYNTAX score) using the count-part regressor

Results Summary

Incorporating genetic markers (SNPs) improved prediction of CAD presence (AUC = 0.760 vs 0.731 for clinical-only model).

Performance for SYNTAX score prediction varied by risk level:
Low risk: genetic model improved accuracy (lower MAE)
Intermediate/high risk: clinical model remained more stable.

SHAP-based explainability identified both clinical and genetic predictors influencing model decisions.

# Citation

If you use or adapt this code, please cite the corresponding publication: XXX (ClinicalTrials.gov Identifier: NCT03150680)
