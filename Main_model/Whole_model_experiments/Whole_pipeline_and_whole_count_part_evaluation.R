load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part.Rda")



######### Fold-level performance ###########
folds<-c()
MAE_clinical<-c()
MAE_snp<-c()
test_all<-res$Target
folds_arr<-res$Fold
clinical_preds<-res$Clinical.Predicted
snp_preds<-res$SNP.Predicted
for(fold in 1:10){
  ae_clinical<-abs(test_all[which(folds_arr==fold)]-clinical_preds[which(folds_arr==fold)])
  mae_clinical<-median(ae_clinical)
  
  MAE_clinical<-c(MAE_clinical,mae_clinical)
  
  
  ae_snp<-abs(test_all[which(folds_arr==fold)]-snp_preds[which(folds_arr==fold)])
  mae_snp<-median(ae_snp)
  
  MAE_snp<-c(MAE_snp,mae_snp)
  
  folds<-c(folds,fold)
}
res<-data.frame(Fold=folds,MAE.clinical=MAE_clinical,MAE.SNP=MAE_snp)
quantile(res$MAE.SNP) 






y_true<-res$Target
pred_model_A<-res$Clinical.Predicted
pred_model_B<-res$SNP.Predicted

n_boot <- 5000
boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)

complete_index <- which(
  complete.cases(y_true, pred_model_A, pred_model_B)
)

for (i in seq_len(n_boot)) {
  
  boot_index <- sample(
    complete_index,
    length(complete_index),
    replace = TRUE
  )
  
  boot_medae_A[i] <- median(
    abs(y_true[boot_index] - pred_model_A[boot_index])
  )
  
  boot_medae_B[i] <- median(
    abs(y_true[boot_index] - pred_model_B[boot_index])
  )
  
  boot_delta_medae[i] <-
    boot_medae_B[i] - boot_medae_A[i]
}

ci_medae_A <- quantile(
  boot_medae_A,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_medae_B <- quantile(
  boot_medae_B,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_delta_medae <- quantile(
  boot_delta_medae,
  c(0.025, 0.975),
  na.rm = TRUE
)


median(boot_medae_B)
ci_delta_medae



load("C:/Users/dtryg/Desktop/GESS_experiments_code/Full_model.Rda")

######### Fold-level performance ###########
folds<-c()
MAE_clinical<-c()
MAE_snp<-c()
test_all<-res$Target
folds_arr<-res$Fold
clinical_preds<-res$Clinical.Predicted
snp_preds<-res$SNP.Predicted
for(fold in 1:10){
  ae_clinical<-abs(test_all[which(folds_arr==fold)]-clinical_preds[which(folds_arr==fold)])
  mae_clinical<-median(ae_clinical)
  
  MAE_clinical<-c(MAE_clinical,mae_clinical)
  
  
  ae_snp<-abs(test_all[which(folds_arr==fold)]-snp_preds[which(folds_arr==fold)])
  mae_snp<-median(ae_snp)
  
  MAE_snp<-c(MAE_snp,mae_snp)
  
  folds<-c(folds,fold)
}
res<-data.frame(Fold=folds,MAE.clinical=MAE_clinical,MAE.SNP=MAE_snp)
quantile(res$MAE.clinical) 


y_true<-res$Target
pred_model_A<-res$Clinical.Predicted
pred_model_B<-res$SNP.Predicted

n_boot <- 5000
boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)

complete_index <- which(
  complete.cases(y_true, pred_model_A, pred_model_B)
)

for (i in seq_len(n_boot)) {
  
  boot_index <- sample(
    complete_index,
    length(complete_index),
    replace = TRUE
  )
  
  boot_medae_A[i] <- median(
    abs(y_true[boot_index] - pred_model_A[boot_index])
  )
  
  boot_medae_B[i] <- median(
    abs(y_true[boot_index] - pred_model_B[boot_index])
  )
  
  boot_delta_medae[i] <-
    boot_medae_B[i] - boot_medae_A[i]
}

ci_medae_A <- quantile(
  boot_medae_A,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_medae_B <- quantile(
  boot_medae_B,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_delta_medae <- quantile(
  boot_delta_medae,
  c(0.025, 0.975),
  na.rm = TRUE
)


median(boot_medae_B)
ci_delta_medae

mean(abs(res$Target-res$Clinical.Predicted))
mean(abs(res$Target-res$SNP.Predicted))

sqrt(mean((res$Target - res$Clinical.Predicted)^2, na.rm = TRUE))
sqrt(mean((res$Target - res$SNP.Predicted)^2, na.rm = TRUE))




res_complete <- res[
  complete.cases(
    res$Target,
    res$Clinical.Predicted,
    res$SNP.Predicted
  ),
]

# Observed performance
mae_clinical <- mean(
  abs(res_complete$Target - res_complete$Clinical.Predicted)
)

mae_snp <- mean(
  abs(res_complete$Target - res_complete$SNP.Predicted)
)

rmse_clinical <- sqrt(
  mean(
    (res_complete$Target - res_complete$Clinical.Predicted)^2
  )
)

rmse_snp <- sqrt(
  mean(
    (res_complete$Target - res_complete$SNP.Predicted)^2
  )
)

# Paired differences: Model B (SNP) minus Model A (clinical)
delta_mae <- mae_snp - mae_clinical
delta_rmse <- rmse_snp - rmse_clinical




set.seed(123)

n_boot <- 5000
n_participants <- nrow(res_complete)

boot_results <- matrix(
  NA_real_,
  nrow = n_boot,
  ncol = 6
)

colnames(boot_results) <- c(
  "MAE_Clinical",
  "MAE_SNP",
  "Delta_MAE",
  "RMSE_Clinical",
  "RMSE_SNP",
  "Delta_RMSE"
)

for (i in seq_len(n_boot)) {
  
  # Resample participants while retaining paired predictions
  boot_index <- sample(
    seq_len(n_participants),
    size = n_participants,
    replace = TRUE
  )
  
  boot_data <- res_complete[boot_index, ]
  
  # MAE
  boot_mae_clinical <- mean(
    abs(
      boot_data$Target -
        boot_data$Clinical.Predicted
    )
  )
  
  boot_mae_snp <- mean(
    abs(
      boot_data$Target -
        boot_data$SNP.Predicted
    )
  )
  
  # RMSE
  boot_rmse_clinical <- sqrt(
    mean(
      (
        boot_data$Target -
          boot_data$Clinical.Predicted
      )^2
    )
  )
  
  boot_rmse_snp <- sqrt(
    mean(
      (
        boot_data$Target -
          boot_data$SNP.Predicted
      )^2
    )
  )
  
  boot_results[i, ] <- c(
    boot_mae_clinical,
    boot_mae_snp,
    boot_mae_snp - boot_mae_clinical,
    boot_rmse_clinical,
    boot_rmse_snp,
    boot_rmse_snp - boot_rmse_clinical
  )
}




bootstrap_ci <- apply(
  boot_results,
  2,
  quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
)

bootstrap_ci



results_table <- data.frame(
  Metric = c(
    "MAE: Model A (Clinical)",
    "MAE: Model B (Clinical+SNPs)",
    "Delta MAE: Model B - Model A",
    "RMSE: Model A (Clinical)",
    "RMSE: Model B (Clinical+SNPs)",
    "Delta RMSE: Model B - Model A"
  ),
  Estimate = c(
    mae_clinical,
    mae_snp,
    delta_mae,
    rmse_clinical,
    rmse_snp,
    delta_rmse
  ),
  CI_lower = bootstrap_ci[1, ],
  CI_upper = bootstrap_ci[2, ]
)

results_table[, 2:4] <- round(
  results_table[, 2:4],
  digits = 3
)

results_table







################## Classification accuracies around 23 ###############
library(caret)

predict_clinical<-ifelse(res$Clinical.Predicted<23,0,1)
predict_snp<-ifelse(res$SNP.Predicted<23,0,1)
target_clinical<-ifelse(res$Target<23,0,1)

caret::confusionMatrix(factor(predict_snp),factor(target_clinical))







############### uncertainty intervals ###########
# Keep identical participants for both models
boot_data <- res[
  complete.cases(
    res$Target,
    res$Clinical.Predicted,
    res$SNP.Predicted
  ),
]

# Function that calculates all metrics in one dataset
calculate_metrics <- function(data) {
  
  observed <- data$Target
  clinical <- data$Clinical.Predicted
  snp <- data$SNP.Predicted
  
  # Continuous prediction metrics
  mae_clinical <- mean(abs(observed - clinical))
  mae_snp <- mean(abs(observed - snp))
  
  rmse_clinical <- sqrt(mean((observed - clinical)^2))
  rmse_snp <- sqrt(mean((observed - snp)^2))
  
  # Classification around the threshold of 23
  observed_high <- observed >= 23
  clinical_high <- clinical >= 23
  snp_high <- snp >= 23
  
  accuracy_clinical <- mean(clinical_high == observed_high)
  accuracy_snp <- mean(snp_high == observed_high)
  
  sensitivity_clinical <- mean(
    clinical_high[observed_high]
  )
  
  sensitivity_snp <- mean(
    snp_high[observed_high]
  )
  
  specificity_clinical <- mean(
    !clinical_high[!observed_high]
  )
  
  specificity_snp <- mean(
    !snp_high[!observed_high]
  )
  
  c(
    MAE_Clinical = mae_clinical,
    MAE_SNP = mae_snp,
    Delta_MAE = mae_snp - mae_clinical,
    
    RMSE_Clinical = rmse_clinical,
    RMSE_SNP = rmse_snp,
    Delta_RMSE = rmse_snp - rmse_clinical,
    
    Accuracy_Clinical = accuracy_clinical,
    Accuracy_SNP = accuracy_snp,
    Delta_Accuracy = accuracy_snp - accuracy_clinical,
    
    Sensitivity_Clinical = sensitivity_clinical,
    Sensitivity_SNP = sensitivity_snp,
    
    Specificity_Clinical = specificity_clinical,
    Specificity_SNP = specificity_snp
  )
}

# Observed estimates
observed_estimates <- calculate_metrics(boot_data)



set.seed(123)

n_boot <- 5000
n_participants <- nrow(boot_data)

bootstrap_estimates <- replicate(
  n_boot,
  {
    # Resample participant rows with replacement
    boot_index <- sample(
      seq_len(n_participants),
      size = n_participants,
      replace = TRUE
    )
    
    calculate_metrics(
      boot_data[boot_index, ]
    )
  }
)
bootstrap_ci <- t(
  apply(
    bootstrap_estimates,
    1,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
)

colnames(bootstrap_ci) <- c(
  "CI_lower",
  "CI_upper"
)

bootstrap_results <- data.frame(
  Metric = names(observed_estimates),
  Estimate = as.numeric(observed_estimates),
  CI_lower = bootstrap_ci[, "CI_lower"],
  CI_upper = bootstrap_ci[, "CI_upper"],
  row.names = NULL
)

bootstrap_results[, c(
  "Estimate",
  "CI_lower",
  "CI_upper"
)] <- round(
  bootstrap_results[, c(
    "Estimate",
    "CI_lower",
    "CI_upper"
  )],
  3
)

bootstrap_results


percentage_metrics <- grepl(
  "Accuracy|Sensitivity|Specificity",
  bootstrap_results$Metric
)

bootstrap_results_percent <- bootstrap_results

bootstrap_results_percent[
  percentage_metrics,
  c("Estimate", "CI_lower", "CI_upper")
] <- 100 * bootstrap_results_percent[
  percentage_metrics,
  c("Estimate", "CI_lower", "CI_upper")
]

bootstrap_results_percent


############ Calibration slope #############
# Retain identical complete cases
calibration_data <- res[
  complete.cases(
    res$Target,
    res$Clinical.Predicted,
    res$SNP.Predicted
  ),
]

# Model A: clinical-only
calibration_clinical <- lm(
  Target ~ Clinical.Predicted,
  data = calibration_data
)

# Model B: clinical+SNPs
calibration_snp <- lm(
  Target ~ SNP.Predicted,
  data = calibration_data
)

summary(calibration_clinical)
summary(calibration_snp)


clinical_coefficients <- coef(
  calibration_clinical
)

snp_coefficients <- coef(
  calibration_snp
)

clinical_ci <- confint(
  calibration_clinical,
  level = 0.95
)

snp_ci <- confint(
  calibration_snp,
  level = 0.95
)

calibration_results <- data.frame(
  Model = rep(
    c(
      "Model A: Clinical",
      "Model B: Clinical+SNPs"
    ),
    each = 2
  ),
  Parameter = rep(
    c(
      "Calibration intercept",
      "Calibration slope"
    ),
    times = 2
  ),
  Estimate = c(
    clinical_coefficients[1],
    clinical_coefficients[2],
    snp_coefficients[1],
    snp_coefficients[2]
  ),
  CI_lower = c(
    clinical_ci[1, 1],
    clinical_ci[2, 1],
    snp_ci[1, 1],
    snp_ci[2, 1]
  ),
  CI_upper = c(
    clinical_ci[1, 2],
    clinical_ci[2, 2],
    snp_ci[1, 2],
    snp_ci[2, 2]
  )
)

calibration_results[, 3:5] <- round(
  calibration_results[, 3:5],
  3
)

calibration_results


calculate_calibration <- function(data) {
  
  fit_clinical <- lm(
    Target ~ Clinical.Predicted,
    data = data
  )
  
  fit_snp <- lm(
    Target ~ SNP.Predicted,
    data = data
  )
  
  c(
    Intercept_Clinical = coef(fit_clinical)[1],
    Slope_Clinical = coef(fit_clinical)[2],
    Intercept_SNP = coef(fit_snp)[1],
    Slope_SNP = coef(fit_snp)[2]
  )
}

# Observed calibration estimates
observed_calibration <- calculate_calibration(
  calibration_data
)

set.seed(123)

n_boot <- 5000
n_participants <- nrow(calibration_data)

bootstrap_calibration <- replicate(
  n_boot,
  {
    boot_index <- sample(
      seq_len(n_participants),
      size = n_participants,
      replace = TRUE
    )
    
    calculate_calibration(
      calibration_data[boot_index, ]
    )
  }
)

# Percentile 95% bootstrap CIs
calibration_ci <- t(
  apply(
    bootstrap_calibration,
    1,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
)

calibration_bootstrap_results <- data.frame(
  Parameter = names(observed_calibration),
  Estimate = as.numeric(observed_calibration),
  CI_lower = calibration_ci[, 1],
  CI_upper = calibration_ci[, 2],
  row.names = NULL
)

calibration_bootstrap_results[, 2:4] <- round(
  calibration_bootstrap_results[, 2:4],
  3
)

calibration_bootstrap_results


cal_model <- lm(res$Target ~ res$Clinical.Predicted)
summary(cal_model)
# 2. Visual Calibration Plot (Observed vs. Predicted)
plot(res$Clinical.Predicted, res$Target, 
     xlab = "Predicted Values", ylab = "Observed Values (test_all)",
     main = "Regression Calibration Plot", pch = 16, col = rgb(0,0,1,0.5))
# Add perfect calibration line (LOESS/Smooth trend)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2) # Perfect line
lines(lowess(res$Clinical.Predicted, res$Target), col = "blue", lwd = 2) # Actual model trend

legend("topleft", legend = c("Perfect Calibration", "Model Trend"),
       col = c("red", "blue"), lty = c(2, 1), lwd = 2)
