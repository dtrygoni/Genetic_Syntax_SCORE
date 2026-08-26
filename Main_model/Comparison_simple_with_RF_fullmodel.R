load("C:/Users/dtryg/Desktop/GESS_experiments_code/Full_model.Rda")


res$SNP.Predicted.Simple<-pred_model_B
res$Clinical.Predicted.Simple<-pred_model_A

save(res,file="C:/Users/dtryg/Desktop/GESS_experiments_code/Full_model_with_SimpleClassification.Rda")
y_true<-res$Target
pred_model_A<-res$SNP.Predicted.Simple
pred_model_B<-res$SNP.Predicted

n_boot <- 5000
boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)


for (i in seq_len(n_boot)) {
  
  boot_index <- sample(
    complete_index,
    length(complete_index),
    replace = TRUE
  )
  
  boot_medae_A[i] <- median(
    abs(y_true[boot_index] - pred_model_A[boot_index]),na.rm=TRUE
  )
  
  boot_medae_B[i] <- median(
    abs(y_true[boot_index] - pred_model_B[boot_index]),na.rm=TRUE
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

################# genai ###############
# ============================================================
# Paired comparison of median absolute error
# Model B minus Model A
# ============================================================

# Observed outcomes and out-of-fold predictions
y_true <- res$Target
pred_model_A <- res$Clinical.Predicted.Simple
pred_model_B <- res$Clinical.Predicted

# ------------------------------------------------------------
# 1. Retain identical complete cases for both models
# ------------------------------------------------------------

complete_index <- which(
  complete.cases(
    y_true,
    pred_model_A,
    pred_model_B
  )
)

y_complete <- y_true[complete_index]
pred_A_complete <- pred_model_A[complete_index]
pred_B_complete <- pred_model_B[complete_index]

n_complete <- length(y_complete)

# Check the final sample size
n_complete

# ------------------------------------------------------------
# 2. Calculate participant-level absolute errors
# ------------------------------------------------------------

absolute_error_A <- abs(
  y_complete - pred_A_complete
)

absolute_error_B <- abs(
  y_complete - pred_B_complete
)

# ------------------------------------------------------------
# 3. Calculate observed MedAE estimates
# ------------------------------------------------------------

medae_A_observed <- median(
  absolute_error_A
)

medae_B_observed <- median(
  absolute_error_B
)

# Signed difference: Model B minus Model A
delta_medae_observed <-
  medae_B_observed -
  medae_A_observed

medae_A_observed
medae_B_observed
delta_medae_observed

# ------------------------------------------------------------
# 4. Paired participant-level bootstrap
# ------------------------------------------------------------

set.seed(123)

n_boot <- 5000

boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  
  # Resample participants with replacement
  boot_index <- sample(
    seq_len(n_complete),
    size = n_complete,
    replace = TRUE
  )
  
  # MedAE for each model in the same bootstrap sample
  boot_medae_A[i] <- median(
    absolute_error_A[boot_index]
  )
  
  boot_medae_B[i] <- median(
    absolute_error_B[boot_index]
  )
  
  # Paired difference: Model B minus Model A
  boot_delta_medae[i] <-
    boot_medae_B[i] -
    boot_medae_A[i]
}

# ------------------------------------------------------------
# 5. Percentile bootstrap 95% confidence intervals
# ------------------------------------------------------------

ci_medae_A <- quantile(
  boot_medae_A,
  probs = c(0.025, 0.975),
  na.rm = TRUE,
  names = FALSE
)

ci_medae_B <- quantile(
  boot_medae_B,
  probs = c(0.025, 0.975),
  na.rm = TRUE,
  names = FALSE
)

ci_delta_medae <- quantile(
  boot_delta_medae,
  probs = c(0.025, 0.975),
  na.rm = TRUE,
  names = FALSE
)

# ------------------------------------------------------------
# 6. Paired permutation test
# ------------------------------------------------------------

set.seed(123)

n_perm <- 5000
delta_medae_perm <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  
  # Randomly exchange the two errors within each participant
  swap <- rbinom(
    n = n_complete,
    size = 1,
    prob = 0.5
  ) == 1
  
  perm_error_A <- ifelse(
    swap,
    absolute_error_B,
    absolute_error_A
  )
  
  perm_error_B <- ifelse(
    swap,
    absolute_error_A,
    absolute_error_B
  )
  
  # Permuted difference: Model B minus Model A
  delta_medae_perm[i] <-
    median(perm_error_B) -
    median(perm_error_A)
}

# Two-sided permutation p-value
permutation_p_value <- (
  1 +
    sum(
      abs(delta_medae_perm) >=
        abs(delta_medae_observed)
    )
) / (
  n_perm + 1
)

# ------------------------------------------------------------
# 7. Compile the results
# ------------------------------------------------------------

results_medae <- data.frame(
  N = n_complete,
  
  MedAE_Model_A = medae_A_observed,
  MedAE_Model_A_CI_lower = ci_medae_A[1],
  MedAE_Model_A_CI_upper = ci_medae_A[2],
  
  MedAE_Model_B = medae_B_observed,
  MedAE_Model_B_CI_lower = ci_medae_B[1],
  MedAE_Model_B_CI_upper = ci_medae_B[2],
  
  Difference_B_minus_A = delta_medae_observed,
  Difference_CI_lower = ci_delta_medae[1],
  Difference_CI_upper = ci_delta_medae[2],
  
  Test_statistic = abs(delta_medae_observed),
  Permutation_p_value = permutation_p_value
)

results_medae <- round(
  results_medae,
  digits = 4
)

results_medae

# ------------------------------------------------------------
# 8. Publication-ready output
# ------------------------------------------------------------

cat(
  sprintf(
    paste0(
      "Model A MedAE = %.2f ",
      "(95%% CI %.2f to %.2f); ",
      "Model B MedAE = %.2f ",
      "(95%% CI %.2f to %.2f); ",
      "paired difference, Model B minus Model A = %.2f ",
      "(95%% CI %.2f to %.2f); ",
      "paired permutation test statistic = %.2f; ",
      "p = %.4f; n = %d.\n"
    ),
    
    medae_A_observed,
    ci_medae_A[1],
    ci_medae_A[2],
    
    medae_B_observed,
    ci_medae_B[1],
    ci_medae_B[2],
    
    delta_medae_observed,
    ci_delta_medae[1],
    ci_delta_medae[2],
    
    abs(delta_medae_observed),
    permutation_p_value,
    n_complete
  )
)
