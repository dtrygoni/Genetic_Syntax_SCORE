rm(list=ls())
######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Models/Models_RF_AUC_None.Rda')



load('C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Test_clinical_zero.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_clinical_zero.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Test_snp_zero.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_snp_zero.Rda')
### Prepare clinical sets
Preprocessing_fnc<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',factor_col_opt='One-hot',numerical_col_opt='None'){
  ########### Factor preprocess
  if(factor_col_opt=='One-hot'){
    dmy <- dummyVars(" ~ .", data = clinical_train)
    trsf_train <- data.frame(predict(dmy, newdata = clinical_train))
    trsf_test<- data.frame(predict(dmy, newdata = clinical_test))
    
    
    dmy_snp <- dummyVars(" ~ .", data = snp_train)
    trsf_trainsnp <- data.frame(predict(dmy_snp, newdata = snp_train))
    trsf_testsnp<- data.frame(predict(dmy_snp, newdata = snp_test))
    
  }else if(factor_col_opt=="L-1_encoding"){
    rec <- recipe(~ ., data = clinical_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_train <- data.frame(bake(prep_rec, new_data = clinical_train))
    trsf_train<-preprocess_baked_structure(trsf_train)
    trsf_test <- data.frame(bake(prep_rec, new_data = clinical_test))
    trsf_test<-preprocess_baked_structure(trsf_test)
    
    
    
    rec <- recipe(~ ., data = snp_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_trainsnp <- data.frame(bake(prep_rec, new_data = snp_train))
    trsf_trainsnp<-preprocess_baked_structure(trsf_trainsnp)
    
    trsf_testsnp <- data.frame(bake(prep_rec, new_data = snp_test))
    trsf_testsnp<-preprocess_baked_structure(trsf_testsnp)
    
    
  }
  
  ########## Numerical preprocess / Scale
  if(numerical_col_opt=='None'){
    clinical_train<-trsf_train
    clinical_test<-trsf_test
    
    snp_train<-trsf_trainsnp
    snp_test<-trsf_testsnp
    scalerClinical_obj<-NULL
    scalerSNP_obj<-NULL
    
  }else if(numerical_col_opt=='MinMax'){
    pre_proc <- preProcess(trsf_train, method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp, method = c("range"))
    # Transform the dataset
    snp_train <- predict(pre_procsnp, newdata = trsf_trainsnp)
    snp_test <- predict(pre_procsnp, newdata = trsf_testsnp)
    
    scalerSNP_obj<-pre_procsnp
    
  }
  
  
  return(list(`ClinicalTrain`=clinical_train,`ClinicalTest`=clinical_test,`SNPtrain`=snp_train,`SNPtest`=snp_test,`ClinicalScaler`=scalerClinical_obj,
              `ClinicalSNPScaler`=scalerSNP_obj))
}
preprocess_baked_structure<-function(data){
  df_clean <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x)) as.numeric(as.character(x))  # convert factors to numeric (via character to keep level values)
    else if (is.matrix(x)) as.numeric(x)           # flatten any matrix columns
    else x
  }))
  
  # Remove all attributes (like "assign", "contrasts")
  df_clean <- as.data.frame(lapply(df_clean, function(x) {
    attributes(x) <- NULL
    return(x)
  }))
  return(df_clean)
  
}

factor_col_opt='One-hot' #One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
train_clinical_ls_processed<-list()
test_clinical_ls_processed<-list()
for(i in 1:10){
  res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[i]],test_clinical_ls[[i]],train_snp_ls[[i]],test_snp_ls[[i]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  train_set<-res_preprocessing$ClinicalTrain
  test_set<-res_preprocessing$ClinicalTest
  
  train_clinical_ls_processed[[i]]<-train_set
  test_clinical_ls_processed[[i]]<-test_set
}


### Prepare snps sets

interaction_analysis <-function(train,test){
  cols_snp<-which(grepl("rs",colnames(train)))
  ## Interaction
  if(length(cols_snp)>=2){
    
    snp_cols<-colnames(train)[cols_snp]
    
    for(snp1 in snp_cols){
      rest_cols<-snp_cols[!snp_cols %in% snp1]
      
      for(snp2 in rest_cols){
        interaction_string<-paste0(snp1,'_',snp2)
        
        interaction_col<-interaction(train[,snp1],train[,snp2])
        interaction_col_test<-interaction(test[,snp1],test[,snp2])
        
        train[[interaction_string]]<-interaction_col
        test[[interaction_string]]<-interaction_col_test
        
        
        
        snp_cols<-snp_cols[!snp_cols %in% snp1]
        
      }
    }
    
  }else{
    return(list(`Train`=train,`Test`=test))
  }
  
  
  
  return(list(`Train`=train,`Test`=test))
  
}
train_snp_ls_processed<-list()
test_snp_ls_processed<-list()
for(i in 1:10){
  res_interaction<-interaction_analysis(train_snp_ls[[i]],test_snp_ls[[i]])
  snp_train<-res_interaction$Train
  snp_test<-res_interaction$Test
  res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[i]],test_clinical_ls[[i]],snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  
  train_snp_ls_processed[[i]]<-res_preprocessing$SNPtrain
  test_snp_ls_processed[[i]]<-res_preprocessing$SNPtest
}

test_all<-c()
probs_clinical<-c()
probs_snp<-c()
preds_clinical<-c()
preds_snp<-c()
threshold_clinical<-c()
threshold_snp<-c()
folds_arr<-c()
for(jj in 1:10){
  #Clinical
  test_clinical <-test_clinical_ls_processed[[jj]]
  
  model_clinical<-final_models[[paste0('Clinical_RF_',jj)]]
  
  train_clinical<-train_clinical_ls_processed[[jj]]
  probs<-predict(model_clinical,train_clinical,type="prob")[,2]
  roc_curve<-pROC::roc(train_clinical$target,probs)
  best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
  threshold<-as.numeric(best_coords$threshold)
  
  probs_cl<-predict(model_clinical,test_clinical,type="prob")[,2]
  preds_cl<-as.factor(ifelse(probs_cl>=threshold,1,0))
  
  preds_clinical<-c(preds_clinical,preds_cl)
  probs_clinical<-c(probs_clinical,probs_cl)
  threshold_clinical<-c(threshold_clinical,threshold)
  #SNP
  test_snp <-test_snp_ls_processed[[jj]] 
  model_snp<-final_models[[paste0('SNP_RF_',jj)]]
  
  train_snp<-train_snp_ls_processed[[jj]]
  probs<-predict(model_snp,train_snp,type="prob")[,2]
  roc_curve<-pROC::roc(train_snp$target,probs)
  best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
  threshold<-as.numeric(best_coords$threshold)
  threshold_snp<-c(threshold_snp,threshold)
  
  probs_snp1<-predict(model_snp,test_snp,type="prob")[,2]
  preds_snp1<-as.factor(ifelse(probs_snp1>=threshold,1,0))
  
  preds_snp<-c(preds_snp,preds_snp1)
  probs_snp<-c(probs_snp,probs_snp1)
  
  
  test_all<-c(test_all,test_clinical$target)
  folds_arr<-c(folds_arr,rep(jj,length(test_clinical$target)))
}

#Metrics over all folds
preds_clinical<-as.factor(ifelse(preds_clinical==2,1,0))
preds_snp<-as.factor(ifelse(preds_snp==2,1,0))
clinical <- caret::confusionMatrix(preds_clinical,as.factor(test_all),positive='1')
snp<-caret::confusionMatrix(preds_snp,as.factor(test_all),positive='1')







############## Comment 48 ####################
evaluation_data <- data.frame(
  outcome = as.numeric(as.character(test_all)),
  prob_clinical = as.numeric(probs_clinical),
  prob_snp = as.numeric(probs_snp),
  Predicted.SNP=preds_snp,
  Predicted.Clinical=preds_clinical
)

library(DescTools)
BrierScore(evaluation_data$outcome,evaluation_data$prob_clinical)
BrierScore(evaluation_data$outcome,evaluation_data$prob_snp)


####################### Calibration ##############################
clip_probability <- function(p, epsilon = 1e-6) {
  pmin(
    pmax(p, epsilon),
    1 - epsilon
  )
}

evaluation_data$lp_clinical <- qlogis(
  clip_probability(
    evaluation_data$prob_clinical
  )
)

evaluation_data$lp_snp <- qlogis(
  clip_probability(
    evaluation_data$prob_snp
  )
)

calibration_clinical <- glm(
  outcome ~ lp_clinical,
  family = binomial,
  data = evaluation_data
)

calibration_snp <- glm(
  outcome ~ lp_snp,
  family = binomial,
  data = evaluation_data
)

summary(calibration_clinical)
summary(calibration_snp)
extract_calibration <- function(model, model_name) {
  
  estimates <- coef(model)
  standard_errors <- sqrt(
    diag(vcov(model))
  )
  
  data.frame(
    Model = model_name,
    Parameter = c(
      "Calibration intercept",
      "Calibration slope"
    ),
    Estimate = estimates,
    CI_lower = estimates -
      1.96 * standard_errors,
    CI_upper = estimates +
      1.96 * standard_errors,
    row.names = NULL
  )
}

calibration_results <- rbind(
  extract_calibration(
    calibration_clinical,
    "Model A: Clinical"
  ),
  extract_calibration(
    calibration_snp,
    "Model B: Clinical+SNPs"
  )
)

calibration_results
############################### Calibration plot ######################
library(ggplot2)

prepare_calibration_data <- function(
    outcome,
    probability,
    model_name,
    n_groups = 10
) {
  
  ranking <- rank(
    probability,
    ties.method = "first"
  )
  
  calibration_group <- ceiling(
    ranking /
      length(probability) *
      n_groups
  )
  
  calibration_group[
    calibration_group > n_groups
  ] <- n_groups
  
  grouped <- aggregate(
    cbind(
      Predicted = probability,
      Observed = outcome
    ) ~ calibration_group,
    FUN = mean
  )
  
  grouped$Model <- model_name
  grouped
}

calibration_points <- rbind(
  prepare_calibration_data(
    evaluation_data$outcome,
    evaluation_data$prob_clinical,
    "Model A: Clinical"
  ),
  prepare_calibration_data(
    evaluation_data$outcome,
    evaluation_data$prob_snp,
    "Model B: Clinical+SNPs"
  )
)

long_predictions <- rbind(
  data.frame(
    Outcome = evaluation_data$outcome,
    Probability =
      evaluation_data$prob_clinical,
    Model = "Model A: Clinical"
  ),
  data.frame(
    Outcome = evaluation_data$outcome,
    Probability =
      evaluation_data$prob_snp,
    Model = "Model B: Clinical+SNPs"
  )
)

calibration_plot <- ggplot() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_smooth(
    data = long_predictions,
    aes(
      x = Probability,
      y = Outcome,
      color = Model
    ),
    method = "loess",
    formula = y ~ x,
    se = TRUE
  ) +
  geom_point(
    data = calibration_points,
    aes(
      x = Predicted,
      y = Observed,
      color = Model
    ),
    size = 2.5
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    x = "Predicted probability",
    y = "Observed event proportion",
    color = "Model"
  ) +
  theme_minimal(base_size = 12)

calibration_plot



####################### Confusion matrices #############
clinical_cm <- confusionMatrix(
  data = factor(as.numeric(as.character(evaluation_data$Predicted.Clinical)),levels=c(0,1)),
  reference = factor(evaluation_data$outcome,levels=c(0,1)),
  positive = "1"
)

snp_cm <- confusionMatrix(
  data = factor(as.numeric(as.character(evaluation_data$Predicted.SNP)),levels=c(0,1)),
  reference = factor(evaluation_data$outcome,levels=c(0,1)),
  positive = "1"
)

#Plot here
# ============================================================
# 2. Prepare plotting data
# ============================================================

prepare_confusion_data <- function(
    confusion_matrix,
    model_name
) {
  
  plot_data <- as.data.frame(
    confusion_matrix$table
  )
  
  names(plot_data) <- c(
    "Predicted",
    "Observed",
    "Count"
  )
  
  # Readable outcome labels
  plot_data$Predicted <- factor(
    plot_data$Predicted,
    levels = c("0", "1"),
    labels = c(
      "Zero",
      "Non-zero"
    )
  )
  
  plot_data$Observed <- factor(
    plot_data$Observed,
    levels = c("0", "1"),
    labels = c(
      "Zero",
      "Non-zero"
    )
  )
  
  # Percentage within each observed outcome group
  plot_data$Percentage <- ave(
    plot_data$Count,
    plot_data$Observed,
    FUN = function(x) {
      100 * x / sum(x)
    }
  )
  
  # Cell label containing n and percentage
  plot_data$Label <- sprintf(
    "%d\n(%.1f%%)",
    plot_data$Count,
    plot_data$Percentage
  )
  
  plot_data$Model <- model_name
  
  plot_data
}

clinical_plot_data <- prepare_confusion_data(
  clinical_cm,
  "Model A"
)

snp_plot_data <- prepare_confusion_data(
  snp_cm,
  "Model B"
)

confusion_plot_data <- rbind(
  clinical_plot_data,
  snp_plot_data
)



# ============================================================
# 3. Plot confusion matrices
# ============================================================

confusion_matrix_plot <- ggplot(
  confusion_plot_data,
  aes(
    x = Predicted,
    y = Observed,
    fill = Count
  )
) +
  geom_tile(
    color = "white",
    linewidth = 1.2
  ) +
  geom_text(
    aes(label = Label),
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +
  facet_wrap(
    ~ Model,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = "#F0F6FC",
    high = "#2878B5",
    name = "Participants"
  ) +
  coord_equal() +
  labs(
    x = "Predicted SYNTAX-score group",
    y = "Observed SYNTAX-score group"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    legend.position = "bottom",
    
    panel.spacing = grid::unit(
      1.5,
      "lines"
    )
  )

confusion_matrix_plot


######################## DCA ##################
calculate_net_benefit <- function(
    outcome,
    probability,
    thresholds,
    model_name
) {
  
  n <- length(outcome)
  
  results <- lapply(
    thresholds,
    function(threshold) {
      
      predicted_positive <-
        probability >= threshold
      
      true_positive <- sum(
        predicted_positive &
          outcome == 1
      )
      
      false_positive <- sum(
        predicted_positive &
          outcome == 0
      )
      
      net_benefit <-
        true_positive / n -
        false_positive / n *
        threshold / (1 - threshold)
      
      data.frame(
        Threshold = threshold,
        Net_benefit = net_benefit,
        Strategy = model_name
      )
    }
  )
  
  do.call(rbind, results)
}
decision_thresholds <- seq(
  0.05,
  0.80,
  by = 0.01
)
dca_clinical <- calculate_net_benefit(
  outcome = evaluation_data$outcome,
  probability =
    evaluation_data$prob_clinical,
  thresholds = decision_thresholds,
  model_name = "Model A: Clinical"
)

dca_snp <- calculate_net_benefit(
  outcome = evaluation_data$outcome,
  probability =
    evaluation_data$prob_snp,
  thresholds = decision_thresholds,
  model_name = "Model B: Clinical+SNPs"
)


prevalence <- mean(
  evaluation_data$outcome == 1
)

dca_treat_all <- data.frame(
  Threshold = decision_thresholds,
  
  Net_benefit =
    prevalence -
    (1 - prevalence) *
    decision_thresholds /
    (1 - decision_thresholds),
  
  Strategy = "Treat all"
)

dca_treat_none <- data.frame(
  Threshold = decision_thresholds,
  Net_benefit = 0,
  Strategy = "Treat none"
)

dca_results <- rbind(
  dca_clinical,
  dca_snp,
  dca_treat_all,
  dca_treat_none
)
decision_curve_plot <- ggplot(
  dca_results,
  aes(
    x = Threshold,
    y = Net_benefit,
    color = Strategy,
    linetype = Strategy
  )
) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = 0,
    color = "grey70"
  ) +
  labs(
    x = "Threshold probability",
    y = "Net benefit",
    color = "Strategy",
    linetype = "Strategy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

decision_curve_plot
