######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)
outer_fold<-1
load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Models/Models_RF_AUC_None.Rda'))



load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_clinical_CV_BA_DM.Rda'))
load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_clinical_CV_BA_DM.Rda'))
load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_snp_CV_BA_DM.Rda'))
load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_snp_CV_BA_DM.Rda'))

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
#MLmetrics::F1_Score(as.factor(test_all),preds_snp,positive=1)  
#MLmetrics::F1_Score(as.factor(test_all),preds_clinical,positive=1)  

  
######### AUC ROC and AUC PR per fold
rocs_snp<-c()
rocs_clinical<-c()
pr_snp<-c()
pr_clinical<-c()
data_roc<-data.frame(Fold=folds_arr,Actual=test_all,Probs.Clinical=probs_clinical,Probs.SNP=probs_snp)
for(fold in 1:10){
  data_fold<-data_roc[data_roc$Fold==fold,]
  
  roc_clinical<-roc(data_fold$Actual,data_fold$Probs.Clinical)
  
  roc_snp<-roc(data_fold$Actual,data_fold$Probs.SNP)
  
  rocs_snp<-c(rocs_snp,round(as.numeric(roc_snp$auc),3))
  
  rocs_clinical<-c(rocs_clinical,round(as.numeric(roc_clinical$auc),3))
  
  
  pr_cl<-pr.curve(scores.class0 = data_fold$Probs.Clinical,weights.class0 = data_fold$Actual)
  pr_clinical<-c(pr_clinical,round(as.numeric(pr_cl$auc.integral),3))
  
  pr_sn<-pr.curve(scores.class0=data_fold$Probs.SNP,weights.class0 = data_fold$Actual)
  pr_snp<-c(pr_snp,round(as.numeric(pr_sn$auc.integral),3))
  }
AUC_per_fold<-data.frame(Fold=1:10,AUCROC.Clinical=rocs_clinical,AUCROC.SNP=rocs_snp,AUCPR.Clinical=pr_clinical,AUCPR.SNP=pr_snp)
#save(AUC_per_fold,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Results/AUC_per_fold.Rda')
#libraries
library(PRROC)




roc1 <- pROC::roc(test_all,probs_clinical)
roc2 <- pROC::roc(test_all,probs_snp)
roc1$auc
roc2$auc
######## Delong
res<-roc.test(roc1,roc2,method="delong")

#################### Custom code for ROC test #################
library(pROC)

# y_true: binary vector (0/1), length N
# pred_A, pred_B: pooled predictions from stratified 10-fold CV, length N

# Observed difference
y_true<-test_all
pred_A<-probs_snp
pred_B<-probs_clinical
roc_A <- roc(y_true, pred_A, quiet = TRUE)
roc_B <- roc(y_true, pred_B, quiet = TRUE)
delta_obs <- auc(roc_A) - auc(roc_B)

# Permutation test
set.seed(123)
n_perm <- 5000
delta_perm <- numeric(n_perm)

for (i in 1:n_perm) {
  swap <- rbinom(length(y_true), 1, 0.5) == 1
  perm_A <- ifelse(swap, pred_B, pred_A)
  perm_B <- ifelse(swap, pred_A, pred_B)
  
  auc_A_perm <- auc(roc(y_true, perm_A, quiet = TRUE))
  auc_B_perm <- auc(roc(y_true, perm_B, quiet = TRUE))
  delta_perm[i] <- auc_A_perm - auc_B_perm
}

# Two-sided p-value
p_val <- (1 + sum(abs(delta_perm) >= abs(delta_obs))) / (1 + n_perm)

list(
  AUC_A = auc(roc_A),
  AUC_B = auc(roc_B),
  delta_obs = delta_obs,
  p_value = p_val
)





n_boot <- 5000
delta_boot <- numeric(n_boot)

# Identify the two outcome groups
outcome_groups <- split(seq_along(y_true), y_true)

for (i in seq_len(n_boot)) {
  
  # Resample participants with replacement within each outcome group
  boot_index <- unlist(
    lapply(
      outcome_groups,
      function(index) sample(index, length(index), replace = TRUE)
    ),
    use.names = FALSE
  )
  
  y_boot <- y_true[boot_index]
  pred_A_boot <- pred_A[boot_index]  # SNP-enhanced model
  pred_B_boot <- pred_B[boot_index]  # Clinical-only model
  
  roc_A_boot <- roc(
    y_boot,
    pred_A_boot,
    levels = roc_A$levels,
    direction = roc_A$direction,
    quiet = TRUE
  )
  
  roc_B_boot <- roc(
    y_boot,
    pred_B_boot,
    levels = roc_B$levels,
    direction = roc_B$direction,
    quiet = TRUE
  )
  
  delta_boot[i] <- as.numeric(
    auc(roc_A_boot) - auc(roc_B_boot)
  )
}

# Percentile bootstrap 95% CI for the paired AUC difference
delta_ci <- quantile(
  delta_boot,
  probs = c(0.025, 0.975),
  na.rm = TRUE,
  names = FALSE
)
results <- data.frame(
  auc_snp = as.numeric(auc(roc_A)),
  auc_clinical = as.numeric(auc(roc_B)),
  paired_difference = as.numeric(delta_obs),
  ci_lower = delta_ci[1],
  ci_upper = delta_ci[2],
  test_statistic = abs(as.numeric(delta_obs)),
  permutation_p_value = p_val
)

#$AUC_A
#Area under the curve: 0.7599

#$AUC_B
#Area under the curve: 0.7311

#$delta_obs
# 0.0288564

#$p_value
#0.03779244

######## For PR
pr_curve_cl<-PRROC::pr.curve(probs_clinical,weights.class0=test_all,curve=TRUE,max.compute = TRUE,min.compute = TRUE)
pr_curve_snp<-PRROC::pr.curve(probs_snp,weights.class0=test_all,curve=TRUE)
y_true<-test_all
pred_A<-probs_snp
pred_B<-probs_clinical

pr_A <- pr.curve(pred_A,weights.class0=y_true)
pr_B <- pr.curve(pred_B,weights.class0=y_true)
delta_obs <- pr_A$auc.integral - pr_B$auc.integral

# Permutation test
set.seed(123)
n_perm <- 5000
delta_perm <- numeric(n_perm)

for (i in 1:n_perm) {
  swap <- rbinom(length(y_true), 1, 0.5) == 1
  perm_A <- ifelse(swap, pred_B, pred_A)
  perm_B <- ifelse(swap, pred_A, pred_B)
  
  auc_A_perm <- pr.curve(perm_A,weights.class0 = y_true)$auc.integral
  auc_B_perm <-pr.curve(perm_B,weights.class0=y_true)$auc.integral
  delta_perm[i] <- auc_A_perm - auc_B_perm
}

# Two-sided p-value
p_val <- (1 + sum(abs(delta_perm) >= abs(delta_obs))) / (1 + n_perm)

list(
  AUC_A = pr_A$auc.integral,
  AUC_B = pr_B$auc.integral,
  delta_obs = delta_obs,
  p_value = p_val
)

#$AUC_A
#0.8828067

#$AUC_B
#[1] 0.8639102

#$delta_obs
#[1] 0.01889644

#$p_value
#0.04219156

########################### End of custom code ####################

pr_curve_cl<-PRROC::pr.curve(probs_clinical,weights.class0=test_all,curve=TRUE,max.compute = TRUE,min.compute = TRUE)
pr_curve_snp<-PRROC::pr.curve(probs_snp,weights.class0=test_all,curve=TRUE)
pr_curve_cl$auc.integral
pr_curve_snp$auc.integral

roc11<-pROC::roc(test_all,probs_clinical)
pROC::ci.auc(roc11,method="bootstrap")

roc22<-pROC::roc(test_all,probs_snp)
pROC::ci.auc(roc22,method="bootstrap")


source('C:/Users/30697/Desktop/Project_Indices/Code/pr_test.R')
res<-pr.test(test_all,probs_clinical,probs_snp,boot.n=2000)


library(precrec)
class <- test_all
prob <- probs_clinical

inds<-seq(1, 799, by = 158)-1
inds1<-1:inds[2]

next_index<-inds[2]+1
inds2<-next_index:inds[3]


next_index<-inds[3]+1
inds3<-next_index:inds[4]

next_index<-inds[4]+1
inds4<-next_index:inds[5]

next_index<-inds[5]+1
inds5<-next_index:inds[6]

next_index<-inds[6]+1
inds6<-next_index:inds[7]

next_index<-inds[7]+1
inds7<-next_index:inds[8]

next_index<-inds[8]+1
inds8<-next_index:inds[9]

next_index<-inds[9]+1
inds9<-next_index:inds[10]

next_index<-inds[10]+1
inds10<-next_index:inds[11]


ind<-list(inds1,inds2,inds3,inds4,inds5)#,inds6,inds7,inds8,inds9,inds10)

scores1<-prob[inds1]
labels1<-class[inds1]

scores2<-prob[inds2]
labels2<-class[inds2]

scores3<-prob[inds3]
labels3<-class[inds3]


scores4<-prob[inds4]
labels4<-class[inds4]


scores5<-prob[inds5]
labels5<-class[inds5]

scores6<-prob[inds6]
labels6<-class[inds6]

scores7<-prob[inds7]
labels7<-class[inds7]

scores8<-prob[inds8]
labels8<-class[inds8]

scores9<-prob[inds9]
labels9<-class[inds9]

#scores10<-prob[inds10]
#labels10<-class[inds10]

scores<-join_scores(scores1,scores2,scores3,scores4,scores5)
#,
                    #scores6,scores7,scores8,scores9)#,scores10)
labels<-join_labels(labels1,labels2,labels3,labels4,labels5)
#,
                    #labels6,labels7,labels8,labels9)#,labels10)



# Create the mmdata object
mdat <- mmdata(scores, labels, dsids = c(1,2,3,4,5))#,6,7,8,9))#,10))
curves<-evalmod(mdat)
r<-auc_ci(curves,alpha=0.05)
r


############## Per fold ##############
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Results/Classifier_results_RF.Rda')


wanted<-res_all[res_all$Sampling=='None',]
wanted_1<-wanted[21:40,]
median(wanted_1[wanted_1$Data_Type=='Clinical','Precision'])

median(wanted_1[wanted_1$Data_Type=='ClinicalSNP','Precision'])









################# Boruta selected ####################
all_clinical<-c()
all_snp<-c()
predicted<-list()
for(fold in 1:10){
  test<-train_snp_list[[fold]]
  predictors_all <- colnames(test)[!colnames(test) %in% c('target')]
  rs_col<-which(grepl("rs",predictors_all))
  clinical_predicted<-predictors_all[-rs_col]
  snp_predicted<-predictors_all[rs_col]
  all_clinical<-c(all_clinical,clinical_predicted)
  all_snp<-c(all_snp,snp_predicted)
  predicted[[fold]]<-predictors_all
  }

unique(all_clinical)



load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Train_sets.Rda')
all_clinical<-c()
all_snp<-c()
predicted<-list()
for(fold in 1:10){
  test<-train_list[[paste0('SNP_AUC_None_',fold)]]
  predictors_all <- colnames(test)[!colnames(test) %in% c('Binary.SS','SYNTAX.SCORE')]
  rs_col<-which(grepl("rs",predictors_all))
  clinical_predicted<-predictors_all[-rs_col]
  snp_predicted<-predictors_all[rs_col]
  all_clinical<-c(all_clinical,clinical_predicted)
  all_snp<-c(all_snp,snp_predicted)
  predicted[[fold]]<-predictors_all
}

unique(all_clinical)

