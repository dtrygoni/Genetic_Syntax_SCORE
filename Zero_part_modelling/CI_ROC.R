######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)
load('./Models_RF_AUC_None.Rda')
load('./Test_clinical_AUC_None.Rda')
load('./Train_clinical_AUC_None.Rda')

load('./Test_snp_AUC_None.Rda')
load('./Train_snp_AUC_None.Rda')


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
  test_clinical <-test_clinical_list[[jj]]
  
  model_clinical<-final_models[[paste0('Clinical_RF_',jj)]]
  
  train_clinical<-train_clinical_list[[jj]]
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
  test_snp <-test_snp_list[[jj]] 
  model_snp<-final_models[[paste0('SNP_RF_',jj)]]
  
  train_snp<-train_snp_list[[jj]]
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
save(AUC_per_fold,file='./Results/AUC_per_fold.Rda')
#libraries
library(PRROC)




roc1 <- roc(test_all,probs_clinical)
roc2 <- roc(test_all,probs_snp)
roc1$auc
roc2$auc
######## Delong
res<-roc.test(roc1,roc2,method="boot")

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



load('./Train_sets.Rda')
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


