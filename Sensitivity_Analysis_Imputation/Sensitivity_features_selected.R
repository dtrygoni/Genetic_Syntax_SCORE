rm(list=ls())
######## Sensitivity analysis for Features Selected among imputation / original pipeline
library(caret)
###1. Zero-part
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_Imputation/Zero_part_classification/Train_snp_zero.Rda")

features_snp<-c()
features_clinical<-c()

for(i in 1:10){
  predictors_i<-colnames(train_snp_ls[[i]])
  features_snp<-c(features_snp,predictors_i)
}


load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_snp_zero.Rda")
features_snp_original<-c()
for(i in 1:10){
  predictors_i<-colnames(train_snp_ls[[i]])
  features_snp_original<-c(features_snp_original,predictors_i)
}

res_orig<-table(features_snp_original)
res_imp<-table(features_snp)



names_all<-unique(c(names(res_orig),names(res_imp)))
Orig_num<-c()
Imp_num<-c()
for(name in names_all){
  if(name %in% names(res_orig)){
    Orig_num<-c(Orig_num,as.numeric(res_orig[name]))
  }else{
    Orig_num<-c(Orig_num,0)
  }
  
  if(name %in% names(res_imp)){
    Imp_num<-c(Imp_num,as.numeric(res_imp[name]))
  }else{
    Imp_num<-c(Imp_num,0)
  }
}
df<-data.frame(Predictors=names_all,Number.Selected.Original=Orig_num,Number.Selected.Imputation=Imp_num)


rm(list=ls())
###2. Count-part
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_Imputation/Count_part_regression/Train_snp_count.Rda")

features_snp<-c()
features_clinical<-c()

for(i in 1:10){
  predictors_i<-colnames(train_snp_ls[[i]])
  features_snp<-c(features_snp,predictors_i)
}


load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Train_snp_count.Rda")
features_snp_original<-c()
for(i in 1:10){
  predictors_i<-colnames(train_snp_ls[[i]])
  features_snp_original<-c(features_snp_original,predictors_i)
}

res_orig<-table(features_snp_original)
res_imp<-table(features_snp)



names_all<-unique(c(names(res_orig),names(res_imp)))
Orig_num<-c()
Imp_num<-c()
for(name in names_all){
  if(name %in% names(res_orig)){
    Orig_num<-c(Orig_num,as.numeric(res_orig[name]))
  }else{
    Orig_num<-c(Orig_num,0)
  }
  
  if(name %in% names(res_imp)){
    Imp_num<-c(Imp_num,as.numeric(res_imp[name]))
  }else{
    Imp_num<-c(Imp_num,0)
  }
}
df<-data.frame(Predictors=names_all,Number.Selected.Original=Orig_num,Number.Selected.Imputation=Imp_num)
