rm(list=ls())
library(caret)
library(SNPassoc)
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Models/Models_RF_AUC_None.Rda")



#Pre-processing 
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_clinical_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Test_clinical_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_snp_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Test_snp_zero.Rda")



colnames(train_clinical_ls[[1]])
predictors(final_models$Clinical_RF_1)


########## Clinical::: For clinical we want min-max and one-hot-encoding ################
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
res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[1]],test_clinical_ls[[1]],train_snp_ls[[1]],test_snp_ls[[1]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)

clinical_train<-res_preprocessing$ClinicalTrain
clinical_test<-res_preprocessing$ClinicalTest
snp_train<-res_preprocessing$SNPtrain
snp_test<-res_preprocessing$SNPtest

# Read for prediction 









############## SNP:: the same pre-processing (since we have the inheritance encoding) but also interaction analysis #############
############## So the only thing that we want for SNP extra is the encoding -- inheritance modes -- ######################
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

res_interaction<-interaction_analysis(train_snp_ls[[1]],test_snp_ls[[1]])
snp_train<-res_interaction$Train
snp_test<-res_interaction$Test
res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[1]],test_clinical_ls[[1]],snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)

clinical_train<-res_preprocessing$ClinicalTrain
clinical_test<-res_preprocessing$ClinicalTest
snp_train<-res_preprocessing$SNPtrain
snp_test<-res_preprocessing$SNPtest







############### In case of not having the numerical encoding ##################
snp_numeric_convert <-function(dataset,snp_columns,levels){
  dataset_no_snp<-dataset[! colnames(dataset) %in% snp_columns]
  dataset_snp<-dataset[colnames(dataset) %in% snp_columns]
  for(snp in snp_columns){
    levels_snp<-levels[[snp]]
    if(length(levels_snp)==2){
      level_0<-strsplit(levels_snp[1],split="-")[[1]]
      level_1<-strsplit(levels_snp[2],split="-")[[1]]
      dataset_snp[[snp]]<-as.factor(ifelse(dataset_snp[[snp]] %in% level_0,0,1))
    }else{
      level_0<-levels_snp[1]
      level_1<-levels_snp[2]
      dataset_snp[[snp]]<-as.factor(ifelse(dataset_snp[[snp]] %in% level_0,0,ifelse(dataset_snp[[snp]] %in% level_1,1,2)))
    }
    
  }
  
  df<-cbind(dataset_no_snp,dataset_snp)
  return(df)
}














############### For encoding - SNP level extraction for zero and count part ##################


optimal_encoding<-function(train,snp_cols,target_col){
  snp_cols_wanted<-which(colnames(train) %in% snp_cols)
  trainSNP<-setupSNP(data=train,colSNPs=snp_cols_wanted,sep="")
  
  levels_modelling<-list()
  snps<-c()
  genetic_modelling<-"Optimal"
  
  
  ind_categorical<-9
  ind_numerical<-7
  
  if(is.numeric(train[[target_col]])){
    wanted<-ind_numerical
  }else{
    wanted<-ind_categorical
  }
  
  
  
  for(snp in snp_cols_wanted){
    form<-as.formula(paste(target_col,' ~',colnames(train)[snp]))
    res<-association(form,data=trainSNP)
    if(length(res[,1])==15){
      snps<-c(snps,snp)
      levels_co<-names(res[,1])[2:4]
      levels_do<-names(res[,1])[6:7]
      levels_re<-names(res[,1])[9:10]
      levels_over<-names(res[,1])[12:13]
      
      if(genetic_modelling=='Codominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_co
      }else if(genetic_modelling=='Dominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_do
      }else if(genetic_modelling=='Recessive'){
        levels_modelling[[colnames(train)[snp]]]<-levels_re
      }else if(genetic_modelling=='Optimal'){
        index<-which.min(res[1:12,wanted])
        if(index==2){
          levels_opt<-levels_co
        }else if(index==6){
          levels_opt<-levels_do
        }else if(index==9){
          levels_opt<-levels_re
        }else if(index==12){
          levels_opt<-levels_over
        }
        levels_modelling[[colnames(train)[snp]]]<-levels_opt
      }
      
      
    }
  }
  
  return(levels_modelling)
}

### SNP level extraction from folds ###
#Load whole dataset
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_dataset.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Test_inds.Rda")

#Load final CV --- Zero part
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_snp_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Test_snp_zero.Rda")

#1. Zero levels 
levels_snp_zero<-list()
for(i in 1:10){
  train_set<-dataset[-list_indexes_test[[i]],]
  train_set<-train_set[!is.na(train_set$SYNTAX.SCORE),]
  
  
  rs_wanted<-colnames(train_snp_ls[[i]])[(grepl("rs",colnames(train_snp_ls[[i]])))]
  
  res<-optimal_encoding(train_set,rs_wanted,target_col='Binary.SS')
  
  levels_snp_zero[[i]]<-res
}
save(levels_snp_zero,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Levels_SNP_zero.Rda')



#Load final CV --- Count part
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Train_snp_count.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Test_snp_count.Rda")



#2. Count levels
levels_snp_count<-list()
for(i in 1:10){
  train_set<-dataset[-list_indexes_test[[i]],]
  train_set<-train_set[!is.na(train_set$SYNTAX.SCORE),]
  train_set<-train_set[train_set$SYNTAX.SCORE>0,]
  
  rs_wanted<-colnames(train_snp_ls[[i]])[(grepl("rs",colnames(train_snp_ls[[i]])))]
  
  res<-optimal_encoding(train_set,rs_wanted,target_col='SYNTAX.SCORE')
  
  levels_snp_count[[i]]<-res
}
save(levels_snp_count,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Levels_SNP_count.Rda')





