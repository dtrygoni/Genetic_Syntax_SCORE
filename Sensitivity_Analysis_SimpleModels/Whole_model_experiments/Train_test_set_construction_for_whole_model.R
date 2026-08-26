########### Check that zero test/train are the same for clinical/clinical+SNP ####

rm(list=ls())
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Train_clinical_count.Rda")
train_clinical_count<-train_clinical_ls

load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_clinical_zero.Rda")
train_clinical_zero<-train_clinical_ls


load("C:/Users/dtryg/Desktop/GESS_experiments_code/Count_part_regression/Train_snp_count.Rda")
train_snp_count<-train_snp_ls

load("C:/Users/dtryg/Desktop/GESS_experiments_code/Zero_part_classification/Train_snp_zero.Rda")
train_snp_zero<-train_snp_ls

########### Construction of train and test sets based on whole model (zero + count) #############

#Features of the Dataset
feature_names<-list()
for(i in 1:length(train_clinical_zero)){
  feats_zero<-colnames(train_clinical_zero[[i]])
  feats_count<-colnames(train_clinical_count[[i]])
  feats_snp_zero<-colnames(train_snp_zero[[i]])
  feats_snp_count<-colnames(train_snp_count[[i]])
  
  
  all_feat_names<-unique(c(feats_zero,feats_count,feats_snp_zero,feats_snp_count))
  all_feat_names<-all_feat_names[! all_feat_names %in% 'target']
  feature_names[[i]]<-all_feat_names
}


#Features of clinical models
feature_names_clinical<-list()
for(i in 1:length(train_clinical_zero)){
  feats_zero<-colnames(train_clinical_zero[[i]])
  feats_count<-colnames(train_clinical_count[[i]])

  
  all_feat_names<-unique(c(feats_zero,feats_count))
  all_feat_names<-all_feat_names[! all_feat_names %in% 'target']
  feature_names_clinical[[i]]<-all_feat_names
}

#Features of snp models
feature_names_snp<-list()
for(i in 1:length(train_snp_zero)){

  feats_snp_zero<-colnames(train_snp_zero[[i]])
  feats_snp_count<-colnames(train_snp_count[[i]])
  
  all_feat_names<-unique(c(feats_snp_zero,feats_snp_count))
  all_feat_names<-all_feat_names[! all_feat_names %in% 'target']
  feature_names_snp[[i]]<-all_feat_names
}

#Inds of the train sets
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Test_inds.Rda")
#Initial dataset
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_dataset.Rda")


train_set_list<-list()
test_set_list<-list()
for(i in 1:10){
  train_set<-dataset[-list_indexes_test[[i]],feature_names[[i]]]
  train_set_f<-train_set[complete.cases(train_set),]
  
  test_set<-dataset[list_indexes_test[[i]],feature_names[[i]]]
  test_set_f<-test_set[complete.cases(test_set),]
  
  train_set_list[[i]]<-train_set_f
  test_set_list[[i]]<-test_set_f
}


########### Clinical sets #############
train_clinical_list<-list()
test_clinical_list<-list()
for(i in 1:10){
  wanted_train<-train_set_list[[i]][feature_names_clinical[[i]]]
  wanted_test<-test_set_list[[i]][feature_names_clinical[[i]]]
  
  
  train_clinical_list[[i]]<-wanted_train
  test_clinical_list[[i]]<-wanted_test
}



########## SNP sets ################
train_snp_list<-list()
test_snp_list<-list()
for(i in 1:10){
  wanted_train<-train_set_list[[i]][feature_names_snp[[i]]]
  wanted_test<-test_set_list[[i]][feature_names_snp[[i]]]
  
  
  train_snp_list[[i]]<-wanted_train
  test_snp_list[[i]]<-wanted_test
}


save(train_clinical_list,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_model_experiments/Train_clinical.Rda')
save(test_clinical_list,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_model_experiments/Test_clinical.Rda')

save(train_snp_list,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_model_experiments/Train_snp.Rda')
save(test_snp_list,file='C:/Users/dtryg/Desktop/GESS_experiments_code/Whole_model_experiments/Test_snp.Rda')

