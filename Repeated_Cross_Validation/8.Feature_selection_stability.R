


############### For zero/non-zero component predictors ###############

predictors<-c()
for(exp in 1:10){
  load(paste0('./Repeated_Cross_Validation/Fold',exp,'/Zero_non_zero_component/Test_snp_zero.Rda'))
  
  for(fold in 1:10){
    predictors<-c(predictors,colnames(test_snp_ls[[fold]]))
  }
}

print("Zero/non-zero predictors")
table(predictors)



############### For positive score component predictors ###############
predictors_pos<-c()
for(exp in 1:10){
  load(paste0('./Repeated_Cross_Validation/Fold',exp,'/Positive_score_component/Test_snp_count.Rda'))
  
  for(fold in 1:10){
    predictors_pos<-c(predictors_pos,colnames(test_snp_ls[[fold]]))
  }
}

print("Positive score predictors")
table(predictors_pos)
