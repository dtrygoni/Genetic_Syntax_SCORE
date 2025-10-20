load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Results/Delong_results_RF.Rda')
important<-df[df$P.value < 0.05,]
save(important,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments_Count_Classifier/Regression_results/Important_diff_Delong.Rda')


load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments_Count_Classifier/Regression_results/WilcoxonResults_HighRiskDMBA.Rda')
wil_all$Diff<-wil_all$Clinical-wil_all$SNP
important<-wil_all[wil_all$P.value<0.05,]

wanted<-NULL
for(i in 1:nrow(important)){
  diff<-important[i,]
  if(diff$Clinical>diff$SNP){
    wanted<-rbind(wanted,diff)
  }
}
wanted_under<-wanted[wanted$Resampling=='under',]
save(wanted_under,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments_Count_Classifier/Regression_results/Wanted_LowRisk.Rda')
