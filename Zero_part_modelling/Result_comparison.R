load('./Delong_results_RF.Rda')
important<-df[df$P.value < 0.05,]
save(important,file='./Important_diff_Delong.Rda')


load('./WilcoxonResults_HighRiskDMBA.Rda')
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
save(wanted_under,file='./Wanted_LowRisk.Rda')

