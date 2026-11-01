rm(list=ls())
load("./Whole_dataset.Rda")
dataset<-dataset[!is.na(dataset$SYNTAX.SCORE),]
predictors_all<-c("Gender","HYPERTENSION","SMOKING","GLU","NEU%",
                  "RATIO1","RATIO3","RATIO4","CHEST PAIN","ATYPICAL ANGINA SYMPTOMS",
                  "ATRIAL FIBRILLATION","RATIO2","EOS%","BASO%",
                  "rs11984041","rs1330344","rs2023938","rs2107595","rs216172","rs25882",
                  "rs3732379","rs4516035","rs56336142","rs9818870",
                  "DIABETES MELLITUS","CRUSADE SCORE","QRS DURATION ms","GLU",
                  "UREA","GALECTIN","AUTOIMMUNE DISEASE","GFR",
                  
                  
                  "rs10871777","rs12526453","rs1332844","rs1333042","rs17514846",
                  "rs2046934","rs2487928","rs4593108","rs4845625","rs55730499",
                  "rs56289821","rs663129","rs6689306","rs6801273","rs870142")

dataset_wanted<-dataset[c(predictors_all,"SYNTAX.SCORE")]

dataset_zero<-dataset_wanted[dataset_wanted$SYNTAX.SCORE==0,]
dataset_nonzero<-dataset_wanted[dataset_wanted$SYNTAX.SCORE>0,]
dataset_low<-dataset_nonzero[dataset_nonzero$SYNTAX.SCORE<23,]
dataset_high<-dataset_nonzero[dataset_nonzero$SYNTAX.SCORE>=23,]

counts<-c()
perc<-c()

zeros<-c()
nonzeros<-c()
lows<-c()
highs<-c()
for(predictor in predictors_all){
  zeros<-c(zeros,paste0(sum(is.na(dataset_zero[predictor])),"(",round(sum(is.na(dataset_zero[predictor]))/nrow(dataset_zero)*100,2),"%)"))
  nonzeros<-c(nonzeros,paste0(sum(is.na(dataset_nonzero[predictor])),"(",round(sum(is.na(dataset_nonzero[predictor]))/nrow(dataset_nonzero)*100,2),"%)"))
  lows<-c(lows,paste0(sum(is.na(dataset_low[predictor])),"(",round(sum(is.na(dataset_low[predictor]))/nrow(dataset_low)*100,2),"%)"))
  highs<-c(highs,paste0(sum(is.na(dataset_high[predictor])),"(",round(sum(is.na(dataset_high[predictor]))/nrow(dataset_high)*100,2),"%)"))
}
df_missing<-data.frame(Predictors=predictors_all,Zero=zeros,NonZero=nonzeros,LowRisk=lows,HighRisk=highs)
