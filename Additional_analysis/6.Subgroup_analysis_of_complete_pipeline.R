rm(list=ls())
load("./Whole_dataset.Rda")
dataset<-dataset[!is.na(dataset$SYNTAX.SCORE),]


inds<-list()
for(i in 1:10){
  library(dplyr)
  full_with_id <- dataset %>% mutate(full_row = row_number()) 
  matched <- test_clinical_ls[[i]] %>% left_join( full_with_id, by = intersect(names( test_clinical_ls[[i]]), names(dataset)) )
  inds[[i]]<-matched$full_row
}



inds_wanted<-unlist(inds)

dataset_inds<-dataset[inds_wanted,]

dataset_inds$Gender
dataset_inds$`DIABETES MELLITUS`
dataset_inds$AGE_group<-ifelse(dataset_inds$AGE>=65,1,0)


load("./Full_model_results.Rda")
#### Gender ###
male_inds<-which(dataset_inds$Gender==1)
wanted_res<-res[male_inds,]
median(abs(wanted_res$Target-wanted_res$SNP.Predicted))
mean(abs(wanted_res$Target-wanted_res$SNP.Predicted))
median(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted))



male_inds<-which(dataset_inds$Gender==0)
wanted_res<-res[male_inds,]
median(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
median(abs(wanted_res$Target-wanted_res$SNP.Predicted))
mean(abs(wanted_res$Target-wanted_res$SNP.Predicted))


#### AGE ###
male_inds<-which(dataset_inds$AGE_group==1)
wanted_res<-res[male_inds,]
median(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
median(abs(wanted_res$Target-wanted_res$SNP.Predicted))
mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
mean(abs(wanted_res$Target-wanted_res$SNP.Predicted))



male_inds<-which(dataset_inds$AGE_group==0)
wanted_res<-res[male_inds,]
median(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
median(abs(wanted_res$Target-wanted_res$SNP.Predicted))
mean(abs(wanted_res$Target-wanted_res$SNP.Predicted))



### Diabetes mellitus ###
male_inds<-which(dataset_inds$`DIABETES MELLITUS`==1)
wanted_res<-res[male_inds,]
mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
mean(abs(wanted_res$Target-wanted_res$SNP.Predicted))
median(abs(wanted_res$Target-wanted_res$Clinical.Predicted))
median(abs(wanted_res$Target-wanted_res$SNP.Predicted))



male_inds<-which(dataset_inds$`DIABETES MELLITUS`==0)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)


#### Smoking ###
male_inds<-which(dataset_inds$SMOKING==1)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)


male_inds<-which(dataset_inds$SMOKING==0)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)



### Hypertension ###
male_inds<-which(dataset_inds$HYPERTENSION==1)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)

male_inds<-which(dataset_inds$HYPERTENSION==0)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)



### Dyslipidaemia ###
male_inds<-which(dataset_inds$DYSLIPIDAEMIA==1)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)


male_inds<-which(dataset_inds$DYSLIPIDAEMIA==0)
wanted_res<-res[male_inds,]
round(mean(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(mean(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$Clinical.Predicted)),3)
round(median(abs(wanted_res$Target-wanted_res$SNP.Predicted)),3)
