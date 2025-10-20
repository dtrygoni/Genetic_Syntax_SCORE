############################### FULL MODEL PREDICTION #####################################
library(recipes)
library(caret)
library(SNPassoc)

### Utility function
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

GESS_prediction<-function(dataset){
######## Create two datasets: zero and count

#####   Zero part #####
load('./Binary_res/SNP_feats_selected.Rda')
load('./Binary_res/Clinical_feats_selected.Rda')
dataset_zero<-dataset[colnames(dataset) %in% c(cl_preds,snp_preds)]

#SNP convert levels
load('./Binary_res/SNP_All_levels.Rda')
snp_cols<-which(grepl("rs",colnames(dataset_zero)))
dataset_zero<-setupSNP(dataset_zero,snp_cols,sep="")
dataset_zero_converted<-snp_numeric_convert(dataset_zero,colnames(dataset_zero)[snp_cols],levels_all_snp)

#One-hot-encoding
load('./Binary_res/Zero_part_OneHotEncoder.Rda')
dataset_zero_converted$target<-rep(1,nrow(dataset_zero_converted))
dataset_zero_onehot<-data.frame(predict(one_hot_encoder,newdata=dataset_zero_converted))
#Numerical scale
load('./Binary_res/Zero_part_NumericalScaler.Rda')
dataset_zero_scaled <- predict(numerical_scaler, newdata = dataset_zero_onehot)


dataset_zero_scaled<-dataset_zero_scaled[!colnames(dataset_zero_scaled) %in% c('target')]
#Zero-part prediction
load('./Binary_res/Zero_part_model.Rda')




#zeroResponse<-predict(zero.part_model,dataset_zero)

#####  Count part #####
load('./Count_res/SNP_feats_selected.Rda')
load('./Count_res/Clinical_feats_selected.Rda')
dataset_count<-dataset[colnames(dataset) %in% c(cl_preds,snp_preds)]

#SNP convert levels
load('./Count_res/SNP_All_levels.Rda')
snp_cols<-which(grepl("rs",colnames(dataset_count)))
dataset_count<-setupSNP(dataset_count,snp_cols,sep="")
dataset_count_converted<-snp_numeric_convert(dataset_count,colnames(dataset_count)[snp_cols],levels_all_snp)

#One-hot-encoding
load('./Count_res/Count_part_OneHotEncoder.Rda')
dataset_count_converted$target<-rep(1,nrow(dataset_count_converted))
dataset_count_converted$SYNTAX.SCORE<-rep(1,nrow(dataset_count_converted))
dataset_count_onehot<-data.frame(predict(one_hot_encoder,newdata=dataset_count_converted))

#Numerical scale
load('./Count_res/Count_part_NumericalScaler.Rda')
dataset_count_scaled <- predict(numerical_scaler, newdata = dataset_count_onehot)
dataset_count_scaled<-dataset_count_scaled[!colnames(dataset_count_scaled) %in% c('target','SYNTAX.SCORE')]


#Count-part models
load('./Count_res/Count_part_classifier.Rda')
load('./Count_res/Count_part_regressor_low.Rda')
load('./Count_res/Count_part_regressor_high.Rda')




#Common dataset
count_inds<-which(complete.cases(dataset_count_scaled))
zero_inds<-which(complete.cases(dataset_zero_scaled))


dataset_zero<-dataset_zero_scaled[count_inds[count_inds %in% zero_inds],]
dataset_count<-dataset_count_scaled[zero_inds[zero_inds %in% count_inds],]



zeroResponse<-predict(zero.part_model,dataset_zero)

countClassifier<-predict(count.part_classifier_model,dataset_count)
regressorlow<-predict(count.part_regressor_low,dataset_count)[[1]]
regressorhigh<-predict(count.part_regressor_high,dataset_count)[[1]]


finalResponse<-round(ifelse(zeroResponse == 'neg',0,ifelse(countClassifier=='neg',regressorlow,regressorhigh)))
return(finalResponse)
}
