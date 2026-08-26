############################################## Model training -- Count Part Random Forest Classifier / Regressors ########################
##Libraries
library(missForestPredict)
library(Boruta)
library(caret)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(xgboost)
library(PRROC)
library(precrec)
library(mltools)
library(precrec)
library(pROC)
library(SNPassoc)
library(tibble)
library(recipes)
library(MLmetrics)
set.seed(123)
xgb.set.config(verbosity = 0)

model_arr<-c()
type_arr<-c()
r2_arr<-c()
mae_arr<-c()
mre_arr<-c()
rmse_arr<-c()
fold_arr<-c()
final_models<-list()
train_list<-list()
test_list<-list()
###### load
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/SNP_dataset.Rda')


##### Functions
### Helper functions
summary_fnc <- function(data, lev = NULL, model = NULL) {
  # Ensure prediction and true labels are factors
  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = lev)
  
  cm <- caret::confusionMatrix(pred, obs)
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  acc <- MLmetrics::Accuracy(pred, obs)
  f1 <- FBeta_Score(obs, pred, positive = "pos", beta = 1)
  f2 <- FBeta_Score(obs, pred, positive = "pos", beta = 2)
  f05 <- FBeta_Score(obs, pred, positive = "pos", beta = 0.5)
  
  # AUC Calculation (assumes the column name for positive class is 'pos')
  if ("pos" %in% colnames(data)) {
    auc <- pROC::roc(response = obs, predictor = data[["pos"]], levels = rev(lev), direction = ">")$auc
  } else {
    auc <- NA
  }
  
  out <- c(
    bal.Acc = bal_acc,
    F1Score = f1,
    Spec = spec,
    Sens = sens,
    Accuracy = acc,
    F2Score = f2,
    F05Score = f05,
    AUC = auc
  )
  return(out)
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
Preprocessing_fnc_whole<-function(clinical_train,target_col='target',factor_col_opt='One-hot',numerical_col_opt='None'){
  ########### Factor preprocess
  if(factor_col_opt=='One-hot'){
    dmy <- dummyVars(" ~ .", data = clinical_train)
    trsf_train <- data.frame(predict(dmy, newdata = clinical_train))
    one_hot_encoder<-dmy
    
    
    
  }else{
    trsf_train<-clinical_train
    one_hot_encoder<-NA
  }
  
  ########## Numerical preprocess / Scale
  if(numerical_col_opt=='MinMax'){
    pre_proc <- preProcess(trsf_train[!colnames(trsf_train) %in% c('SYNTAX.SCORE')], method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_train$SYNTAX.SCORE<-trsf_train$SYNTAX.SCORE
    scalerClinical_obj<-pre_proc
    
    
  }
  
  
  return(list(`Dataset`=clinical_train,`Scaler`=scalerClinical_obj,`One.hot.encoder`=one_hot_encoder))
}
Model_training_whole<-function(clinical_train,model_selection='RF',number_iterations=100,optimization_metric='balAcc'){
  set.seed(123)
  
  input_x_clinical <- dplyr::select(clinical_train, -c("target","SYNTAX.SCORE"))
  input_y_clinical <- clinical_train$target
  input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  
  
  
  
  ctrl <- trainControl(
    method = "cv",
    number = 4,
    search = "random",
    summaryFunction = summary_fnc,
    classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  
  # Control for final model training
  tune_control <- trainControl(
    method = "cv",
    number = 4,
    classProbs = TRUE,
    summaryFunction = summary_fnc,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),          
      method = "rf",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid<- expand.grid(
      mtry=rf_model_clinical$bestTune$mtry
    ))
    
    (model_Clinical <- caret::train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
  }
  
  
  return(list(`ClinicalModel`=model_Clinical))
  
}
Model_training_count_whole<-function(clinical_train,model_selection='RF',number_iterations=100,optimization_metric='MAE',
                                     target_col=c('SYNTAX.SCORE','target')){
  set.seed(123)
  input_x_clinical <- dplyr::select(clinical_train, -target_col)
  input_y_clinical <- clinical_train$SYNTAX.SCORE
  #input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  input_x_snp<-as.matrix(dplyr::select(snp_train, -target_col))
  input_y_snp <- snp_train$SYNTAX.SCORE
  #input_y_snp <- factor(input_y_snp, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  
  ctrl <- trainControl(
    method = "cv",
    number = 4,
    search = "random",
    #summaryFunction = summary_fnc,
    #classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  
  # Control for final model training
  tune_control <- trainControl(
    method = "cv",
    number = 4,
    #classProbs = TRUE,
    #summaryFunction = summary_fnc,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  
  
  if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = input_y_clinical,          
      method = "rf",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid<- expand.grid(
      mtry=rf_model_clinical$bestTune$mtry
      
    ))
    
    (model_Clinical <- caret::train(
      x = input_x_clinical,
      y =input_y_clinical,
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
   }
    
  
  
  return(list(`ClinicalModel`=model_Clinical))
  
}


### Parameters
target_col<-'SYNTAX.SCORE'
factor_col_opt='None'#One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
optimization_metric="bal.Acc" #balAcc,F1Score,Spec.Specificity, Sens.Sensitivity  Accuracy, F2Score, F05Score,"AUC"
model_selection='RF'
number_iterations=50

############## Classifiers model train ###############
snp_train<-snp_set
res_preprocessing<-Preprocessing_fnc_whole(snp_train,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
dataset<-res_preprocessing$Dataset
numerical_scaler<-res_preprocessing$Scaler
one_hot_encoder<-res_preprocessing$One.hot.encoder

res_model<-Model_training_whole(dataset,model_selection=model_selection,number_iterations=number_iterations,optimization_metric = optimization_metric)
count.part_classifier_model<-res_model$ClinicalModel

dataset_high<-dataset[dataset$target==1,]
dataset_low<-dataset[dataset$target==0,]    

### Clinical set only ###
rs_cols<-colnames(dataset_high)[which(grepl("rs",colnames(dataset_high)))]
dataset_clinical<-dataset_high[!colnames(dataset_high) %in% rs_cols]
count.part_regressor_high<-Model_training_count_whole(dataset_clinical,number_iterations = number_iterations,optimization_metric="MAE")

#count.part_regressor_low<-Model_training_count_whole(dataset_low,number_iterations = number_iterations,optimization_metric="MAE")
#count.part_regressor_high<-Model_training_count_whole(dataset_high,number_iterations = number_iterations,optimization_metric="MAE")
    

save(count.part_classifier_model,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/Count_part_classifier.Rda')
save(count.part_regressor_low,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/Count_part_regressor_low.Rda')
save(count.part_regressor_high,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/Count_part_regressor_high.Rda')
save(numerical_scaler,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/Count_part_NumericalScaler.Rda')
save(one_hot_encoder,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/Count_part_OneHotEncoder.Rda')


############ Shapley values ################
library(kernelshap)
library(shapviz)
library(caret)
X <- dataset_low[, -c(8,13), drop = FALSE]  
bg_idx <- sample.int(nrow(X), min(100, nrow(X)))
pred_fun <- function(object, newdata) {
  nd <- as.data.frame(newdata)
  pr <- predict(object, newdata = nd)
  
  v  <- as.numeric(pr[[1]])
  #return(pr)
  unname(v)  # plain numeric vector
}

s <- kernelshap::kernelshap(
  object       = count.part_regressor_low,
  X            = X,
  pred_fun=pred_fun,# rows to explain
  bg_X = X[bg_idx,]
)

save(s,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_low.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_low.Rda')
sv <- shapviz(s)

########################## Visualization with same fonts ##############
name_map <- c(
  `QRS DURATION ms`="QRS DURATION (ms)",
  rs2487928_rs4845625= "rs2487928 x rs4845625",
  rs2487928_rs663129="rs2487928 x rs663129",
  rs2487928_rs6689306="rs2487928 x rs6689306",
  rs4845625_rs663129=" rs4845625 x rs663129",
  rs4845625_rs6689306="rs4845625 x rs6689306",
  rs663129_rs6689306="rs663129 x rs6689306"
)


p <- sv_importance(sv, kind = "bar", max_display = Inf) +
  scale_y_discrete(labels = name_map) +
  
  labs(
    title = "Global SHAP Feature Importance",
    x = "Mean |SHAP value|",
    y = "Predictor"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",                         # no legend in bar plot
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )


p

tiff('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_values_CountPartLowRegressor_bar.jpeg',units="in",width=10,height=6,res=300)
p
dev.off()






########################################################
sv_importance(sv, kind = "bee")

dataset<-snp_set
rs_cols<-colnames(dataset)[which(grepl("rs",colnames(dataset)))]
clinical_cols<-colnames(dataset)[!colnames(dataset) %in% c(rs_cols,'target')]
sv_dependence(sv, v = c('GLU','UREA','CRUSADE SCORE'))
sv_dependence(sv,v=rs_cols)

library(tidyverse)
tiff('Shapley_values_CountPartLowRegressor_bee.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bee",max_display = Inf)
dev.off()

##### Depedence plots --- Clinical ####
tiff('Shapley_values_CountPartLowRegressor_Depedence_Clinical.jpeg',units="in",width=12,height=6,res=300)
clinical_cols<-c('GLU','UREA','CRUSADE SCORE','GALECTIN','DIABETES MELLITUS')
sv_dependence(sv, v = clinical_cols)
dev.off()



### Dependence plots --- SNP ###########
tiff('Shapley_values_CountPartLowRegressor_Depedence_SNP.jpeg',units="in",width=10,height=4,res=300)
snp_cols<-c('rs2487928','rs4845625','rs663129','rs6689306')
sv_dependence(sv, v = snp_cols)
dev.off()


### Dependence plots --- SNP with Interactions ###########
tiff('Shapley_values_CountPartLowRegressor_Depedence_SNPInts.jpeg',units="in",width=10,height=4,res=300)
snp_cols<-c('rs4845625_rs663129','rs2487928_rs4845625','rs4845625','rs6689306',
            'rs663129_rs6689306')
sv_dependence(sv, v = snp_cols)
dev.off()




X <- dataset_high[, -c(8,13), drop = FALSE]  
bg_idx <- sample.int(nrow(X), min(100, nrow(X)))
pred_fun <- function(object, newdata) {
  nd <- as.data.frame(newdata)
  pr <- predict(object, newdata = nd)
  
  v  <- as.numeric(pr[[1]])
  #return(pr)
  unname(v)  # plain numeric vector
}

s <- kernelshap::kernelshap(
  object       = count.part_regressor_high,
  X            = X,
  pred_fun=pred_fun,# rows to explain
  bg_X = X[bg_idx,]
)
save(s,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_high.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_high.Rda')
sv <- shapviz(s)
sv_importance(sv, kind = "bee")

library(tidyverse)
#setwd('C:/Users/30697/Desktop')
tiff('Shapley_values_CountPartHighRegressor_bar.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bar",max_display = Inf)
dev.off()


tiff('Shapley_values_CountPartHighRegressor_bee.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bee",max_display = Inf)
dev.off()



############################### Clinical model --- Regressor high ##########


X <- dataset_clinical[, -c(8,9), drop = FALSE]  
bg_idx <- sample.int(nrow(X), min(100, nrow(X)))
pred_fun <- function(object, newdata) {
  nd <- as.data.frame(newdata)
  pr <- predict(object, newdata = nd)
  
  v  <- as.numeric(pr[[1]])
  #return(pr)
  unname(v)  # plain numeric vector
}

s <- kernelshap::kernelshap(
  object       = count.part_regressor_high,
  X            = X,
  pred_fun=pred_fun,# rows to explain
  bg_X = X[bg_idx,]
)
save(s,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_high_clinical.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_countpart_regressor_high_clinical.Rda')
sv <- shapviz(s)

################## High Regressor ###############
name_map <-name_map <- c(
  `QRS DURATION ms`="QRS DURATION (ms)"
)
p <- sv_importance(sv, kind = "bar", max_display = Inf) +
  scale_y_discrete(labels = name_map) +
  
  labs(
    title = "Global SHAP Feature Importance",
    x = "Mean |SHAP value|",
    y = "Predictor"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",                         # no legend in bar plot
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

tiff('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_values_HighRegressor_bar.jpeg',units="in",width=10,height=6,res=300)
p
dev.off()


####################################################
sv_importance(sv, kind = "bar")

library(tidyverse)
#setwd('C:/Users/30697/Desktop')
tiff('Shapley_values_CountPartHighRegressor_bar_clinical.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bar",max_display = Inf)
dev.off()


tiff('Shapley_values_CountPartHighRegressor_bee.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bee",max_display = Inf)
dev.off()



##### Depedence plots --- Clinical ####
tiff('Shapley_values_CountPartHighRegressor_Depedence_Clinical.jpeg',units="in",width=12,height=6,res=300)
clinical_cols<-c('UREA','CRUSADE SCORE','GALECTIN','DIABETES MELLITUS','QRS DURATION ms')
sv_dependence(sv, v = clinical_cols)
dev.off()



### Dependence plots --- SNP ###########
tiff('Shapley_values_CountPartHighRegressor_Depedence_SNP.jpeg',units="in",width=10,height=4,res=300)
snp_cols<-c('rs2487928','rs4845625','rs663129','rs6689306')
sv_dependence(sv, v = snp_cols)
dev.off()


### Dependence plots --- SNP with Interactions ###########
tiff('Shapley_values_CountPartHighRegressor_Depedence_SNPInts.jpeg',units="in",width=10,height=4,res=300)
snp_cols<-c('rs663129_rs6689306','rs2487928_rs4845625','rs6689306','rs2487928_rs6689306',
            'rs4845625_rs663129')
sv_dependence(sv, v = snp_cols)
dev.off()

###################### LR for interaction and Binary.Risk ######################
snp_cols<-colnames(snp_set)[which(grepl("rs",colnames(snp_set)))]
interaction_cols<-snp_cols[which(grepl("_",snp_cols))]
snp_alone<-setdiff(snp_cols,interaction_cols)
pvals<-c()
imp_interactions<-c()
betas<-list()
for(interaction_col in interaction_cols){
  interaction<-snp_set[,interaction_col]
  target<-snp_set[,"SYNTAX.SCORE"]
  target<-ifelse(target>23,1,0)
  
  df<-data.frame(Interaction=interaction,Target=as.factor(target))
  res<-glm(Target ~ Interaction, family='binomial',data=df)
  res_null<-glm(Target ~ 1, family='binomial',data=df)
  
  res_anova<- anova(res_null,res,test="Chisq")
  
  pvals<-c(pvals,res_anova$`Pr(>Chi)`[2])
  imp_interactions<-c(imp_interactions,interaction_col)
  
  
  coef_table <- coef(summary(res))        # log-odds coefficients
  odds_ratios <- exp(coef(res))
  CI <- exp(confint(res))
  
  level_cols<-strsplit(row.names(CI)[2:length(row.names(CI))],split="Interaction")
  levels<-c()
  for(lev in level_cols){
    levels<-c(levels,lev[2])
  }
  
  info<-c()
  for(num in 2:length(odds_ratios)){
    odds_ratio<-as.numeric(odds_ratios[num])
    CI_num<-as.numeric(CI[num,])
    information<-paste0(round(odds_ratio,3),' (',round(CI_num[1],3),'-',round(CI_num[2],3),')')
    info<-c(info,information)
  }
  temp_list<-list()
  for(num in 1:length(levels)){
    temp_list[levels[num]]<-info[num]
  }
  betas[[interaction_col]]<-temp_list
}
res<-data.frame(Interaction=imp_interactions,P.Value=pvals)
res$P.Value.Adjusted<-p.adjust(res$P.Value,method='fdr')


load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Count_res/SNP_All_levels.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/COunt_res/SNP_All_Mode_selected.Rda')

levels_wanted<-levels_all_snp[names(levels_all_snp) %in% snp_alone]


#####################################################################
library(dplyr)
library(tibble)
library(stringr)
library(purrr)

# -----------------------------
# Function to build results per interaction
# -----------------------------
make_table_for_interaction <- function(interaction_name, betas_list, res_df, levels_wanted) {
  
  # Extract SNP names (e.g. "rs11984041_rs1330344")
  snps <- str_split(interaction_name, "_")[[1]]
  snp1 <- snps[1]
  snp2 <- snps[2]
  
  # Levels for each SNP
  lev1 <- levels_wanted[[snp1]]
  lev2 <- levels_wanted[[snp2]]
  
  # Overall row (interaction p-values)
  interaction_row <- tibble(
    Interaction = paste(snp1, "x", snp2),
    OR = "-",
    `95%CI` = "-",
    pvalue = as.character(signif(res_df$P.Value[res_df$Interaction == interaction_name], 3)),
    pvalue_adjusted = as.character(signif(res_df$P.Value.Adjusted[res_df$Interaction == interaction_name], 3))
  )
  
  # Reference row
  ref_row <- tibble(
    Interaction = paste0(lev1[1], " x ", lev2[1], " (Reference)"),
    OR = "Reference",
    `95%CI` = "-",
    pvalue = "-",
    pvalue_adjusted = "-"
  )
  
  # Other categories (ORs + CI from betas)
  betas_int <- betas_list[[interaction_name]]
  
  other_rows <- tibble(
    Interaction = c(paste0(lev1[2], " x ", lev2[1]),
                    paste0(lev1[1], " x ", lev2[2]),
                    paste0(lev1[2], " x ", lev2[2])),
    OR = str_extract(betas_int, "^[0-9.]+"),
    `95%CI` = str_extract(betas_int, "\\(.*\\)") %>%
      str_remove_all("[()]"),
    pvalue = "-",
    pvalue_adjusted = "-"
  )
  
  bind_rows(interaction_row, ref_row, other_rows)
}

# -----------------------------
# Build the full results table
# -----------------------------
final_table <- map_dfr(
  names(betas), 
  ~ make_table_for_interaction(.x, betas, res, levels_wanted)
)

# View
final_table



library(officer)
library(flextable)
library(dplyr)

# make a copy so we don't overwrite
table_export <- final_table %>%
  mutate(
    # format p-values if numeric-looking
    pvalue = ifelse(pvalue != "-" & !is.na(as.numeric(pvalue)), 
                    formatC(as.numeric(pvalue), format = "e", digits = 2), 
                    pvalue),
    pvalue_adjusted = ifelse(pvalue_adjusted != "-" & !is.na(as.numeric(pvalue_adjusted)), 
                             formatC(as.numeric(pvalue_adjusted), format = "e", digits = 2), 
                             pvalue_adjusted)
  )

# build flextable
ft <- flextable(table_export) %>%
  theme_vanilla() %>%          # clean theme
  autofit() %>%                # auto column sizes
  align(align = "center", part = "all") %>%  # center everything
  bold(part = "header") %>%    # bold headers
  fontsize(size = 11, part = "all") %>%
  fontsize(size = 12, part = "header") %>%
  set_table_properties(layout = "autofit")   # keep table compact

# export to Word
read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = "final_table_publication_BinaryRisk.docx")

library(dplyr)
library(tibble)
library(ggplot2)
importance <- varImp(count.part_classifier_model, scale = TRUE)

# Convert to dataframe
df_imp <- importance$importance
df_imp$Feature <- rownames(df_imp)
df_imp <- df_imp[order(-df_imp$Overall), ][1:17, ] # keep top 28

# Plot with ggplot2
p1<-ggplot(df_imp, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Count part model: Classifier Feature Importance",
       x = "Features",
       y = "Importance (scaled)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text   = element_text(size = 12)
  )



importance <- varImp(count.part_regressor_low$ClinicalModel, scale = TRUE)

# Convert to dataframe
df_imp <- importance$importance
df_imp$Feature <- rownames(df_imp)
df_imp <- df_imp[order(-df_imp$Overall), ][1:17, ] # keep top 28

# Plot with ggplot2
p2<-ggplot(df_imp, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Count part model: Regressor Low Risk Feature Importance",
       x = "Features",
       y = "Importance (scaled)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text   = element_text(size = 12)
  )



importance <- varImp(count.part_regressor_high$ClinicalModel, scale = TRUE)

# Convert to dataframe
df_imp <- importance$importance
df_imp$Feature <- rownames(df_imp)
df_imp <- df_imp[order(-df_imp$Overall), ][1:17, ] # keep top 28

# Plot with ggplot2
p3<-ggplot(df_imp, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Count part model: Regressor High Risk Feature Importance",
       x = "Features",
       y = "Importance (scaled)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text   = element_text(size = 12)
  )





setwd('C:/Users/30697/Desktop')
tiff('FS_CountPartRegressorHigh.jpeg',units="in",width=10,height=6,res=300)
p3
dev.off()


plot_grid(p1, p2,p3, nrow=3, align="h")






######################## Sanity check ########################

data<-snp_set[,c(1:13)]
strings<-''
res<-glm(target ~ Gender + HYPERTENSION + SMOKING +`ATRIAL FIBRILLATION`+ `CHEST PAIN`+              
           GLU+`NEU%`+`EOS%`+`BASO%`+`RATIO1`+RATIO3 +RATIO4+`ATYPICAL ANGINA SYMPTOMS`+
           rs11984041+rs1330344+rs2023938+rs2107595+rs25882+
           rs11984041*rs1330344+rs11984041*rs2023938+rs11984041*rs2107595+rs11984041*rs25882+
           rs1330344*rs2023938+rs1330344*rs2107595+rs1330344*rs25882+
           rs2023938*rs2107595+rs2023938*rs25882+rs2107595*rs25882,
         family='binomial',data=dataset)
r<-anova(res,test="Chisq")

res<-glm(target ~ rs1330344*rs2023938,data=dataset,family='binomial')

