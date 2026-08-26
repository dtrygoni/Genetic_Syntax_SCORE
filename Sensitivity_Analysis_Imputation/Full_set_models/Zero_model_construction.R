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
library(reticulate)
set.seed(123)
xgb.set.config(verbosity = 0)

model_arr<-c()
type_arr<-c()
acc_arr<-c()
prec_arr<-c()
rec_arr<-c()
f1_arr<-c()
aucpr_arr<-c()
aucroc_arr<-c()
fold_arr<-c()
final_models<-list()
###### load
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/SNP_dataset.Rda')


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
    balAcc = bal_acc,
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
    pre_proc <- preProcess(trsf_train, method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)

    scalerClinical_obj<-pre_proc
    
    
   }else{
    clinical_train<-trsf_train
    scalerClinical_obj<-NA
  }
  
  
  return(list(`Dataset`=clinical_train,`Scaler`=scalerClinical_obj,`One.hot.encoder`=one_hot_encoder))
}
Model_training_whole<-function(clinical_train,model_selection='RF',number_iterations=100,optimization_metric='balAcc'){
  set.seed(123)
  
  input_x_clinical <-dplyr::select(clinical_train, -c("target"))
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

 ### Parameters
genetic_modelling<-'Codominant'
HWE_threshold<-0.05
MAF_threshold<-0.05
clinical_cols<-1:47
LR_threshold<-0.10
Boruta_threshold<-0.05
strategy_boruta<-'Seperate' ## All second experiments
target_col<-'Binary'
snp_cols<-48:ncol(snp_set)
factor_col_opt='None'#'One-hot' #One-hot, L-1_encoding
numerical_col_opt='None' #None, MinMax
optimization_metric="AUC" #balAcc,F1Score,Spec.Specificity, Sens.Sensitivity  Accuracy, F2Score, F05Score,"AUC"
model_selection='RF'
number_iterations=50






snp_train<-snp_set
res_preprocessing<-Preprocessing_fnc_whole(snp_train,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
dataset<-res_preprocessing$Dataset
numerical_scaler<-res_preprocessing$Scaler
one_hot_encoder<-res_preprocessing$One.hot.encoder


res_model<-Model_training_whole(dataset,model_selection=model_selection,number_iterations=number_iterations,optimization_metric = optimization_metric)
zero.part_model<-res_model$ClinicalModel

save(zero.part_model,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/Zero_part_model.Rda')
save(numerical_scaler,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/Zero_part_NumericalScaler.Rda')
save(one_hot_encoder,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/Zero_part_OneHotEncoder.Rda')


################################################## SHAPLEY VALUES ##############################################
library(kernelshap)
library(shapviz)
library(caret)
# 1) Pick the positive class column name robustly
X <- dataset[, -19, drop = FALSE]     # predictors only, keep as data.frame

# 1) Detect the positive-class column once
detect_pos_col <- function(model, newdata) {
  p1 <- predict(model, newdata = newdata[1, , drop = FALSE], type = "prob")
  p1 <- as.data.frame(p1)
  # Prefer a column literally named "pos" / "1" etc., otherwise take the last
  candidates <- intersect(colnames(p1), c("pos","1","yes","Y","TRUE","True"))
  if (length(candidates)) candidates[1] else tail(colnames(p1), 1)
}
pos_col <- detect_pos_col(zero.part_model, X)

# 2) Prediction function MUST return numeric
pred_fun <- function(object, newdata) {
  nd <- as.data.frame(newdata)
  pr <- predict(object, newdata = nd, type = "prob")
  pr <- as.data.frame(pr)
  v  <- pr[[pos_col]]
  v  <- as.numeric(v)
  if (!is.numeric(v) || anyNA(v)) stop("pred_fun did not return numeric probs")
  unname(v)  # plain numeric vector
}

# 3) Sanity checks (run these)
stopifnot(is.numeric(pred_fun(zero.part_model, X[1, , drop = FALSE])),
          length(pred_fun(zero.part_model, X[1, , drop = FALSE])) == 1L)
stopifnot(is.numeric(pred_fun(zero.part_model, X)),
          length(pred_fun(zero.part_model, X)) == nrow(X))

# 4) Choose a smaller background for speed/robustness
set.seed(42)
bg_idx <- sample.int(nrow(X), min(100, nrow(X)))
bg_idx2 <- sample.int(nrow(X), min(150, nrow(X)))

# 5) KernelSHAP
set.seed(1)
s <- kernelshap::kernelshap(
  object       = zero.part_model,
  X            = X,                         # rows to explain
  pred_fun     = pred_fun,
  bg_X = X[bg_idx,]
)


save(s,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_zeropart.Rda')

load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_zeropart.Rda')
sv <- shapviz(s)
sv_importance(sv,kind="bar",max_display = Inf)

################### Plot modification ####################
# Create mapping dictionary
name_map <- c(
  Gender="GENDER",
  HYPERTENSION="HYPERTENSION",
  SMOKING="SMOKING",
  rs11984041_rs1330344= "rs11984041 x rs1330344",
  rs11984041_rs2023938="rs11984041 x rs2023938",
  rs11984041_rs2107595="rs11984041 x rs2107595",
  rs11984041_rs25882="rs11984041 x rs25882",
  rs1330344_rs2023938="rs1330344 x rs2023938",
  rs1330344_rs2107595="rs1330344 x rs2107595",
  rs1330344_rs25882= "rs1330344 x rs25882" ,      
  rs2023938_rs2107595="rs2023938 x rs2107595",
  rs2023938_rs25882="rs2023938 x rs25882",
  rs2107595_rs25882="rs2107595 x rs25882"       
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

tiff('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_values_ZeroPart_bar.jpeg',units="in",width=10,height=6,res=300)
p
dev.off()


####################################################

sv_importance(sv, kind = "bee")
sv_importance(sv, kind = "both")

#sv_summary(sv)
sv_waterfall(sv, row_id = 100)

library(tidyverse)
#setwd('C:/Users/30697/Desktop')
tiff('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Shapley_values_ZeroPart_bee.jpeg',units="in",width=10,height=4,res=300)
sv_importance(sv, kind = "bar",max_display= Inf)
dev.off()





dataset<-snp_set
rs_cols<-colnames(dataset)[which(grepl("rs",colnames(dataset)))]
clinical_cols<-colnames(dataset)[!colnames(dataset) %in% c(rs_cols,'target')]


clinical_cols<-c('RATIO1','GLU','NEU%','CHEST PAIN','RATIO3','RATIO4')
tiff('Shapley_values_ZeroPart_Depedence_Clinical.jpeg',units="in",width=10,height=4,res=300)
sv_dependence(sv, v = clinical_cols)
dev.off()


snp_cols<-c("rs11984041_rs2107595", "rs2023938_rs2107595",  "rs1330344_rs2023938",  "rs2107595_rs25882", "rs1330344_rs25882",
            "rs1330344_rs2107595")
tiff('Shapley_values_ZeroPart_Depedence_SNPsInt.jpeg',units="in",width=10,height=6,res=300)
sv_dependence(sv, v = snp_cols)
dev.off()

########## Second approach 
# install.packages(c("fastshap", "shapviz"))
library(fastshap)
library(shapviz)

# Predictors only; keep as data.frame (factors + numerics are fine)
X <- dataset[, -19, drop = FALSE]   # adjust if target is not the 19th column

# caret model must have class probabilities enabled when trained:
# trControl = trainControl(classProbs = TRUE, ...)

# Detect the positive-class column (e.g., "pos")
detect_pos_col <- function(model, newdata) {
  p1 <- predict(model, newdata = newdata[1, , drop = FALSE], type = "prob")
  nm <- colnames(as.data.frame(p1))
  cand <- intersect(nm, c("pos","1","yes","Y","TRUE","True"))
  if (length(cand)) cand[1] else tail(nm, 1)
}
pos_col <- detect_pos_col(zero.part_model, X)

# Numeric probability wrapper for fastshap
pred_prob <- function(object, newdata) {
  as.numeric(as.data.frame(predict(object, newdata = as.data.frame(newdata), type = "prob"))[[pos_col]])
}

# Choose a compact, representative background (reference distribution)
set.seed(42)
bg_n  <- min(150, nrow(X))           # 50–200 usually good
bg_ix <- sample.int(nrow(X), bg_n)
bg_X  <- X[bg_ix, , drop = FALSE]

# Compute SHAP for all rows (or subset if you only want global insights)
set.seed(1)
phi <- fastshap::explain(
  object       = zero.part_model,
  X            = bg_X,               # background set
  pred_wrapper = pred_prob,          # returns numeric probs
  newdata      = X,                  # rows to explain (all)
  nsim         = 120                 # 80–200 typical; increase if noisy
)

# Visualize
sv <- shapviz::shapviz(phi, X = X)
sv_importance(sv)             # global: mean |SHAP| per feature
sv_summary(sv)                # global: beeswarm (distribution of effects)
sv_dependence(sv, "GLU")      # partial effect, colored by strongest interaction
sv_waterfall(sv, row_index=1) # local: a single prediction explanation









############## ROC evaluation #########
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/Zero_part_model.Rda')
dd<-dataset[!colnames(dataset) %in% c('target')]
probs<-predict(zero.part_model$finalModel,newdata=dataset[1,])


###################### LR for interaction and Binary.SS ######################
snp_cols<-colnames(snp_set)[which(grepl("rs",colnames(snp_set)))]
interaction_cols<-snp_cols[which(grepl("_",snp_cols))]
snp_alone<-setdiff(snp_cols,interaction_cols)
pvals<-c()
imp_interactions<-c()
betas<-list()
for(interaction_col in interaction_cols){
  interaction<-snp_set[,interaction_col]
  target<-snp_set[,"target"]
  
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


load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/SNP_All_levels.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Final_model/Binary_res/SNP_All_Mode_selected.Rda')

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
  print(target = "final_table_publication.docx")
########################################################

# Print importance table
print(importance)

# Plot importance
plot(importance, top = 28)  # sho

importance <- varImp(zero.part_model, scale = TRUE)


library(dplyr)
library(tibble)
library(ggplot2)
# Convert to dataframe
df_imp <- importance$importance
df_imp$Feature <- rownames(df_imp)
df_imp <- df_imp[order(-df_imp$Overall), ][1:28, ] # keep top 28

# Plot with ggplot2
p<-ggplot(df_imp, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Zero part model Feature Importance",
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
tiff('FS_ZeroPart.jpeg',units="in",width=10,height=6,res=300)
p
dev.off()









##############################  ANOVA glm
data<-dataset[,c(1:19)]
res<-glm(target ~ Gender + HYPERTENSION + SMOKING +`ATRIAL FIBRILLATION`+ `CHEST PAIN`+              
          GLU+`NEU%`+`EOS%`+`BASO%`+`RATIO1`+RATIO3 +RATIO4+`ATYPICAL ANGINA SYMPTOMS`+
         rs11984041+rs1330344+rs2023938+rs2107595+rs25882+
           rs11984041*rs1330344+rs11984041*rs2023938+rs11984041*rs2107595+rs11984041*rs25882+
           rs1330344*rs2023938+rs1330344*rs2107595+rs1330344*rs25882+
           rs2023938*rs2107595+rs2023938*rs25882+rs2107595*rs25882,
         family='binomial',data=dataset)
r<-anova(res,test="Chisq")
res<-glm(target ~ rs1330344*rs2023938,data=dataset,family='binomial')
anova(res,test="Chisq")
dataset<-snp_set

library(officer)
library(flextable)


# Turn into data.frame
anova_df <- as.data.frame(r)
anova_df_round<-round(anova_df[1:4],3)
anova_df_round$`Pr(>Chi)`<-round(anova_df$`Pr(>Chi)`,4)
anova_df_round$Predictor<-row.names(anova_df)
anova_df<-anova_df_round[,c(6,1:5)]
# Create a Word document
doc <- read_docx() %>%
  body_add_flextable(flextable(anova_df))

# Save
print(doc, target = "anova_table_ModelB_Zero.docx")


rs_col<-colnames(dataset)[which(grepl("rs",colnames(dataset)))]
clinical_cols<-colnames(dataset)[!colnames(dataset) %in% rs_col]
dataset_clinical<-dataset[,clinical_cols]
res_clinical<-glm(target~.,family='binomial',data=dataset_clinical)
r<-anova(res_clinical,test="Chisq")
# Turn into data.frame
anova_df <- as.data.frame(r)
anova_df_round<-round(anova_df[1:4],3)
anova_df_round$`Pr(>Chi)`<-round(anova_df$`Pr(>Chi)`,4)
anova_df_round$Predictor<-row.names(anova_df)
anova_df<-anova_df_round[,c(6,1:5)]
# Create a Word document
doc <- read_docx() %>%
  body_add_flextable(flextable(anova_df))

# Save
print(doc, target = "anova_table_ModelA_Zero.docx")