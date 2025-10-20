library(pROC)
library(PRROC)
library(ggplot2)
library(cowplot)

roc1 <- roc(test_all,probs_clinical)
roc2 <- roc(test_all,probs_snp)

# --- ROC Curves ---
#roc1 <- roc(labels, pred_clinical)
#roc2 <- roc(labels, pred_clinical_snps)

df_roc1 <- data.frame(FPR = 1 - roc1$specificities,
                      TPR = roc1$sensitivities,
                      Model = "Clinical")

df_roc2 <- data.frame(FPR = 1 - roc2$specificities,
                      TPR = roc2$sensitivities,
                      Model = "Clinical+SNPs")

df_roc <- rbind(df_roc1, df_roc2)


# --- PR Curves ---
pr1 <- PRROC::pr.curve(probs_clinical,weights.class0=test_all,curve=TRUE,max.compute = TRUE,min.compute = TRUE)

pr2 <- PRROC::pr.curve(probs_snp,weights.class0=test_all,curve=TRUE)

df_pr1 <- data.frame(Recall = pr1$curve[,1],
                     Precision = pr1$curve[,2],
                     Model = "Clinical")

df_pr2 <- data.frame(Recall = pr2$curve[,1],
                     Precision = pr2$curve[,2],
                     Model = "Clinical+SNPs")

df_pr <- rbind(df_pr1, df_pr2)

# --- ROC Curves ---
p1 <- ggplot(df_roc, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title = "ROC Curves",
       x = "False Positive Rate (1-Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(
    legend.position = c(0.7, 0.2),   # inside plot
    legend.background = element_rect(fill="white", color="black"),
    legend.text = element_text(size=14, face="bold"),   # bigger legend text
    legend.title = element_blank(),                     # remove legend title
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=12, face="bold"),
    plot.title = element_text(size=16, face="bold", hjust=0.5)
  )

# --- PR Curves ---
p2 <- ggplot(df_pr, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(size=1) +
  geom_hline(yintercept = 0.7, linetype="dashed") +
  labs(title = "PR Curves",
       x = "Recall",
       y = "Precision") +
  theme_bw() +
  theme(
    legend.position = c(0.7, 0.2),
    legend.background = element_rect(fill="white", color="black"),
    legend.text = element_text(size=14, face="bold"),
    legend.title = element_blank(),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=12, face="bold"),
    plot.title = element_text(size=16, face="bold", hjust=0.5)
  )
p2 <- p2 + ylim(0.0, 1)
p2 <- p2 + xlim(0,1)
# --- Combine plots side by side ---



library(tidyverse)
#setwd('C:/Users/30697/Desktop')
tiff('ROC_PR_corrected.jpeg',units="in",width=10,height=4,res=300)
plot_grid(p1, p2, ncol=2, align="h")
dev.off()