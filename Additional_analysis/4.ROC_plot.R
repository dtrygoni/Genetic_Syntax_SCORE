roc1 <- roc(test_all,probs_clinical)
roc2 <- roc(test_all,probs_snp)

# --- ROC Curves ---


df_roc1 <- data.frame(FPR = 1 - roc1$specificities,
                      TPR = roc1$sensitivities,
                      Model = "Clinical (AUCROC[95% CI] = 0.731[0.690 - 0.767];N=799)")

df_roc2 <- data.frame(FPR = 1 - roc2$specificities,
                      TPR = roc2$sensitivities,
                      Model = "Clinical+SNPs (AUCROC[95%CI]= 0.760 [0.725 - 0.793];N=799)")

df_roc <- rbind(df_roc1, df_roc2)


specificity_grid <- seq(0, 1, by = 0.01)

ci1 <- pROC::ci.se(
  roc1,
  specificities = specificity_grid,
  boot.n = 2000,
  conf.level = 0.95,
  stratified = TRUE,
  progress = "none"
)

ci2 <- pROC::ci.se(
  roc2,
  specificities = specificity_grid,
  boot.n = 2000,
  conf.level = 0.95,
  stratified = TRUE,
  progress = "none"
)

df_ci1 <- data.frame(
  FPR = 1 - as.numeric(rownames(ci1)),
  Lower = as.numeric(ci1[, 1]),
  Upper = as.numeric(ci1[, 3]),
  Model = "Clinical"
)

df_ci2 <- data.frame(
  FPR = 1 - as.numeric(rownames(ci2)),
  Lower = as.numeric(ci2[, 1]),
  Upper = as.numeric(ci2[, 3]),
  Model = "Clinical+SNPs"
)

df_ci <- rbind(df_ci1, df_ci2)
df_ci <- df_ci[order(df_ci$Model, df_ci$FPR), ]

# Original plot with confidence bands added
p1 <- ggplot(df_roc, aes(x = FPR, y = TPR, color = Model)) +
  geom_ribbon(
    data = df_ci,
    aes(
      x = FPR,
      ymin = Lower,
      ymax = Upper,
      fill = Model
    ),
    inherit.aes = FALSE,
    alpha = 0.20,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "ROC Curves",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.7, 0.2),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

p1
library(tidyverse)
setwd('C:/Users/.../Desktop')
tiff('Figure_2.jpeg',units="in",width=25,height=10,res=300)
plot_grid(p1, ncol=2, align="h")
dev.off()
