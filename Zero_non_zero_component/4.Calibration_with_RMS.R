library(rms)
r<-val.prob(probs_snp, test_all,pl=TRUE)
r



#slope 0.789
#intercept -0.049


#slope 0.945
#intercept 0.033



library(CalibrationCurves)
c<-valProbggplot(p=probs_snp,y=test_all,smooth="rcs",CL.smooth = FALSE)
c

# 1. Δημιουργία τυχαίων δεδομένων (Παράδειγμα)
set.seed(123)
y <- test_all
pred_model1 <- probs_clinical
pred_model2 <- probs_snp

# 2. Ένωση των δεδομένων σε ένα "tidy" data frame
df1 <- data.frame(y = y, pred = pred_model1, Model = "Clinical")
df2 <- data.frame(y = y, pred = pred_model2, Model = "Clinical+SNPs")
plot_data <- rbind(df1, df2)






library(dplyr)
library(ggplot2)

# Ομαδοποίηση ανά 10% (deciles) και υπολογισμός μέσων όρων
binned_data <- plot_data %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(Model, bin) %>%
  summarize(
    mean_pred = mean(pred),
    mean_obs = mean(y),
    .groups = 'drop'
  )

# Σχεδίαση binned γραμμών
p<-ggplot(binned_data, aes(x = mean_pred, y = mean_obs, color = Model, group = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray") +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +
  xlim(0, 1) + ylim(0, 1) +
  labs(
    x = "Observed proportion",
    y = "Predicted probability"
  ) +
  theme(
    legend.position = c(0.7, 0.2),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  theme_bw()

library(tidyverse)
setwd('C:/Users/.../Desktop')
tiff('Calibration.jpeg',units="in",width=10,height=5,res=300)
p
dev.off()
