


library(dcurves)
#install.packages("dcurves")


dca_data <- data.frame(
  # Convert factor outcome back to a 0/1 numeric variable for dcurves
  outcome_numeric = test_all,
  Clinical = probs_clinical,
  `Clinical+SNPs` = probs_snp
)
colnames(dca_data)<-c("Target","Clinical","Clinical+SNPs")

dca_results <- dca(
  Target ~ Clinical+`Clinical+SNPs`,
  data = dca_data,
  thresholds = seq(0.10, 0.80, by = 0.01) # Define threshold range of interest (e.g., 0% to 50%)
)
dca_results
?dca()


p<-dca(
  Target ~ Clinical+`Clinical+SNPs`,
  data = dca_data,
  thresholds = seq(0.05, 0.80, by = 0.01) # Define threshold range of interest (e.g., 0% to 50%)
) %>%
  plot(smooth=TRUE)+
  ggplot2::labs(
    x="Threshold Probability",
    y="Net Benetif"
  )+
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
p
library(tidyverse)
setwd('C:/Users/.../Desktop')
tiff('DCA_plot.jpeg',units="in",width=10,height=5,res=300)
p
dev.off()
