### Final analysis based on square-root-transformed absolute errors (AEs) ####

# Load the required R packages
library(reshape2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)
library(performance)
library(emmeans)
library(ggExtra)

# -------------------------------------------------------------------------
# Model specification and comparison
# -------------------------------------------------------------------------
load('./Count_part.Rda')
# Fit the reduced (main-effects) model using maximum likelihood (ML).
# Fixed effects include Model, Risk, and Fold.
# A participant-level random intercept (1 | ID) accounts for the
# paired observations obtained by evaluating both predictive models
# on the same participant.
# The outcome is the square root of the absolute error (AE).

model_reduced <- lmer(sqrt(AE) ~ Model+Risk  + Fold + (1|ID), data = res_df,
                      REML=FALSE)

# Fit the full model, including the Model × Risk interaction.
# This interaction evaluates whether the difference in predictive
# performance between Model A and Model B varies across risk strata.

model_full <- lmer(sqrt(AE) ~ Model*Risk+ Fold + (1|ID), data = res_df,
                   REML=FALSE)

# Compare the nested models using a likelihood ratio (LR) test.
# Both models are fitted using ML to ensure a valid comparison
# of their fixed-effects structures.
#
# A statistically significant result supports retaining the
# Model × Risk interaction.

anova(model_reduced, model_full)

# -------------------------------------------------------------------------
# Final model estimation using restricted maximum likelihood (REML)
# -------------------------------------------------------------------------

# Since the interaction was retained following the LR test,
# refit the full interaction model using REML.
# REML estimation is used for the final estimation of the variance
# components and fixed effects.

model_interaction <- lmer(sqrt(AE) ~ Model*Risk + Fold + (1|ID), data = res_df, 
                          REML=TRUE)

# Display the estimated fixed effects, random-effects variance
# components, and associated inferential statistics.
# The lmerTest package provides Satterthwaite-approximated
# degrees of freedom for the fixed-effects tests.

summary(model_interaction)

# Obtain Type III tests for the fixed effects.
# These tests assess the contribution of each fixed effect,
# including the Model × Risk interaction, conditional on the
# other terms included in the model.
anova(model_interaction)


# -------------------------------------------------------------------------
# Diagnostic assessment of model assumptions
# -------------------------------------------------------------------------

# Extract fitted values and conditional residuals from the final model.
# Conditional residuals are calculated after accounting for the
# estimated participant-specific random intercepts.

diag_df <- data.frame(
  Fitted = fitted(model_interaction),
  Residuals = resid(model_interaction)
)

# Diagnostic plots

p1 <- ggplot(diag_df, aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted values",
    y = "Conditional residuals"
  ) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    axis.title.x = element_text(
      size = 20),
    axis.title.y = element_text(
      size = 20
    ),
    plot.title = element_text(
      size = 24
    ))

p2 <- ggplot(diag_df, aes(sample = Residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Normal Q-Q Plot of Random Intercepts",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    axis.title.x = element_text(
      size = 20),
    axis.title.y = element_text(
      size = 20
    ),
    plot.title = element_text(
      size = 24
    ))


p3 <- ggplot(
  diag_df,
  aes(x = Fitted, y = sqrt(abs(scale(Residuals))))
) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    title = "Scale-Location Plot",
    x = "Fitted values",
    y = expression(sqrt("|Standardized residuals|"))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    axis.title.x = element_text(
      size = 20),
    axis.title.y = element_text(
      size = 20
    ),
    plot.title = element_text(
      size = 24
    ))


re_df <- data.frame(
  RandomIntercepts = ranef(model_interaction)$ID[, 1]
)

p4 <- ggplot(re_df, aes(sample = RandomIntercepts)) +
  stat_qq() +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(
    title = "Normal Q-Q Plot of Random Intercepts",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    # Axis tick labels
    axis.text.x = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    axis.title.x = element_text(
      size = 20),
    axis.title.y = element_text(
      size = 20
    ),
    plot.title = element_text(
      size = 24
    ))

# Combine the four diagnostic plots into a single figure.
(p1 | p3) / (p2 | p4)


# -------------------------------------------------------------------------
# Estimated marginal means and pairwise comparisons
# -------------------------------------------------------------------------

# Calculate estimated marginal means (EMMs) for Model separately
# within each Risk stratum.
#
# This allows the comparison of the two predictive models
# within each SYNTAX score complexity stratum.

emm_change <- emmeans(model_interaction, ~ Model | Risk)

# Perform pairwise comparisons between the predictive models
# within each Risk stratum.
# Tukey adjustment is applied to account for multiple comparisons.


pairs(emm_change, adjust = "tukey")

# Back-transform the estimated marginal means and confidence
# intervals to the original AE scale.
# Since the fitted model uses sqrt(AE) as the response,
# regrid() with response transformation allows the results
# to be presented on the original measurement scale.

emm_original <- regrid(emm_change, transform ="response")

# Convert the back-transformed EMMs into a data frame
# for graphical presentation.

mm_df <- as.data.frame(emm_original)

# Assign descriptive labels to the risk strata and predictive models.

mm_df$Risk <- factor (mm_df$Risk,labels=c("Low complexity stratum", "Intermediate/High complexity stratum"))
mm_df$Model <- factor(mm_df$Model, labels = c("Clinical", "Clinical+SNPs"))

# -------------------------------------------------------------------------
# Graphical presentation of estimated marginal means
# -------------------------------------------------------------------------

# Create a faceted plot displaying the estimated marginal means
# and their corresponding confidence intervals for each predictive
# model, separately for the two SYNTAX score complexity strata.
#
# The vertical axis represents the estimated marginal AE on
# the original measurement scale.

EMM_stratum <- ggplot(mm_df, aes(x = Model, y = response)) +
  facet_wrap(.~ Risk)+
  
  geom_point(position = position_dodge(0.3), size = 3) +
  #geom_line(position = position_dodge(0.3), linewidth=1.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2,linewidth=1.5, position = position_dodge(0.3)) +
  theme_bw()+
  removeGrid(x=TRUE, y=TRUE)+
  theme(legend.position= c("bottom"),
        #legend.justification = c("bottom"),
        legend.title = element_text(size=24, face="bold"),
        legend.text = element_text(size=24, face="bold"),
        strip.text.x = element_text(size = 22, hjust = 0.5, vjust = 0.5, face = 'bold'),
        axis.text.x = element_text(face="bold", color="black", size=22),
        axis.text.y = element_text(face="bold", color="black", size=22),
        axis.title.x = element_text(face="bold", color="black", size=24),
        axis.title.y = element_text(face="bold", color="black", size=24),
        plot.title = element_text(face="bold", color="black", size=22, hjust = 0.5))+
  labs(y = "Estimated Marginal Means for AE", x = "Model")


# The same code can also be applied to the end-to-end performance evaluation.
# In this case, the Risk factor should include all three SYNTAX score
# complexity strata (Zero SYNTAX score, Low complexity, and Intermediate/High complexity),
# and the corresponding labels should be adjusted accordingly.