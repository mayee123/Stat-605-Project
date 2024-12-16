library(car)
library(emmeans)
library(broom)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)

#read in data
data <- read.csv("plasmid_data.csv")
data$rank_f <- factor(data$observed_host_range_ncbi_rank)
data$mobility_f_ <- factor(data$predicted_mobility)
data$size_quart <- factor(
  cut(
    data$size, 
    breaks = quantile(data$size, probs = seq(0, 1, 0.25), na.rm = TRUE), 
    include.lowest = TRUE, 
    labels = c("First", "Second", "Third", "Fourth")
)
)

#One way ANOVA
mod=lm(arg_count ~ rank_f, data=data)
summary(mod)
Anova(mod)


#Pairwise comparison
means=emmeans(mod, specs=c("rank_f"))
pairs(means)


#logistic regression
mod.log=glm(range~arg_count+mobility_f_+size_quart,data=data, family=binomial)
summary(mod.log)


#Figure 1
data$observed_host_range_ncbi_rank_f <- factor(data$observed_host_range_ncbi_rank, 
                                               levels = c("genus", "family", "order", "class", "phylum", "multi-phylla"))
data$log_arg_count <- log1p(data$arg_count)
ggplot(data, aes(x = observed_host_range_ncbi_rank_f, y = log_arg_count)) +
  geom_boxplot(fill = "lightblue", color = "black") + # Boxplot with custom colors
  labs(x = "Host Range (NCBI Rank)", y = "Log-transformed Antimicrobial Resistant Gene Count", title = "Box Plot of ARG Count by Host Range Rank") +
  theme_minimal() + # Use a clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Figure 2
coef_df <- tidy(mod.log, conf.int = TRUE)
coef_df$term <- factor(coef_df$term, levels = c("size_quartFourth", "size_quartThird", "size_quartSecond", "mobility_f_non-mobilizable", "mobility_f_mobilizable", "arg_count", "(Intercept)"))
coef_df$significance <- ifelse(coef_df$p.value < 0.05, "Significant", "Not Significant")
ggplot(coef_df, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, color = significance)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  scale_color_manual(values = c("Significant" = "blue4", "Not Significant" = "red")) +
  labs(
    title = "Logistic Regression Coefficients",
    x = "Coefficient Estimate",
    y = "Predictor",
    color = "Significance"
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black")



# Supplemental Figure 1 - 
probabilities <- predict(mod.log, type = "response")
mydata <- data[, 'arg_count', drop=FALSE]
mydata$prob <- probabilities
mydata$logodds <- log(probabilities/(1-probabilities))
ggplot(mydata, aes(x=arg_count, y=logodds))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess", se=FALSE) + 
  theme_bw() + 
  labs(
    title = "Log Odds vs Antimicrobial Resistance Gene Count",
    x = "Antimicrobial Resistance Gene Count",
    y = "Log Odds")
