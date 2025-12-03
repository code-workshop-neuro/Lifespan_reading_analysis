library(pwr)
library(lme4)
library(mgcv)
library(dplyr)
library(ggseg)
library(ppcor)
library(tidyr)
library(readxl)
library(scales)
library(cowplot)
library(emmeans)
library(ggplot2)
library(stringr)
library(writexl)
library(R.matlab)
library(gridExtra)
library(patchwork)
library(paletteer)
library(tidyverse)
library(colorspace)

########################################################################################################################################################################################################
df1 <- read_excel(".../bhv_157_long_updated_group.xls",sheet = "Sheet2")
df_area <- read_excel(".../merged_data_aparc_area.xlsx")
df_thickness <- read_excel(".../merged_data_aparc_thickness.xlsx")
df_volume <- read_excel(".../merged_data_aparc_volume.xlsx")
df_lgi <- read_excel(".../merged_data_aparc_lgi.xlsx")

merged_df_area <- merge(df1, df_area, by = "id")  
cleaned_df_area <- merged_df_area[merged_df_area$incl == 1,]

merged_df_thickness <- merge(df1, df_thickness, by = "id") 
cleaned_df_thickness <- merged_df_thickness[merged_df_thickness$incl == 1,]

merged_df_volume <- merge(df1, df_volume, by = "id")  
cleaned_df_volume <- merged_df_volume[merged_df_volume$incl == 1,]

merged_df_lgi <- merge(df1, df_lgi, by = "id")  
cleaned_df_lgi <- merged_df_lgi[merged_df_lgi$incl == 1,]

columns_to_process <- names(cleaned_df_thickness)[25:92]

df <- read_excel(".../euler_summary.xlsx")
df_long <- df %>%
  pivot_longer(cols = c(lh_defect, rh_defect),
               names_to = "hemisphere",
               values_to = "defect")
df <- df %>% mutate(total = lh_defect + rh_defect)
threshold <- mean(df$total) + sd(df$total)

p <- ggplot(df_long, aes(x = defect, y = factor(subject), fill = hemisphere)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 0.5) +
  scale_fill_manual(
    values = c("lh_defect" = "#1f77b4", "rh_defect" = "#ff7f0e"),
    labels = c("lh_defect" = "Left hemisphere",
               "rh_defect" = "Right hemisphere"),
    name = "Hemisphere"
  ) +
  labs(
    x = "Euler Number",
    y = "Subject ID"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

########################################################################################################################################################################################################
coefficients_list <- data.frame(
  column = integer(),
  label = character(),
  t_value_area = numeric(),
  p_value_area = numeric(),
  
  t_value_thickness = numeric(),
  p_value_thickness = numeric(),
  
  t_value_volume = numeric(),
  p_value_volume = numeric(),
  
  t_value_lgi = numeric(),
  p_value_lgi = numeric()
)

for(i in 1:length(columns_to_process)){
  col_name <- columns_to_process[i];
  print(paste("Processing column:", col_name));
  
  model_area <- lm(paste(col_name , " ~ age_at_test + sex + eTIV + euler_number"), data = cleaned_df_area);
  model_thickness <- lm(paste(col_name , " ~ age_at_test + sex + eTIV + euler_number"), data = cleaned_df_thickness);
  model_volume <- lm(paste(col_name , " ~ age_at_test + sex + eTIV + euler_number"), data = cleaned_df_volume);
  model_lgi <- lm(paste(col_name , " ~ age_at_test + sex + eTIV + euler_number"), data = cleaned_df_lgi);

  model_area_summary <- summary(model_area);
  model_thickness_summary <- summary(model_thickness);
  model_volume_summary <- summary(model_volume);
  model_lgi_summary <- summary(model_lgi);

  t_values_area <- coef(model_area_summary)["age_at_test", "t value"]
  p_values_area <- coef(model_area_summary)["age_at_test", "Pr(>|t|)"]
  
  t_values_thickness <- coef(model_thickness_summary)["age_at_test", "t value"]
  p_values_thickness <- coef(model_thickness_summary)["age_at_test", "Pr(>|t|)"]
  
  t_values_volume <- coef(model_volume_summary)["age_at_test", "t value"]
  p_values_volume <- coef(model_volume_summary)["age_at_test", "Pr(>|t|)"]
  
  t_values_lgi <- coef(model_lgi_summary)["age_at_test", "t value"]
  p_values_lgi <- coef(model_lgi_summary)["age_at_test", "Pr(>|t|)"]
  
  coefficients_list <- rbind(coefficients_list,data.frame(
    column = i,
    label = col_name,
    t_value_area = t_values_area,
    p_value_area = p_values_area,
    
    t_value_thickness = t_values_thickness,
    p_value_thickness = p_values_thickness,
    
    t_value_volume = t_values_volume,
    p_value_volume = p_values_volume,
    
    t_value_lgi = t_values_lgi,
    p_value_lgi = p_values_lgi
  ))
}

coefficients_list$adjusted_p_value_area <- p.adjust(coefficients_list$p_value_area, method = "fdr")
coefficients_list$corrected_area <- ifelse(coefficients_list$adjusted_p_value_area*4 < 0.05, coefficients_list$t_value_area, NA)

coefficients_list$adjusted_p_value_thickness <- p.adjust(coefficients_list$p_value_thickness, method = "fdr")
coefficients_list$corrected_thickness <- ifelse(coefficients_list$adjusted_p_value_thickness*4 < 0.05, coefficients_list$t_value_thickness, NA)

coefficients_list$adjusted_p_value_volume <- p.adjust(coefficients_list$p_value_volume, method = "fdr")
coefficients_list$corrected_volume <- ifelse(coefficients_list$adjusted_p_value_volume*4 < 0.05, coefficients_list$t_value_volume, NA)

coefficients_list$adjusted_p_value_lgi <- p.adjust(coefficients_list$p_value_lgi, method = "fdr")
coefficients_list$corrected_lgi <- ifelse(coefficients_list$adjusted_p_value_lgi*4 < 0.05, coefficients_list$t_value_lgi, NA)

colors.morph = paletteer_d('rcartocolor::TealRose')
lim = 15
lim_max = 15

plot_area <- ggseg(coefficients_list, mapping = aes(fill = corrected_area),color='black',size=0.1) +
  labs(title = "Area", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

plot_thickness <- ggseg(coefficients_list, mapping = aes(fill = corrected_thickness),color='black',size=0.1) +
  labs(title = "Thickness", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

plot_volume <- ggseg(coefficients_list, mapping = aes(fill = corrected_volume),color='black',size=0.1) +
  labs(title = "Volume", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

plot_lgi <- ggseg(coefficients_list, mapping = aes(fill = corrected_lgi),color='black',size=0.1) +
  labs(title = "LGI", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())#终于去掉讨厌的“hemisphere”字样！！

plot_area_with_legend <- plot_area + theme(legend.position = 'bottom')
combined_plot <- (plot_area_with_legend | plot_thickness) /
  (plot_volume | plot_lgi) +
  plot_layout(ncol = 1, nrow = 2, guides = 'collect') &
  theme(legend.position = 'right')
print(combined_plot)


colors.morph = paletteer_d('rcartocolor::TealRose')
process_and_plot <- function(df, group_label) {
  count_linear <- 0
  count_moderate <- 0
  count_high <- 0
  linear_like <- list()
  moderate_non_linearity <- list()
  high_non_linearity <- list()
  all_derivatives <- list()   
  target_ages <- c(11, 22, 55) 
  target_derivatives <- list() 
  for (i in seq_along(indices)) {
    index <- indices[i]
    col_name <- columns_to_process[index]
    print(paste("Processing column:", col_name))
    # model <- gam(as.formula(paste(col_name, " ~ s(age_at_test) + sex")), data = df, method = "REML")
    model <- gam(as.formula(paste(col_name, " ~ s(age_at_test) + sex + eTIV + euler_number")), data = df, method = "REML")
    edf_value <- summary(model)$s.table["s(age_at_test)", "edf"]
    age_values <- seq(min(df$age_at_test), max(df$age_at_test), length.out = 100)
    new_data <- data.frame(age_at_test = age_values, sex = unique(df$sex)[1])
    smooth_values <- predict(model, newdata = new_data, type = "terms")[, "s(age_at_test)"]
    lpmatrix <- predict(model, newdata = new_data, type = "lpmatrix")
    smooth_cols <- grep("s\\(age_at_test\\)", colnames(lpmatrix)) 
    deriv_values <- lpmatrix[, smooth_cols] %*% coef(model)[smooth_cols]
    all_derivatives[[i]] <- list(age_values = age_values, deriv_values = deriv_values)
    specific_derivatives <- sapply(target_ages, function(target_age) {
      closest_index <- which.min(abs(age_values - target_age))
      deriv_values[closest_index]
    })
    target_derivatives[[col_name]] <- data.frame(
      Age = target_ages,
      Derivative = specific_derivatives
    )
    if (edf_value < 1.5) {
      linear_like[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
      count_linear <- count_linear + 1
    } else if (edf_value >= 1.5 & edf_value <= 3) {
      moderate_non_linearity[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
      count_moderate <- count_moderate + 1
    } else {
      high_non_linearity[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
      count_high <- count_high + 1
    }
  }
  print(paste("Linear-like group count:", count_linear))
  print(paste("Moderate non-linearity group count:", count_moderate))
  print(paste("High non-linearity group count:", count_high))
  
  all_smooth_values <- unlist(c(sapply(linear_like, function(x) x$smooth_values),
                                sapply(moderate_non_linearity, function(x) x$smooth_values),
                                sapply(high_non_linearity, function(x) x$smooth_values)))
  y_range <- range(all_smooth_values)
  all_deriv_values <- unlist(sapply(all_derivatives, function(x) x$deriv_values))
  y_range_deriv <- range(all_deriv_values)
  plot_group <- function(group_data, title) {
    plot(NULL, xlim = range(df$age_at_test), ylim = y_range, 
         xlab = "Age", ylab = "s(Age)", main = title)
    for (i in seq_along(group_data)) {
      if (!is.null(group_data[[i]])) {
        lines(group_data[[i]]$age_values, group_data[[i]]$smooth_values, 
              col = colors[i %% length(colors)], lty = 1, lwd = 2)
      }
    }
  }
  plot_group(linear_like, paste(group_label, "Linear-like (edf < 1.5)"))
  plot_group(moderate_non_linearity, paste(group_label, "Moderate non-linearity (1.5 ≤ edf ≤ 3)"))
  plot_group(high_non_linearity, paste(group_label, "High non-linearity (edf > 3)"))
}
par(mfrow = c(4, 3))
process_and_plot(cleaned_df_area, group_label = "Area")
process_and_plot(cleaned_df_thickness, group_label = "Thickness")
process_and_plot(cleaned_df_volume, group_label = "Volume")
process_and_plot(cleaned_df_lgi, group_label = "Local Gyrification Index")
par(mfrow = c(1, 1))
target_ages <- c(11, 28, 65)
target_derivatives <- list()
traget_feature = cleaned_df_lgi

for (i in seq_along(indices)) {
  index <- indices[i]
  col_name <- columns_to_process[index]
  print(paste("Processing column:", col_name))
  # model_area_age <- gam(as.formula(paste(col_name, " ~ s(age_at_test) + sex")), data = traget_feature, method = "REML")
  model_area_age <- gam(as.formula(paste(col_name, " ~ s(age_at_test) + sex + eTIV + euler_number")), data = traget_feature, method = "REML")
  edf_value <- summary(model_area_age)$s.table["s(age_at_test)", "edf"]  
  age_values <- seq(min(traget_feature$age_at_test), max(traget_feature$age_at_test), length.out = 100)  
  new_data <- data.frame(age_at_test = age_values, sex = unique(traget_feature$sex)[1]) 
  smooth_values <- predict(model_area_age, newdata = new_data, type = "terms")[, "s(age_at_test)"]
  deriv_values <- diff(smooth_values) / diff(age_values)
  deriv_values <- c(NA, deriv_values)
  all_derivatives[[i]] <- list(age_values = age_values, deriv_values = deriv_values)
  specific_derivatives <- sapply(target_ages, function(target_age) {
    closest_index <- which.min(abs(age_values - target_age))
    deriv_values[closest_index]
  })
  target_derivatives[[col_name]] <- data.frame(
    Age = target_ages,
    Derivative = specific_derivatives
  )
  if (edf_value < 1.5) {
    linear_like[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
    count_linear <- count_linear + 1
  } else if (edf_value >= 1.5 & edf_value <= 3) {
    moderate_non_linearity[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
    count_moderate <- count_moderate + 1  
  } else {
    high_non_linearity[[i]] <- list(age_values = age_values, smooth_values = smooth_values)
    count_high <- count_high + 1  
  }
}

plot_age_derivative <- function(age_to_plot, target_derivatives, colors.morph) {
  derivatives_at_age <- data.frame(
    label = names(target_derivatives),
    derivative = sapply(target_derivatives, function(x) x$Derivative[x$Age == age_to_plot])
  )
  plot <- ggseg(derivatives_at_age, mapping = aes(fill = derivative), color = "black",size=0.1) +
    labs(title = paste(age_to_plot, " years old"), fill = " ") +
    scale_fill_gradientn(colors = colors.morph, limits = c(-0.10, 0.10), oob = squish, na.value = "white") +
    theme_minimal(base_family = "Gill Sans") +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 10), axis.title = element_blank()) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  return(plot)
}

colors.morph = paletteer_d('rcartocolor::TealRose')
plot_1 <- plot_age_derivative(11, target_derivatives, colors.morph)
plot_2 <- plot_age_derivative(28, target_derivatives, colors.morph)
plot_3 <- plot_age_derivative(65, target_derivatives, colors.morph)
plot_1_no_legend <- plot_1 + theme(legend.position = 'none')
plot_2_no_legend <- plot_2 + theme(legend.position = 'none')
plot_3_no_legend <- plot_3 + theme(legend.position = 'none')

get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}
legend <- get_legend(plot_3)
combined_plot <- (
  plot_1_no_legend / 
    plot_2_no_legend /
    plot_3  
) + plot_layout(ncol = 1)
print(combined_plot)


########################################################################################################################################################################################################
# coefficients_list <- data.frame(
#   column = integer(),
#   label = character(),
#   t_value_area = numeric(),
#   p_value_area = numeric(),
#   t_value_thickness = numeric(),
#   p_value_thickness = numeric(),
#   t_value_volume = numeric(),
#   p_value_volume = numeric(),
#   t_value_lgi = numeric(),
#   p_value_lgi = numeric()
# )
# language = "list"#"sentence"


behav_base <- cleaned_df_thickness[, c("id", "list", "age_at_test", "sex", "nonverbal","eTIV","euler_number")]
model_read <- lm(list ~ age_at_test , data = behav_base)
behav_base$list_resid <- resid(model_read)
cleaned_df_area      <- merge(cleaned_df_area,      behav_base[, c("id", "list_resid")], by = "id")
cleaned_df_thickness <- merge(cleaned_df_thickness, behav_base[, c("id", "list_resid")], by = "id")
cleaned_df_volume    <- merge(cleaned_df_volume,    behav_base[, c("id", "list_resid")], by = "id")
cleaned_df_lgi       <- merge(cleaned_df_lgi,       behav_base[, c("id", "list_resid")], by = "id")

coefficients_list <- data.frame(
  column = integer(),
  label = character(),
  t_value_area = numeric(),
  p_value_area = numeric(),
  t_value_thickness = numeric(),
  p_value_thickness = numeric(),
  t_value_volume = numeric(),
  p_value_volume = numeric(),
  t_value_lgi = numeric(),
  p_value_lgi = numeric()
)
language <- "list_resid"


for(i in 1:length(columns_to_process)){
  col_name <- columns_to_process[i];
  print(paste("Processing column:", col_name));
  
  # model_area <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal"), data = cleaned_df_area);
  # model_thickness <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal"), data = cleaned_df_thickness);
  # model_volume <- lm(paste(col_name , " ~ ",language," + age_at_test  + sex + nonverbal"), data = cleaned_df_volume);
  # model_lgi <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal"), data = cleaned_df_lgi);
  
  # model_area <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV"), data = cleaned_df_area);
  # model_thickness <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV"), data = cleaned_df_thickness);
  # model_volume <- lm(paste(col_name , " ~ ",language," + age_at_test  + sex + nonverbal + eTIV"), data = cleaned_df_volume);
  # model_lgi <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV"), data = cleaned_df_lgi);
  
  model_area <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_area);
  model_thickness <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_thickness);
  model_volume <- lm(paste(col_name , " ~ ",language," + age_at_test  + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_volume);
  model_lgi <- lm(paste(col_name , " ~ ",language," + age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_lgi);

  model_area_summary <- summary(model_area);
  model_thickness_summary <- summary(model_thickness);
  model_volume_summary <- summary(model_volume);
  model_lgi_summary <- summary(model_lgi);
  
  t_values_area <- coef(model_area_summary)[language, "t value"]
  p_values_area <- coef(model_area_summary)[language, "Pr(>|t|)"]
  t_values_thickness <- coef(model_thickness_summary)[language, "t value"]
  p_values_thickness <- coef(model_thickness_summary)[language, "Pr(>|t|)"]
  t_values_volume <- coef(model_volume_summary)[language, "t value"]
  p_values_volume <- coef(model_volume_summary)[language, "Pr(>|t|)"]
  t_values_lgi <- coef(model_lgi_summary)[language, "t value"]
  p_values_lgi <- coef(model_lgi_summary)[language, "Pr(>|t|)"]
  
  coefficients_list <- rbind(coefficients_list,data.frame(
    column = i,
    label = col_name,
    t_value_area = t_values_area,
    p_value_area = p_values_area,
    t_value_thickness = t_values_thickness,
    p_value_thickness = p_values_thickness,
    t_value_volume = t_values_volume,
    p_value_volume = p_values_volume,
    t_value_lgi = t_values_lgi,
    p_value_lgi = p_values_lgi
  ))
}

coefficients_list$adjusted_p_value_area <- p.adjust(coefficients_list$p_value_area, method = "fdr")
coefficients_list$corrected_area <- ifelse(coefficients_list$adjusted_p_value_area*2 < 0.05, coefficients_list$t_value_area, NA)

coefficients_list$adjusted_p_value_thickness <- p.adjust(coefficients_list$p_value_thickness, method = "fdr")
coefficients_list$corrected_thickness <- ifelse(coefficients_list$adjusted_p_value_thickness*2 < 0.05, coefficients_list$t_value_thickness, NA)

coefficients_list$adjusted_p_value_volume <- p.adjust(coefficients_list$p_value_volume, method = "fdr")
coefficients_list$corrected_volume <- ifelse(coefficients_list$adjusted_p_value_volume*2 < 0.05, coefficients_list$t_value_volume, NA)

coefficients_list$adjusted_p_value_lgi <- p.adjust(coefficients_list$p_value_lgi, method = "fdr")
coefficients_list$corrected_lgi <- ifelse(coefficients_list$adjusted_p_value_lgi*2 < 0.05, coefficients_list$t_value_lgi, NA)

colors.morph = paletteer_d('rcartocolor::TealRose')
lim = 5
lim_max = 5
plot_area <- ggseg(coefficients_list, mapping = aes(fill = corrected_area),color = "white",size=0.05) +
  labs(title = "Area", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish) +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())
plot_thickness <- ggseg(coefficients_list, mapping = aes(fill = corrected_thickness),color = "black",size=0.1) +
  labs(title = "Thickness", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())
plot_volume <- ggseg(coefficients_list, mapping = aes(fill = corrected_volume),color = "black",size=0.1) +
  labs(title = "Volume", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish,na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())
plot_lgi <- ggseg(coefficients_list, mapping = aes(fill = corrected_lgi),color = "black",size=0.1) +
  labs(title = "LGI", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = "white") +
  theme_minimal(base_family = "Gill Sans") +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())
plot_lgi_with_legend <- plot_lgi + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
combined_plot <- (plot_thickness | plot_volume | plot_lgi_with_legend) +
  plot_layout(ncol = 1, nrow = 3)
print(combined_plot)

########################################################################################################################################################################################################
coefficients_area <- data.frame(
  column = integer(),
  colume_name = character(),
  t_area_list = numeric(),
  p_area_list = numeric(),
  t_area_sentence = numeric(),
  p_area_sentence = numeric()
)
coefficients_thickness <- data.frame(
  column = integer(),
  label = character(),
  t_thickness_list = numeric(),
  p_thickness_list = numeric(),
  t_thickness_sentence = numeric(),
  p_thickness_sentence = numeric()
)
coefficients_volume <- data.frame(
  column = integer(),
  colume_name = character(),
  t_volume_list = numeric(),
  p_volume_list = numeric(),
  t_volume_sentence = numeric(),
  p_volume_sentence = numeric()
)
coefficients_lgi <- data.frame(
  column = integer(),
  colume_name = character(),
  t_lgi_list = numeric(),
  p_lgi_list = numeric(),
  t_lgi_sentence = numeric(),
  p_lgi_sentence = numeric()
)
for (i in 1:length(columns_to_process)) {
  col_name <- columns_to_process[i];
  print(paste("Processing column:", col_name));
  # model_area_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal"), data = cleaned_df_area);
  # model_area_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal"), data = cleaned_df_area);
  # model_thickness_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal"), data = cleaned_df_thickness);
  # model_thickness_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal"), data = cleaned_df_thickness);
  # model_volume_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal"), data = cleaned_df_volume);
  # model_volume_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal"), data = cleaned_df_volume);
  # model_lgi_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal"), data = cleaned_df_lgi);
  # model_lgi_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal"), data = cleaned_df_lgi);
  
  model_area_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_area);
  model_area_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_area);
  model_thickness_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_thickness);
  model_thickness_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_thickness);
  model_volume_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverba + eTIV + euler_numberl"), data = cleaned_df_volume);
  model_volume_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_volume);
  model_lgi_list <- lm(paste(col_name , " ~ list * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_lgi);
  model_lgi_sentence <- lm(paste(col_name , " ~ sentence * age_at_test + sex + nonverbal + eTIV + euler_number"), data = cleaned_df_lgi);
  
  model_area_list_summary <- summary(model_area_list);
  model_thickness_list_summary <- summary(model_thickness_list);
  model_volume_list_summary <- summary(model_volume_list);
  model_lgi_list_summary <- summary(model_lgi_list);
  
  model_area_sentence_summary <- summary(model_area_sentence);
  model_thickness_sentence_summary <- summary(model_thickness_sentence);
  model_volume_sentence_summary <- summary(model_volume_sentence);
  model_lgi_sentence_summary <- summary(model_lgi_sentence);
  
  estimate_area_list <- coef(model_area_list_summary)["list:age_at_test","Estimate"]
  t_area_list <- coef(model_area_list_summary)["list:age_at_test", "t value"]      
  p_area_list <- coef(model_area_list_summary)["list:age_at_test", "Pr(>|t|)"]     
  estimate_area_sentence <- coef(model_area_sentence_summary)["sentence:age_at_test","Estimate"]
  t_area_sentence <- coef(model_area_sentence_summary)["sentence:age_at_test", "t value"] 
  p_area_sentence <- coef(model_area_sentence_summary)["sentence:age_at_test", "Pr(>|t|)"] 
  coefficients_area <- rbind(coefficients_area,data.frame(
    column = i,
    column_name = col_name,
    t_area_list = t_area_list,
    p_area_list = p_area_list,
    t_area_sentence = t_area_sentence,
    p_area_sentence = p_area_sentence
  ))
  
  estimate_thickness_list <- coef(model_thickness_list_summary)["list:age_at_test","Estimate"]
  t_thickness_list <- coef(model_thickness_list_summary)["list:age_at_test", "t value"]          
  p_thickness_list <- coef(model_thickness_list_summary)["list:age_at_test", "Pr(>|t|)"]            
  estimate_thickness_sentence <- coef(model_thickness_sentence_summary)["sentence:age_at_test","Estimate"]
  t_thickness_sentence <- coef(model_thickness_sentence_summary)["sentence:age_at_test", "t value"]  
  p_thickness_sentence <- coef(model_thickness_sentence_summary)["sentence:age_at_test", "Pr(>|t|)"]  
  coefficients_thickness <- rbind(coefficients_thickness,data.frame(
    column = i,
    label = col_name,
    t_thickness_list = t_thickness_list,
    p_thickness_list = p_thickness_list,
    t_thickness_sentence = t_thickness_sentence,
    p_thickness_sentence = p_thickness_sentence
  ))
  
  estimate_volume_list <- coef(model_volume_list_summary)["list:age_at_test","Estimate"]
  t_volume_list <- coef(model_volume_list_summary)["list:age_at_test", "t value"]              
  p_volume_list <- coef(model_volume_list_summary)["list:age_at_test", "Pr(>|t|)"]             
  estimate_volume_sentence <- coef(model_volume_sentence_summary)["sentence:age_at_test","Estimate"]
  t_volume_sentence <- coef(model_volume_sentence_summary)["sentence:age_at_test", "t value"]  
  p_volume_sentence <- coef(model_volume_sentence_summary)["sentence:age_at_test", "Pr(>|t|)"]  
  coefficients_volume <- rbind(coefficients_volume,data.frame(
    column = i,
    column_name = col_name,
    t_volume_list = t_volume_list,
    p_volume_list = p_volume_list,
    t_volume_sentence = t_volume_sentence,
    p_volume_sentence = p_volume_sentence
  ))
  
  estimate_lgi_list <- coef(model_lgi_list_summary)["list:age_at_test","Estimate"]
  t_lgi_list <- coef(model_lgi_list_summary)["list:age_at_test", "t value"]             
  p_lgi_list <- coef(model_lgi_list_summary)["list:age_at_test", "Pr(>|t|)"]             
  estimate_lgi_sentence <- coef(model_lgi_sentence_summary)["sentence:age_at_test","Estimate"]
  t_lgi_sentence <- coef(model_lgi_sentence_summary)["sentence:age_at_test", "t value"]  
  p_lgi_sentence <- coef(model_lgi_sentence_summary)["sentence:age_at_test", "Pr(>|t|)"] 
  coefficients_lgi <- rbind(coefficients_lgi,data.frame(
    column = i,
    column_name = col_name,
    t_lgi_list = t_lgi_list,
    p_lgi_list = p_lgi_list,
    t_lgi_sentence = t_lgi_sentence,
    p_lgi_sentence = p_lgi_sentence
  ))
}

coefficients_thickness$adjust_p_value_list <- p.adjust(coefficients_thickness$p_thickness_list,method = "fdr")
coefficients_thickness$corrected_t_list <- ifelse(coefficients_thickness$adjust_p_value_list*4 < 0.05, coefficients_thickness$t_thickness_list, NA)
coefficients_thickness$adjust_p_value_sentence <- p.adjust(coefficients_thickness$p_thickness_sentence,method = "fdr")
coefficients_thickness$corrected_t_sentence <- ifelse(coefficients_thickness$adjust_p_value_sentence*4 < 0.05, coefficients_thickness$t_thickness_sentence, NA)

coefficients_area$adjust_p_value_list <- p.adjust(coefficients_area$p_area_list,method = "fdr")
coefficients_area$corrected_t_list <- ifelse(coefficients_area$adjust_p_value_list*4 < 0.05, coefficients_area$t_area_list, NA)
coefficients_area$adjust_p_value_sentence <- p.adjust(coefficients_area$p_area_sentence,method = "fdr")
coefficients_area$corrected_t_sentence <- ifelse(coefficients_area$adjust_p_value_sentence*4 < 0.05, coefficients_area$t_area_sentence, NA)

coefficients_volume$adjust_p_value_list <- p.adjust(coefficients_volume$p_volume_list,method = "fdr")
coefficients_volume$corrected_t_list <- ifelse(coefficients_volume$adjust_p_value_list*4 < 0.05, coefficients_volume$t_volume_list, NA)
coefficients_volume$adjust_p_value_sentence <- p.adjust(coefficients_volume$p_volume_sentence,method = "fdr")
coefficients_volume$corrected_t_sentence <- ifelse(coefficients_volume$adjust_p_value_sentence*4 < 0.05, coefficients_volume$t_volume_sentence, NA)

coefficients_lgi$adjust_p_value_list <- p.adjust(coefficients_lgi$p_lgi_list,method = "fdr")
coefficients_lgi$corrected_t_list <- ifelse(coefficients_lgi$adjust_p_value_list*4 < 0.05, coefficients_lgi$t_lgi_list, NA)
coefficients_lgi$adjust_p_value_sentence <- p.adjust(coefficients_lgi$p_lgi_sentence,method = "fdr")
coefficients_lgi$corrected_t_sentence <- ifelse(coefficients_lgi$adjust_p_value_sentence*4 < 0.05, coefficients_lgi$t_lgi_sentence, NA)

colors.morph = paletteer_d('rcartocolor::TealRose')
lim = 15
lim_max = 0
ggseg(coefficients_thickness, mapping = aes(fill = corrected_t_list)) +
  labs(title = " ", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-5, 5), oob = squish) +
  theme_minimal(base_family = "Gill Sans") +
  theme( plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

coefficients_list$label <- t_data_area$label
plot_lgi_with_legend <- plot_lgi + 
  theme(legend.position = 'bottom') +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
combined_plot <- (plot_thickness | plot_lgi_with_legend) +
  plot_layout(ncol = 1, nrow = 2)
print(combined_plot)

indices <- c(7, 12, 21, 24, 26,   28, 30, 41, 46, 50,   60, 62)
cleaned_df_thickness$new_group <- as.factor(cleaned_df_thickness$new_group)
plot_list <- vector("list", length(indices))
for (i in seq_along(indices)) {
  col_name <- columns_to_process[indices[i]]
  print(col_name)
  plot_list[[i]] <- ggplot(cleaned_df_thickness, aes_string(x = col_name, y = "list", color = "new_group", group = "new_group")) +
    geom_point(alpha = 0.6) +  
    geom_smooth(method = "lm", se = FALSE) + 
    labs(title = paste(col_name)) + 
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(), 
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10,10,10,10)
    ) +
    scale_color_manual(values = c("#E41A1C", "#4DAF4A","#377EB8"))
  ggsave(filename = paste0(col_name, ".pdf"), plot = plot_list[[i]], width = 2.25, height = 2)
}
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 4))
ggsave(filename = "combined_plot.pdf", plot = combined_plot, width = 12, height = 8)

grouped_data <- split(cleaned_df_thickness, cleaned_df_thickness$group_3)
cor_results <- data.frame(
  col_name = character(),
  group = character(),
  cor = numeric(),
  p_value = numeric(),
  pcor_r = numeric(),
  pcor_p = numeric(),
  stringsAsFactors = FALSE
)
indices <- c(7, 12, 21, 24, 26,   28, 30, 41, 46, 50,   60, 62)
for (i in indices) {
  col_name <- columns_to_process[i]
  for(group in names(grouped_data)) {
    group_data <- grouped_data[[group]]
    cor_test <- cor.test(group_data[[col_name]], group_data$list)
    cor_value <- cor_test$estimate
    p_value <- cor_test$p.value
    pcor_test <- pcor.test(group_data[[col_name]], group_data$list,group_data[,c("age_at_test","sex","nonverbal","eTIV","euler_number")])
    pcor_r_value <- pcor_test$estimate
    pcor_p_value <- pcor_test$p.value
    cor_results <- rbind(cor_results, data.frame(
      col_name = col_name,
      group = group,
      cor = cor_value,
      p_value = p_value,
      pcor_r = pcor_r_value,
      pcor_p = pcor_p_value,
      stringsAsFactors = FALSE
    ))
  }
}
print(cor_results)

########################################################################################################################################################################################################


########################################################################################################################################################################################################
msn_data <- readMat(".../4features_Gmat_z_score_row.mat")
Gmat_no_diag <- msn_data$Gmat
for(i in 1:68) {
  Gmat_no_diag[i, i, ] <- NA
}
msn_degree <- apply(Gmat_no_diag, c(1, 3), mean, na.rm = TRUE)
msn_degree <- t(msn_degree)
msn_data_grouped <- data.frame(
  id = cleaned_df$id, 
  msn_degree                   
) %>%
  inner_join(cleaned_df, by = c("id" = "id"))
group_mean_values <- msn_data_grouped %>%
  group_by(new_group) %>%
  summarise(across(starts_with("X"), mean)) 
t_data_area <- read_freesurfer_table(paste0(path, "aparc_area_t_value.table"), measure = "area")
group_mean_long <- group_mean_values %>%
  pivot_longer(-new_group, names_to = "brain_region", values_to = "mean_value") %>%
  mutate(
    label = rep(t_data_area$label, times = 3)  
  )
colors.morph <- paletteer_d("rcartocolor::TealRose")
lim <- 0.05
lim_max <- 0.05
plot_group <- function(group_data, group_label) {
  ggseg(group_data, mapping = aes(fill = mean_value),color = "black",size=0.2) +
    labs(title = paste("MSN - Group", group_label), fill = "mean_degree") +
    scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = scales::squish,na.value = "white") +
    theme_minimal(base_family = "Gill Sans") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_blank()
    )+
    guides(
      fill = guide_colorbar(
        barwidth = 1,
        barheight = 10
      )
    )+
    theme(legend.position = 'none')
}
group_plots <- group_mean_long %>%
  group_by(new_group) %>%
  group_split() %>%
  lapply(function(group_data) {
    group_label <- unique(group_data$new_group)
    plot_group(group_data, group_label)
  })

grid.arrange(grobs = group_plots, ncol = 1)

residuals_matrix <- matrix(NA, nrow = nrow(msn_data_grouped), ncol = 68)
colnames(residuals_matrix) <- colnames(msn_data_grouped)[2:69]
for (i in 2:69) {
  brain_region <- colnames(msn_data_grouped)[i]
  model <- lm(msn_data_grouped[[brain_region]] ~ sex + eTIV + euler_number, data = msn_data_grouped)
  residuals_matrix[, i - 1] <- residuals(model)
}
residuals_df <- as.data.frame(residuals_matrix)
colnames(residuals_df) <- paste0("residual_", colnames(msn_data_grouped)[2:69])
msn_data_grouped_with_residuals <- cbind(msn_data_grouped, residuals_df)
anova_results <- data.frame(brain_region = character(0), df1=numeric(0), df2=numeric(0), F_value = numeric(0), p_value = numeric(0), eta_squared = numeric(0))
msn_data_grouped$new_group <- as.factor(msn_data_grouped$new_group)

for (i in 1:68) {
  brain_region <- colnames(residuals_df)[i]
  print(brain_region)
  model <- aov(residuals_df[[brain_region]] ~ new_group, data = msn_data_grouped)
  anova_summary <- summary(model)[[1]]
  SS_effect <- anova_summary["new_group", "Sum Sq"]
  SS_resid  <- anova_summary["Residuals", "Sum Sq"]
  SS_total  <- SS_effect + SS_resid
  eta_sq <- SS_effect / SS_total
  anova_results <- rbind(anova_results, data.frame(
    brain_region = brain_region,
    df1 = anova_summary["new_group", "Df"],
    df2 = anova_summary["Residuals", "Df"],
    F_value = anova_summary["new_group", "F value"],
    p_value = anova_summary["new_group", "Pr(>F)"],
    eta_squared = eta_sq
  ))
}

anova_results <- anova_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))%>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"), 
    F_adj = ifelse(p_adj < 0.05, F_value, NA), 
    F_no_corr = ifelse(p_value < 0.05,  F_value, NA)
  )

n1 <- 44
n2 <- 64
n3 <- 37
n_harm <- 3 / (1/n1 + 1/n2 + 1/n3)
anova_results$posthoc_power <- apply(anova_results, 1, function(row) {
  eta_sq <- as.numeric(row["eta_squared"])
  if (eta_sq <= 0) return(0)
  f_value <- sqrt(eta_sq / (1 - eta_sq))
  pwr_res <- pwr.anova.test(
    k = 3,
    n = n_harm,
    f = f_value,
    sig.level = 0.05
  )
  return(pwr_res$power)
})

anova_results <- anova_results %>%
  mutate(label = t_data_area$label)  
colors.morph <- c("#F1EAC8FF", "#E5B9ADFF", "#D98994FF", "#D0587EFF")
anova_plot_no_corr <- ggseg(anova_results, mapping = aes(fill = F_no_corr),color='black',size=0.2) +
  scale_fill_gradientn(colors = colors.morph, name = "F-value",na.value = "white") +
  labs(title = "ANOVA Results", fill = "F-value") +
  theme_minimal(base_family = "Gill Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 10)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      barwidth = 0.5,
      barheight = 6
    )
  )+
  theme(legend.position = 'none')

anova_plot <- ggseg(anova_results, mapping = aes(fill = F_adj),color='black',size=0.2) +
  scale_fill_gradientn(colors = colors.morph, name = "F-value",na.value = "white") +
  labs(title = "ANOVA Results", fill = "F-value") +
  theme_minimal(base_family = "Gill Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 10)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      barwidth = 0.5,
      barheight = 6
    )
  )+
  theme(legend.position = 'none')
all_plots <- c(group_plots, list(anova_plot))
grid.arrange(grobs = all_plots, ncol = 1)

significant_regions <- anova_results %>%
  filter(!is.na(F_no_corr)) %>%
  mutate(ROI = str_remove(brain_region, "residual_")) %>%
  arrange(F_no_corr) 

significant_regions <- anova_results %>%
  filter(!is.na(F_adj)) %>%
  mutate(ROI = str_remove(brain_region, "residual_")) %>%
  arrange(F_adj) 

roi_labels <- significant_regions %>%
  select(ROI, label)
valid_regions <- intersect(roi_labels$ROI, colnames(group_mean_values))
filtered_data <- group_mean_values %>%
  select(new_group, all_of(valid_regions)) %>%
  pivot_longer(cols = -new_group, names_to = "ROI", values_to = "value")
long_data <- filtered_data %>%
  left_join(roi_labels, by = "ROI") %>%
  mutate(ROI = label) %>% 
  select(-label)          
long_data$new_group <- factor(long_data$new_group, levels = c(1, 2, 3))

long_data <- long_data %>%
  mutate(ROI = ROI %>%
           str_replace("^lh_", "L ") %>%
           str_replace("^rh_", "R ") %>%
           str_replace("_", " ") %>% 
           str_to_title() 
  )

ggplot(long_data, aes(y = ROI, x = value, fill = new_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(y = "Region", x = "MSN degree", fill = "Group") +
  scale_fill_manual(
    values = c("#E41A1C", "#4DAF4A","#377EB8"),
    labels = c("Children", "Young Adults", "Middle-to-Older Adults"),
    breaks = c(1, 2, 3)
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle =-45,hjust = 0,vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "right", 
    legend.title = element_text(size = 8, face = "bold", color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 50, b = 5, l = 5)
  )+ coord_flip()

trend_labels <- long_data %>%
  group_by(ROI, new_group) %>%
  summarise(mean_value = mean(value), .groups = 'drop') %>%
  arrange(ROI, new_group) %>%
  group_by(ROI) %>%
  summarise(
    trend = case_when(
      all(diff(mean_value) > 0) ~ "↑",
      all(diff(mean_value) < 0) ~ "↓",
      TRUE ~ "-"
    ),
    .groups = 'drop'
  )
long_data2 <- left_join(long_data, trend_labels, by = "ROI")

long_data2 <- long_data2 %>%
  mutate(
    trend_rank = case_when(
      trend == "↑" ~ 1,
      trend == "↓" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(trend_rank, ROI) %>%
  mutate(ROI = factor(ROI, levels = unique(ROI)))

ggplot(long_data2, aes(y = ROI, x = value, fill = new_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  labs(y = "Region", x = "MSN degree", fill = "Group") +
  scale_fill_manual(
    values = c("#E41A1C", "#4DAF4A", "#377EB8"),
    labels = c("Children", "Young Adults", "Middle-to-Older Adults"),
    breaks = c(1, 2, 3)
  ) +
  geom_text(data = distinct(long_data2, ROI, trend),
            aes(y = ROI, x = max(long_data2$value) + 0.02, label = trend),
            inherit.aes = FALSE, size = 5, hjust = 0, color = "black") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(t = 5, r = 80, b = 5, l = 5)
  )

df1 <- read_excel(".../bhv_157_long_updated.xls")  
cleaned_df <- df1[df1$incl==1,]
msn_data <- readMat(".../4features_Gmat_25rois_z_score.mat")
Gmat_no_diag <- msn_data$Gmat
for(i in 1:25) {
  Gmat_no_diag[i, i, ] <- NA
}

msn_degree <- apply(Gmat_no_diag, c(1, 3), mean, na.rm = TRUE)
msn_degree <- t(msn_degree)

coefficients_list <- data.frame(
  column = integer(),
  msndegree_age_t_values = numeric(),
  msndegree_age_p_values = numeric(),
  
  msndegree_list_t_values = numeric(),
  msndegree_list_p_values = numeric(),
  
  msndegree_sentence_t_values = numeric(),
  msndegree_sentence_p_values = numeric()
)

for (i in 1:ncol(msn_degree)) {
  model <- lm(msn_degree[, i] ~ cleaned_df$age_at_test+cleaned_df$sex+cleaned_df$eTIV+cleaned_df$euler_number)
  model_list <- lm(msn_degree[, i] ~ cleaned_df$list+cleaned_df$age_at_test+cleaned_df$sex+cleaned_df$eTIV+cleaned_df$euler_number+cleaned_df$nonverbal)
  model_sentence <- lm(msn_degree[, i] ~ cleaned_df$sentence+cleaned_df$age_at_test+cleaned_df$sex+cleaned_df$eTIV+cleaned_df$euler_number+cleaned_df$nonverbal)
  
  coefficients_list <- rbind(coefficients_list,data.frame(
    column = i,
    msndegree_age_t_values = summary(model)$coefficients["cleaned_df$age_at_test", "t value"],
    msndegree_age_p_values = summary(model)$coefficients["cleaned_df$age_at_test","Pr(>|t|)"], 
    
    msndegree_list_t_values = summary(model_list)$coefficients["cleaned_df$list","t value"],
    msndegree_list_p_values = summary(model_list)$coefficients["cleaned_df$list","Pr(>|t|)"],
    
    msndegree_sentence_t_values = summary(model_sentence)$coefficients["cleaned_df$sentence","t value"],
    msndegree_sentence_p_values = summary(model_sentence)$coefficients["cleaned_df$sentence","Pr(>|t|)"]
    
  ))
}

coefficients_list$age_p_value_corrrected <- p.adjust(coefficients_list$msndegree_age_p_values,method = "fdr")
coefficients_list$age_t_value_corrected <- ifelse(coefficients_list$age_p_value_corrrected*2 < 0.05, coefficients_list$msndegree_age_t_values, NA)
coefficients_list$list_p_value_corrected <- p.adjust(coefficients_list$msndegree_list_p_values,method = "fdr")
coefficients_list$list_t_value_sig <- ifelse(coefficients_list$msndegree_list_p_values*2 < 0.05, coefficients_list$msndegree_list_t_values, NA)
coefficients_list$sentence_p_value_corrected <- p.adjust(coefficients_list$msndegree_sentence_p_values,method = "fdr")
coefficients_list$sentence_t_value_sig <- ifelse(coefficients_list$msndegree_sentence_p_values*2 < 0.05, coefficients_list$msndegree_sentence_t_values, NA)

colors.morph = paletteer_d('rcartocolor::TealRose') 
lim = 3
lim_max = 3
roi_label <- read_xlsx(".../select_rois.xlsx",sheet = "Sheet7")
coefficients_list$label <- roi_label$label

ggseg(coefficients_list, mapping = aes(fill = age_t_value_corrected),color = "black",size=0.1) +
  labs(title = "MSN", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = 'white') +
  theme_minimal(base_family = "Gill Sans") +
  theme( plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

ggseg(coefficients_list, mapping = aes(fill = list_t_value_sig),color = "black",size=0.1) +
  labs(title = "MSN", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = 'white') +
  theme_minimal(base_family = "Gill Sans") +
  theme( plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())

ggseg(coefficients_list, mapping = aes(fill = sentence_t_value_sig),color = "black",size=0.1) +
  labs(title = "MSN", fill = "t_value") +
  scale_fill_gradientn(colors = colors.morph, limits = c(-lim, lim_max), oob = squish, na.value = 'white') +
  theme_minimal(base_family = "Gill Sans") +
  theme( plot.title = element_text(hjust = 0.5, size = 14),axis.title = element_blank())











