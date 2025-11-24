#解决method的q5 FDR p correct_p
###############################################################################################################################################################################













































#解决method的q1统计效力
###############################################################################################################################################################################
#ANOVA分析
###################################
#仅控制sex 不控制age
# Figure 6b. ANOVA分析 MSN differences across age groups
#除了导入包外，可单独运行
#为了与被试的 MSN 顺序匹配，不能选其他Sheet!!!!!
df1 <- read_excel("/Users/qiantingcheng/newdata_2023.12.26/data_lifespan/data_lifespan/bhv_157_long_updated_group.xls")
cleaned_df <- df1[df1$incl==1,]
msn_data <- readMat("/Users/qiantingcheng/2024.6.5processing/4features_Gmat_z_score_row.mat")  #68*68*145 MF

Gmat_no_diag <- msn_data$Gmat  # 去掉对角元素
for(i in 1:68) {
  Gmat_no_diag[i, i, ] <- NA
}

msn_degree <- apply(Gmat_no_diag, c(1, 3), mean, na.rm = TRUE)# 计算平均值，忽略对角元素
msn_degree <- t(msn_degree)  # 转置

# 将 MSN 数据转为数据框并与分组信息匹配
msn_data_grouped <- data.frame(
  id = cleaned_df$id,  # 直接使用分组信息中的 ID 为什么能犯这么低级的错误！！！！ why？？？？？
  msn_degree                   # MSN 数据，每个脑区的度信息
) %>%
  inner_join(cleaned_df, by = c("id" = "id"))  # 合并分组信息


residuals_matrix <- matrix(NA, nrow = nrow(msn_data_grouped), ncol = 68)# 初始化残差存储矩阵
colnames(residuals_matrix) <- colnames(msn_data_grouped)[2:69]  # 假设脑区列是第 2 列到第 69 列

# 循环计算每个脑区的残差
for (i in 2:69) {  # 假设第 2 列到第 69 列是脑区数据
  brain_region <- colnames(msn_data_grouped)[i]
  model <- lm(msn_data_grouped[[brain_region]] ~ sex, data = msn_data_grouped) #创建线性模型，控制混杂因素
  residuals_matrix[, i - 1] <- residuals(model)# 存储残差
}

# 将残差矩阵合并到原始数据框
residuals_df <- as.data.frame(residuals_matrix)
colnames(residuals_df) <- paste0("residual_", colnames(msn_data_grouped)[2:69])  # 给残差列命名
msn_data_grouped_with_residuals <- cbind(msn_data_grouped, residuals_df)

anova_results <- data.frame(brain_region = character(0), df1=numeric(0), df2=numeric(0), F_value = numeric(0), p_value = numeric(0), eta_squared = numeric(0))
msn_data_grouped$new_group <- as.factor(msn_data_grouped$new_group)

# for (i in 1:68) {
#   brain_region <- colnames(residuals_df)[i]
#   print(brain_region)
#   # ANOVA 分析
#   model <- aov(residuals_df[[brain_region]] ~ new_group, data = msn_data_grouped)
#   anova_summary <- summary(model)[[1]]
#   
#   # 存储结果
#   anova_results <- rbind(anova_results, data.frame(
#     brain_region = brain_region,
#     df1 = anova_summary["new_group", "Df"],
#     df2 = anova_summary["Residuals", "Df"],  
#     F_value = anova_summary["new_group", "F value"],
#     p_value = anova_summary["new_group", "Pr(>F)"]
#   ))
# }

for (i in 1:68) {
  brain_region <- colnames(residuals_df)[i]
  print(brain_region)
  # Fit model
  model <- aov(residuals_df[[brain_region]] ~ new_group, data = msn_data_grouped)
  anova_summary <- summary(model)[[1]]
  
  # Extract Sum Sq
  SS_effect <- anova_summary["new_group", "Sum Sq"]
  SS_resid  <- anova_summary["Residuals", "Sum Sq"]
  SS_total  <- SS_effect + SS_resid
  
  # Compute eta squared
  eta_sq <- SS_effect / SS_total
  
  # Store results
  anova_results <- rbind(anova_results, data.frame(
    brain_region = brain_region,
    df1 = anova_summary["new_group", "Df"],
    df2 = anova_summary["Residuals", "Df"],
    F_value = anova_summary["new_group", "F value"],
    p_value = anova_summary["new_group", "Pr(>F)"],
    eta_squared = eta_sq
  ))
}


# 多重比较校正
anova_results <- anova_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))%>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),  # FDR 校正 p 值
    F_adj = ifelse(p_adj < 0.05, F_value, NA),  # 将不显著的 F 值设为 NA
    F_no_corr = ifelse(p_value < 0.05,  F_value, NA)
  )

#################################
#post hoc power
library(pwr)
# 你的实际样本数
n1 <- 44
n2 <- 64
n3 <- 37
n_harm <- 3 / (1/n1 + 1/n2 + 1/n3)
anova_results$posthoc_power <- apply(anova_results, 1, function(row) {
  eta_sq <- as.numeric(row["eta_squared"])
  if (eta_sq <= 0) return(0)
  # η² 转 f
  f_value <- sqrt(eta_sq / (1 - eta_sq))
  # post hoc power
  pwr_res <- pwr.anova.test(
    k = 3,
    n = n_harm,     # 不等组的最佳近似
    f = f_value,
    sig.level = 0.05
  )
  return(pwr_res$power)
})





path <- "/Users/qiantingcheng/testproject/aparc__t_value/age/"
t_data_area <- read_freesurfer_table(paste0(path, "aparc_area_t_value.table"), measure = "area")
# 假设 significant_regions 包含脑区标签 (label) 和 F 值
anova_results <- anova_results %>%
  mutate(label = t_data_area$label)  # 匹配脑区标签（t_data_area 中已有）
# 保存为 CSV 文件
# write.csv(anova_results, file = "/Users/qiantingcheng/Desktop/word/2025.4.24_picture_final/5.29/anova_results_FDR_filtered.csv", row.names = TRUE)
write.csv(anova_results, file = "/Users/qiantingcheng/Desktop/word/original_data/revised_anova_results_FDR.csv", row.names = TRUE)
# 绘制 F-map
colors.morph <- c("#F1EAC8FF", "#E5B9ADFF", "#D98994FF", "#D0587EFF")

# 绘制 F-map
anova_plot_no_corr <- ggseg(anova_results, mapping = aes(fill = F_no_corr),color='black',size=0.2) +
  scale_fill_gradientn(colors = colors.morph, name = "F-value",na.value = "white") +
  labs(title = "ANOVA Results", fill = "F-value") +
  theme_minimal(base_family = "Gill Sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),        # 图标题居中
    axis.title = element_blank(),                            # 移除坐标轴标题
    legend.position = "right",                              # 图例放置在底部
    legend.title = element_text(size = 12, hjust = 0.5),     # 图例标题字体样式
    legend.text = element_text(size = 10)                    # 图例文字字体样式
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",         # 标题在颜色条上方
      title.hjust = 0.5,              # 标题居中
      label.position = "right",      # 数值在颜色条下方
      barwidth = 0.5,                 # 调整颜色条的宽度
      barheight = 6                 # 调整颜色条的高度
    )
  )+
  theme(legend.position = 'none')

# anova_plot <- ggseg(anova_results, mapping = aes(fill = F_adj),color='black',size=0.2) +
#   scale_fill_gradientn(colors = colors.morph, name = "F-value",na.value = "white") +
#   labs(title = "ANOVA Results", fill = "F-value") +
#   theme_minimal(base_family = "Gill Sans") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 14),        # 图标题居中
#     axis.title = element_blank(),                            # 移除坐标轴标题
#     legend.position = "right",                              # 图例放置在底部
#     legend.title = element_text(size = 12, hjust = 0.5),     # 图例标题字体样式
#     legend.text = element_text(size = 10)                    # 图例文字字体样式
#   ) +
#   guides(
#     fill = guide_colorbar(
#       title.position = "top",         # 标题在颜色条上方
#       title.hjust = 0.5,              # 标题居中
#       label.position = "right",      # 数值在颜色条下方
#       barwidth = 0.5,                 # 调整颜色条的宽度
#       barheight = 6                 # 调整颜色条的高度
#     )
#   )+
#   theme(legend.position = 'none')


# 将 ggseg 图加入到分组图列表中
all_plots <- c(group_plots, list(anova_plot))

# 按一列显示所有图（3组 + 1 个 ggseg 图 = 共 4 个图）
grid.arrange(grobs = all_plots, ncol = 1)
