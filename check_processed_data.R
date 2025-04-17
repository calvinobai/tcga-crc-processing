# 加载必要的包
library(here)

# 设置项目路径
project_dir <- here()
results_dir <- here("results")

# 加载处理后的数据
load(file.path(results_dir, "tcga_crc_processed_data.RData"))

# 查看数据结构
cat("表达矩阵维度:", dim(exp_matrix), "\n")
cat("临床数据维度:", dim(clinical_data), "\n")
cat("生存数据维度:", dim(survival_data), "\n")

# 查看生存数据的前几行
print(head(survival_data))

# 检查生存数据中的缺失值
cat("\n各生存终点的缺失值情况:\n")
cat("OS缺失数量:", sum(is.na(survival_data$OS)), "\n")
cat("OS.time缺失数量:", sum(is.na(survival_data$OS.time)), "\n")
cat("DSS缺失数量:", sum(is.na(survival_data$DSS)), "\n")
cat("DSS.time缺失数量:", sum(is.na(survival_data$DSS.time)), "\n")
cat("DFI缺失数量:", sum(is.na(survival_data$DFI)), "\n")
cat("DFI.time缺失数量:", sum(is.na(survival_data$DFI.time)), "\n")
cat("PFI缺失数量:", sum(is.na(survival_data$PFI)), "\n")
cat("PFI.time缺失数量:", sum(is.na(survival_data$PFI.time)), "\n")

# 检查生存数据的分布
cat("\n生存状态分布:\n")
cat("OS=0 (存活):", sum(survival_data$OS == 0, na.rm = TRUE), "\n")
cat("OS=1 (死亡):", sum(survival_data$OS == 1, na.rm = TRUE), "\n")
cat("PFI=0 (无进展):", sum(survival_data$PFI == 0, na.rm = TRUE), "\n")
cat("PFI=1 (有进展):", sum(survival_data$PFI == 1, na.rm = TRUE), "\n")

# 检查生存时间的分布
cat("\n生存时间分布:\n")
cat("OS.time中位数:", median(survival_data$OS.time, na.rm = TRUE), "\n")
cat("OS.time平均值:", mean(survival_data$OS.time, na.rm = TRUE), "\n")
cat("OS.time最小值:", min(survival_data$OS.time, na.rm = TRUE), "\n")
cat("OS.time最大值:", max(survival_data$OS.time, na.rm = TRUE), "\n")

# 检查表达矩阵的分布
cat("\n表达矩阵分布:\n")
expr_mean <- mean(as.matrix(exp_matrix))
expr_median <- median(as.matrix(exp_matrix))
expr_min <- min(as.matrix(exp_matrix))
expr_max <- max(as.matrix(exp_matrix))
cat("表达值平均数:", expr_mean, "\n")
cat("表达值中位数:", expr_median, "\n")
cat("表达值最小值:", expr_min, "\n")
cat("表达值最大值:", expr_max, "\n")

# 检查表达矩阵中的缺失值
na_count <- sum(is.na(exp_matrix))
cat("表达矩阵中的NA数量:", na_count, "\n")
cat("表达矩阵中的NA比例:", na_count / (nrow(exp_matrix) * ncol(exp_matrix)), "\n")

# 检查表达矩阵和生存数据的样本匹配情况
cat("\n样本匹配情况:\n")
exp_samples <- colnames(exp_matrix)
surv_samples <- survival_data$submitter_id
cat("表达矩阵样本数:", length(exp_samples), "\n")
cat("生存数据样本数:", length(surv_samples), "\n")
cat("共有样本数:", length(intersect(exp_samples, surv_samples)), "\n")

# 保存检查结果
check_file <- here("results", "data_check_summary.txt")
sink(check_file)
cat("TCGA CRC Data Check Summary\n")
cat("====================\n")
cat("\n表达矩阵维度:", dim(exp_matrix))
cat("\n临床数据维度:", dim(clinical_data))
cat("\n生存数据维度:", dim(survival_data))
cat("\n\n各生存终点的有效数据数量:")
cat("\nOS:", sum(!is.na(survival_data$OS) & !is.na(survival_data$OS.time)))
cat("\nDSS:", sum(!is.na(survival_data$DSS) & !is.na(survival_data$DSS.time)))
cat("\nDFI:", sum(!is.na(survival_data$DFI) & !is.na(survival_data$DFI.time)))
cat("\nPFI:", sum(!is.na(survival_data$PFI) & !is.na(survival_data$PFI.time)))
cat("\n\n生存事件数量:")
cat("\nOS死亡事件:", sum(survival_data$OS == 1, na.rm = TRUE))
cat("\nDSS死亡事件:", sum(survival_data$DSS == 1, na.rm = TRUE))
cat("\nDFI复发事件:", sum(survival_data$DFI == 1, na.rm = TRUE))
cat("\nPFI进展事件:", sum(survival_data$PFI == 1, na.rm = TRUE))
cat("\n\n表达矩阵和生存数据匹配情况:")
cat("\n表达矩阵样本数:", length(exp_samples))
cat("\n生存数据样本数:", length(surv_samples))
cat("\n共有样本数:", length(intersect(exp_samples, surv_samples)))
sink()
cat("检查结果已保存到:", check_file, "\n")
