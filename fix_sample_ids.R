# 加载必要的包
library(here)

# 设置项目路径
project_dir <- here()
results_dir <- here("results")

# 加载处理后的数据
load(file.path(results_dir, "tcga_crc_processed_data.RData"))

# 检查表达矩阵和生存数据的样本匹配情况
cat("原始样本匹配情况:\n")
exp_samples <- colnames(exp_matrix)
surv_samples <- survival_data$submitter_id
cat("表达矩阵样本数:", length(exp_samples), "\n")
cat("生存数据样本数:", length(surv_samples), "\n")
cat("共有样本数:", length(intersect(exp_samples, surv_samples)), "\n")

# 查看前几个样本ID
cat("\n表达矩阵前几个样本ID:\n")
print(head(exp_samples))
cat("\n生存数据前几个样本ID:\n")
print(head(surv_samples))

# 提取表达矩阵样本ID的病人部分（前12个字符）
extract_patient_id <- function(tcga_barcodes) {
  # 提取前12个字符作为病人ID
  substr(tcga_barcodes, 1, 12)
}

# 将表达矩阵的列名转换为病人ID
exp_patient_ids <- extract_patient_id(exp_samples)

# 检查转换后的匹配情况
cat("\n转换后的样本匹配情况:\n")
cat("表达矩阵病人ID数:", length(exp_patient_ids), "\n")
cat("生存数据病人ID数:", length(surv_samples), "\n")
cat("共有病人ID数:", length(intersect(exp_patient_ids, surv_samples)), "\n")

# 创建一个映射表，将表达矩阵的样本ID与生存数据的病人ID对应
sample_map <- data.frame(
  sample_id = exp_samples,
  patient_id = exp_patient_ids,
  stringsAsFactors = FALSE
)

# 将生存数据与表达矩阵按病人ID匹配
matched_samples <- merge(sample_map, survival_data, by.x = "patient_id", by.y = "submitter_id")
cat("\n匹配后的样本数:", nrow(matched_samples), "\n")

# 创建一个新的表达矩阵，列名为病人ID
exp_matrix_fixed <- exp_matrix
colnames(exp_matrix_fixed) <- exp_patient_ids

# 检查新表达矩阵和生存数据的匹配情况
cat("\n修复后的样本匹配情况:\n")
exp_fixed_samples <- colnames(exp_matrix_fixed)
cat("修复后表达矩阵样本数:", length(exp_fixed_samples), "\n")
cat("生存数据样本数:", length(surv_samples), "\n")
cat("共有样本数:", length(intersect(exp_fixed_samples, surv_samples)), "\n")

# 保存修复后的数据
fixed_data_file <- here("results", "tcga_crc_processed_data_fixed.RData")
save(exp_matrix_fixed, clinical_data, survival_data, file = fixed_data_file)
cat("\n修复后的数据已保存到:", fixed_data_file, "\n")

# 更新数据检查摘要
check_file <- here("results", "data_check_summary_fixed.txt")
sink(check_file)
cat("TCGA CRC Data Check Summary (Fixed)\n")
cat("====================\n")
cat("\n表达矩阵维度:", dim(exp_matrix_fixed))
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
cat("\n表达矩阵样本数:", length(exp_fixed_samples))
cat("\n生存数据样本数:", length(surv_samples))
cat("\n共有样本数:", length(intersect(exp_fixed_samples, surv_samples)))
sink()
cat("更新的检查结果已保存到:", check_file, "\n")
