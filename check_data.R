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

# 统计生存状态
survival_status_counts <- table(survival_data$survival_status, useNA = "always")
print("生存状态统计:")
print(survival_status_counts)

# 计算有效的生存数据
valid_survival <- !is.na(survival_data$survival_time) & !is.na(survival_data$survival_status)
cat("有效的生存数据数量:", sum(valid_survival), "\n")

# 如果有有效的生存数据，计算中位随访时间和死亡事件数
if(sum(valid_survival) > 0) {
  valid_survival_data <- survival_data[valid_survival, ]
  cat("有效数据中的死亡事件数:", sum(valid_survival_data$survival_status), "\n")
  cat("有效数据中的中位随访时间(天):", median(valid_survival_data$survival_time), "\n")
}

# 更新数据概要文件
if(sum(valid_survival) > 0) {
  info_file <- here("results", "data_summary_updated.txt")
  sink(info_file)
  cat("TCGA CRC Data Summary (Updated)\n")
  cat("====================\n")
  cat("\nNumber of patients:", nrow(survival_data))
  cat("\nNumber of genes:", nrow(exp_matrix))
  cat("\nNumber of patients with valid survival data:", sum(valid_survival))
  cat("\nNumber of death events:", sum(valid_survival_data$survival_status))
  cat("\nMedian follow-up time (days):", median(valid_survival_data$survival_time))
  sink()
  cat("更新的数据概要信息已保存到:", info_file, "\n")
}
