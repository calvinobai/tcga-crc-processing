# 加载必要的包
library(here)
library(dplyr)

# 设置项目路径
project_dir <- here()
results_dir <- here("results")

# 加载处理后的数据
load(file.path(results_dir, "tcga_crc_processed_data.RData"))

# 检查临床数据中的生存状态分布
cat("临床数据中的生存状态分布:\n")
print(table(clinical_data$vital_status, useNA = "ifany"))

# 检查"Alive"病人的days_to_last_follow_up是否大量缺失
alive_patients <- clinical_data$vital_status == "Alive"
cat("\n'Alive'病人的days_to_last_follow_up统计:\n")
print(summary(clinical_data$days_to_last_follow_up[alive_patients]))
cat("'Alive'病人中days_to_last_follow_up为NA的比例:", 
    sum(is.na(clinical_data$days_to_last_follow_up[alive_patients])) / sum(alive_patients), "\n")

# 检查"Dead"病人的days_to_death是否也有缺失
dead_patients <- clinical_data$vital_status == "Dead"
cat("\n'Dead'病人的days_to_death统计:\n")
print(summary(clinical_data$days_to_death[dead_patients]))
cat("'Dead'病人中days_to_death为NA的比例:", 
    sum(is.na(clinical_data$days_to_death[dead_patients])) / sum(dead_patients), "\n")

# 创建一个中间的生存数据框，使用更明确的生存时间计算逻辑
survival_data_improved <- clinical_data %>%
  select(submitter_id, days_to_death, days_to_last_follow_up, vital_status, ajcc_pathologic_stage, age_at_diagnosis) %>%
  mutate(
    # 更明确的生存时间计算逻辑
    survival_time_improved = case_when(
      vital_status == "Dead" & !is.na(days_to_death) ~ days_to_death,
      vital_status == "Alive" & !is.na(days_to_last_follow_up) ~ days_to_last_follow_up,
      TRUE ~ NA_real_
    ),
    # 原始脚本中的生存时间计算逻辑
    survival_time_original = case_when(
      !is.na(days_to_death) ~ days_to_death,
      !is.na(days_to_last_follow_up) ~ days_to_last_follow_up,
      TRUE ~ NA_real_
    ),
    survival_status = case_when(
      vital_status == "Dead" ~ 1L,
      vital_status == "Alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# 比较两种生存时间计算方法的结果
cat("\n两种生存时间计算方法的NA值数量比较:\n")
cat("原始方法中survival_time为NA的数量:", sum(is.na(survival_data_improved$survival_time_original)), "\n")
cat("改进方法中survival_time为NA的数量:", sum(is.na(survival_data_improved$survival_time_improved)), "\n")

# 检查生存状态中NA的数量
cat("\n生存状态中NA的数量:", sum(is.na(survival_data_improved$survival_status)), "\n")

# 检查有效的生存数据（原始方法）
valid_survival_original <- !is.na(survival_data_improved$survival_time_original) & 
                           !is.na(survival_data_improved$survival_status)
cat("\n使用原始方法的有效生存数据数量:", sum(valid_survival_original), "\n")
cat("其中死亡事件数量:", sum(survival_data_improved$survival_status[valid_survival_original] == 1), "\n")

# 检查有效的生存数据（改进方法）
valid_survival_improved <- !is.na(survival_data_improved$survival_time_improved) & 
                           !is.na(survival_data_improved$survival_status)
cat("\n使用改进方法的有效生存数据数量:", sum(valid_survival_improved), "\n")
cat("其中死亡事件数量:", sum(survival_data_improved$survival_status[valid_survival_improved] == 1), "\n")

# 检查为什么有些病人的生存时间为NA
cat("\n检查'Alive'病人中缺失days_to_last_follow_up的情况:\n")
alive_missing_followup <- alive_patients & is.na(clinical_data$days_to_last_follow_up)
cat("'Alive'病人中缺失days_to_last_follow_up的数量:", sum(alive_missing_followup), "\n")
cat("占所有'Alive'病人的比例:", sum(alive_missing_followup) / sum(alive_patients), "\n")

# 提出解决方案：对于"Alive"但缺失days_to_last_follow_up的病人，使用其他可用的时间信息
# 例如，可以尝试使用days_to_last_known_disease_status或其他相关字段
# 检查clinical_data中是否有其他可用的时间字段
time_columns <- grep("days_to|time", names(clinical_data), value = TRUE)
cat("\n临床数据中可能的时间相关列:\n")
print(time_columns)

# 对于每个时间列，检查其在"Alive"病人中的缺失情况
cat("\n各时间列在'Alive'病人中的缺失情况:\n")
for(col in time_columns) {
  missing_rate <- sum(is.na(clinical_data[[col]][alive_patients])) / sum(alive_patients)
  cat(col, ":", missing_rate, "\n")
}

# 创建一个更完善的生存数据处理方案
# 这里我们尝试使用多个时间字段来填补缺失的生存时间
survival_data_complete <- clinical_data %>%
  select(submitter_id, vital_status, days_to_death, days_to_last_follow_up, 
         all_of(time_columns), ajcc_pathologic_stage, age_at_diagnosis) %>%
  mutate(
    # 尝试使用多个时间字段来确定生存时间
    survival_time_complete = case_when(
      vital_status == "Dead" & !is.na(days_to_death) ~ days_to_death,
      vital_status == "Alive" & !is.na(days_to_last_follow_up) ~ days_to_last_follow_up,
      # 如果有其他可用的时间字段，可以在这里添加更多的条件
      TRUE ~ NA_real_
    ),
    survival_status = case_when(
      vital_status == "Dead" ~ 1L,
      vital_status == "Alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# 检查最终的有效生存数据
valid_survival_complete <- !is.na(survival_data_complete$survival_time_complete) & 
                           !is.na(survival_data_complete$survival_status)
cat("\n使用完善方案的有效生存数据数量:", sum(valid_survival_complete), "\n")
cat("其中死亡事件数量:", sum(survival_data_complete$survival_status[valid_survival_complete] == 1), "\n")
cat("死亡率:", sum(survival_data_complete$survival_status[valid_survival_complete] == 1) / 
    sum(valid_survival_complete), "\n")

# 保存改进后的生存数据
improved_data_file <- here("results", "tcga_crc_improved_survival.RData")
save(survival_data_improved, survival_data_complete, file = improved_data_file)
cat("\n改进后的生存数据已保存到:", improved_data_file, "\n")

# 更新数据概要文件
info_file <- here("results", "data_summary_improved.txt")
sink(info_file)
cat("TCGA CRC Data Summary (Improved)\n")
cat("====================\n")
cat("\nNumber of patients:", nrow(clinical_data))
cat("\nNumber of genes:", nrow(exp_matrix))
cat("\nNumber of patients with valid survival data (original method):", sum(valid_survival_original))
cat("\nNumber of death events (original method):", sum(survival_data_improved$survival_status[valid_survival_original] == 1))
cat("\nDeath rate (original method):", sum(survival_data_improved$survival_status[valid_survival_original] == 1) / sum(valid_survival_original))
cat("\n\nNumber of patients with valid survival data (improved method):", sum(valid_survival_improved))
cat("\nNumber of death events (improved method):", sum(survival_data_improved$survival_status[valid_survival_improved] == 1))
cat("\nDeath rate (improved method):", sum(survival_data_improved$survival_status[valid_survival_improved] == 1) / sum(valid_survival_improved))
sink()
cat("改进的数据概要信息已保存到:", info_file, "\n")
