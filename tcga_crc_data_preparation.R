# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 加载必要的包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("TCGAbiolinks", quietly = TRUE))
    BiocManager::install("TCGAbiolinks")
if (!require("SummarizedExperiment", quietly = TRUE))
    BiocManager::install("SummarizedExperiment")
if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr")
if (!require("here", quietly = TRUE))
    install.packages("here")
if (!require("logger", quietly = TRUE))
    install.packages("logger")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)
library(logger)

# 配置参数
PROJECTS <- c("TCGA-COAD", "TCGA-READ")
PRIMARY_TUMOR_CODE <- "01"
MAX_DOWNLOAD_ATTEMPTS <- 3
DOWNLOAD_RETRY_DELAY <- 60  # 秒

# 设置项目结构和数据保存路径
project_dir <- here()
data_dir <- here("data")
results_dir <- here("results")
log_dir <- here("logs")

# 创建必要的目录
for (dir in c(data_dir, results_dir, log_dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# 设置日志
log_file <- file.path(log_dir, format(Sys.time(), "tcga_download_%Y%m%d_%H%M%S.log"))
log_appender(appender_file(log_file))

# 辅助函数：从TCGA条码中提取病人ID
extract_patient_id <- function(tcga_barcodes) {
    # 提取前12个字符作为病人ID
    substr(tcga_barcodes, 1, 12)
}

# 辅助函数：从TCGA条码中提取样本类型
extract_sample_type <- function(tcga_barcodes) {
    # 提取第14-15位字符作为样本类型代码
    substr(tcga_barcodes, 14, 15)
}

# 函数：下载和处理TCGA数据
download_tcga_data <- function(project) {
    log_info(sprintf("开始处理项目: %s", project))

    # 查询RNA-seq数据
    query <- GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    )

    # 下载数据（带重试机制）
    for (attempt in 1:MAX_DOWNLOAD_ATTEMPTS) {
        tryCatch({
            log_info(sprintf("尝试下载 %s 数据 (第 %d 次)", project, attempt))
            GDCdownload(query)
            log_info("下载成功")
            break
        }, error = function(e) {
            log_error(sprintf("下载失败: %s", e$message))
            if (attempt < MAX_DOWNLOAD_ATTEMPTS) {
                log_info(sprintf("等待 %d 秒后重试...", DOWNLOAD_RETRY_DELAY))
                Sys.sleep(DOWNLOAD_RETRY_DELAY)
            } else {
                stop("达到最大重试次数，下载失败")
            }
        })
    }

    # 准备数据
    log_info("准备数据...")
    exp_data <- GDCprepare(query)

    # 获取表达矩阵
    exp_matrix <- assay(exp_data)

    # 筛选原发肿瘤样本
    sample_types <- extract_sample_type(colnames(exp_matrix))
    primary_tumor_samples <- sample_types == PRIMARY_TUMOR_CODE
    exp_matrix <- exp_matrix[, primary_tumor_samples]

    # 查询临床数据
    log_info("获取临床数据...")
    clinical_query <- GDCquery_clinic(project = project)

    # 保存中间结果
    checkpoint_file <- file.path(results_dir, sprintf("%s_checkpoint.RData", project))
    save(exp_matrix, clinical_query, file = checkpoint_file)
    log_info(sprintf("项目 %s 的中间结果已保存", project))

    return(list(expression = exp_matrix, clinical = clinical_query))
}

# 下载和处理所有项目数据
project_data <- list()
for (project in PROJECTS) {
    project_data[[project]] <- download_tcga_data(project)
}

# 合并数据前的检查
log_info("检查和合并数据...")

# 检查基因是否一致
gene_lists <- lapply(project_data, function(x) rownames(x$expression))
common_genes <- Reduce(intersect, gene_lists)
if (length(common_genes) == 0) {
    stop("没有找到共同的基因")
}

# 合并表达数据和临床数据
exp_matrix <- do.call(cbind, lapply(project_data, function(x) {
    x$expression[common_genes, ]
}))

# 处理临床数据列不匹配的问题
log_info("处理临床数据...")

# 获取所有临床数据的列名
clinical_colnames <- lapply(project_data, function(x) colnames(x$clinical))

# 找出共有的列
common_clinical_cols <- Reduce(intersect, clinical_colnames)
log_info(sprintf("共有的临床数据列数: %d", length(common_clinical_cols)))

# 只保留共有的列
clinical_data <- do.call(rbind, lapply(project_data, function(x) {
    x$clinical[, common_clinical_cols, drop = FALSE]
}))

# 样本匹配和去重
log_info("进行样本匹配和去重...")
exp_patient_ids <- extract_patient_id(colnames(exp_matrix))
clinical_patient_ids <- clinical_data$submitter_id

common_patients <- intersect(exp_patient_ids, clinical_patient_ids)
log_info(sprintf("找到 %d 个匹配的病人", length(common_patients)))

if (length(common_patients) == 0) {
    stop("没有找到匹配的病人")
}

# 为每个病人选择第一个样本
exp_df <- data.frame(
    sample_id = colnames(exp_matrix),
    patient_id = exp_patient_ids,
    stringsAsFactors = FALSE
)
# 只保留common_patients中的病人，并且每个病人只取第一个样本
selected_samples <- exp_df %>%
    filter(patient_id %in% common_patients) %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()

log_info(sprintf("去重后的样本数量: %d", nrow(selected_samples)))

# 使用选定的样本筛选表达矩阵
exp_matrix <- exp_matrix[, selected_samples$sample_id]

# 根据匹配的病人ID筛选临床数据
clinical_data <- clinical_data[clinical_data$submitter_id %in% selected_samples$patient_id, ]

# 加载高质量的TCGA临床和生存数据
log_info("加载高质量的TCGA临床和生存数据...")
tcga_cdr_file <- here("clinical", "TCGA-CDR.csv")

if (!file.exists(tcga_cdr_file)) {
    stop("TCGA-CDR.csv文件不存在，请确保文件已放在clinical目录中")
}

# 读取TCGA-CDR数据
tcga_cdr <- read.csv(tcga_cdr_file, stringsAsFactors = FALSE, na.strings = c("NA", "#N/A", ""))

# 筛选结直肠癌数据（COAD和READ）
crc_clinical_curated <- tcga_cdr %>%
    filter(type %in% c("COAD", "READ"))

log_info(sprintf("找到 %d 个结直肠癌病人的高质量临床数据", nrow(crc_clinical_curated)))

# 将病人条码转换为与表达矩阵匹配的格式
# 在表达矩阵中，样本名可能是完整的TCGA条码，而在CDR文件中是病人条码
# 我们需要将表达矩阵的样本名转换为病人条码进行匹配
exp_patient_ids <- extract_patient_id(colnames(exp_matrix))

# 将表达矩阵的列名与病人条码关联
exp_sample_map <- data.frame(
    sample_id = colnames(exp_matrix),
    patient_id = exp_patient_ids,
    stringsAsFactors = FALSE
)

# 找出在高质量临床数据中有对应记录的病人
common_patients <- intersect(exp_patient_ids, crc_clinical_curated$bcr_patient_barcode)
log_info(sprintf("找到 %d 个匹配的病人", length(common_patients)))

if (length(common_patients) == 0) {
    stop("没有找到匹配的病人")
}

# 为每个病人选择第一个样本
selected_samples <- exp_sample_map %>%
    filter(patient_id %in% common_patients) %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()

log_info(sprintf("去重后的样本数量: %d", nrow(selected_samples)))

# 使用选定的样本筛选表达矩阵
exp_matrix <- exp_matrix[, selected_samples$sample_id]

# 根据匹配的病人 ID 筛选高质量临床数据
crc_clinical_filtered <- crc_clinical_curated %>%
    filter(bcr_patient_barcode %in% selected_samples$patient_id)

# 创建生存数据框
survival_data <- crc_clinical_filtered %>%
    select(
        submitter_id = bcr_patient_barcode,
        ajcc_pathologic_stage = ajcc_pathologic_tumor_stage,
        age_at_diagnosis = age_at_initial_pathologic_diagnosis,
        OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time
    ) %>%
    # 将生存状态转换为数值（0=存活，1=死亡）
    mutate(
        OS = as.integer(OS),
        DSS = as.integer(DSS),
        DFI = as.integer(DFI),
        PFI = as.integer(PFI)
    )

# 检查是否存在缺失的生存时间
missing_os <- sum(is.na(survival_data$OS.time))
if (missing_os > 0) {
    log_warn(sprintf("发现 %d 个病人缺失总生存时间信息", missing_os))
}

# 保存最终处理后的数据
processed_data_file <- here("results", "tcga_crc_processed_data.RData")
save(exp_matrix, clinical_data, survival_data,
     file = processed_data_file)

# 保存数据统计信息
info_file <- here("results", "data_summary.txt")
sink(info_file)
cat("TCGA CRC Data Summary\n")
cat("====================\n")
cat("\nNumber of patients:", nrow(survival_data))
cat("\nNumber of genes:", nrow(exp_matrix))

# 总生存统计
valid_os <- !is.na(survival_data$OS) & !is.na(survival_data$OS.time)
cat("\n\nOverall Survival (OS) Statistics:")
cat("\nNumber of patients with valid OS data:", sum(valid_os))
cat("\nNumber of death events (OS):", sum(survival_data$OS[valid_os]))
cat("\nDeath rate (OS):", sum(survival_data$OS[valid_os]) / sum(valid_os))
cat("\nMedian follow-up time (OS, days):", median(survival_data$OS.time[valid_os]))

# 无病生存统计
valid_dfi <- !is.na(survival_data$DFI) & !is.na(survival_data$DFI.time)
cat("\n\nDisease-Free Interval (DFI) Statistics:")
cat("\nNumber of patients with valid DFI data:", sum(valid_dfi))
cat("\nNumber of recurrence events (DFI):", sum(survival_data$DFI[valid_dfi]))
cat("\nRecurrence rate (DFI):", sum(survival_data$DFI[valid_dfi]) / sum(valid_dfi))
cat("\nMedian time to recurrence (DFI, days):", median(survival_data$DFI.time[valid_dfi]))

# 无进展生存统计
valid_pfi <- !is.na(survival_data$PFI) & !is.na(survival_data$PFI.time)
cat("\n\nProgression-Free Interval (PFI) Statistics:")
cat("\nNumber of patients with valid PFI data:", sum(valid_pfi))
cat("\nNumber of progression events (PFI):", sum(survival_data$PFI[valid_pfi]))
cat("\nProgression rate (PFI):", sum(survival_data$PFI[valid_pfi]) / sum(valid_pfi))
cat("\nMedian time to progression (PFI, days):", median(survival_data$PFI.time[valid_pfi]))
sink()

log_info("数据处理完成！")
log_info(sprintf("处理后的数据保存在: %s", processed_data_file))
log_info(sprintf("数据概要信息保存在: %s", info_file))
