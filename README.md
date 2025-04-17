# TCGA CRC Data Processing

这个项目提供了一套用于处理TCGA结直肠癌（CRC）数据的R脚本，包括数据下载、预处理和生存分析数据准备。

## 功能特点

- 自动下载TCGA-COAD和TCGA-READ项目的RNA-seq数据
- 处理基因表达数据和临床数据
- 准备多个生存分析终点（OS、DSS、DFI、PFI）
- 数据质量检查和验证
- 样本ID修复和匹配

## 依赖包

- BiocManager
- TCGAbiolinks
- SummarizedExperiment
- dplyr
- here
- logger

## 安装说明

1. 克隆此仓库：
```bash
git clone https://github.com/calvinobai/tcga-crc-processing.git
cd tcga-crc-processing
```

2. 安装所需的R包：
```R
install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
install.packages(c("dplyr", "here", "logger"))
```

## 使用方法

1. 准备数据目录结构：
```R
source("R/tcga_crc_data_preparation.R")
```

2. 检查数据质量：
```R
source("R/check_data.R")
```

3. 修复样本ID：
```R
source("R/fix_sample_ids.R")
```

4. 检查生存数据：
```R
source("R/check_survival_data.R")
```

## 输出文件

- `results/tcga_crc_processed_data.RData`: 处理后的表达矩阵和临床数据
- `results/data_summary.txt`: 数据处理摘要
- `logs/`: 处理日志文件

## 注意事项

- 需要足够的磁盘空间存储TCGA数据
- 建议使用稳定的网络连接
- 首次运行可能需要较长时间

## 许可证

MIT License

## 引用

如果您在研究中使用了这个项目，请引用：

[你的论文或项目引用信息]
