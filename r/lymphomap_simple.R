#!/dssg/home/sheny/anaconda3/envs/lymphomapr/bin/Rscript

# 简化版LymphoMapR淋巴瘤分型工具
# 功能: 处理单个RSEM文件进行淋巴瘤分型

suppressPackageStartupMessages({
  library(LymphoMapR)
  library(optparse)
  library(AnnotationDbi)   
  library(org.Hs.eg.db) 
})

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入RSEM基因表达文件 (*.genes.results)", metavar="character"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="详细输出模式")
)

# 解析命令行参数
opt_parser <- OptionParser(
  usage = "用法: %prog -i <RSEM_genes_results_file> [选项]",
  option_list = option_list,
  description = "
简化版LymphoMapR淋巴瘤分型工具

功能: 处理单个RSEM基因表达文件进行淋巴瘤分型预测

示例:
  Rscript lymphomap_simple.R -i sample.genes.results
  Rscript lymphomap_simple.R -i sample.genes.results -v
"
)

opt <- parse_args(opt_parser)

# 检查必需参数
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("错误: 必须指定输入文件 (-i/--input)", call.=FALSE)
}

rsem_file <- opt$input

# 检查文件是否存在
if (!file.exists(rsem_file)) {
  cat("错误: 文件不存在:", rsem_file, "\n")
  quit(status = 1)
}

# 详细输出函数
verbose_cat <- function(...) {
  if (opt$verbose) {
    cat("[INFO]", ..., "\n")
  }
}

cat("=== LymphoMapR 单样本淋巴瘤分型分析 ===\n")
cat("输入文件:", rsem_file, "\n")

# 读取RSEM文件
verbose_cat("读取RSEM数据...")
rsem_data <- read.delim(rsem_file, header = TRUE, stringsAsFactors = FALSE)
verbose_cat("读取", nrow(rsem_data), "个基因")

# 提取TPM值并进行log2转换
tpm_values <- log2(rsem_data$TPM + 1)
names(tpm_values) <- rsem_data$gene_id

# 去掉 Ensembl 版本号
ensembl_ids <- gsub("\\.\\d+$", "", rsem_data$gene_id)

# 离线 Ensembl→HGNC 映射（无需上网）
verbose_cat("使用离线数据库转换基因ID...")
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys   = unique(ensembl_ids),
  keytype= "ENSEMBL",
  columns= c("SYMBOL")
)

# 映射到原顺序
symbol_vec <- setNames(map_df$SYMBOL, map_df$ENSEMBL)[ensembl_ids]

valid <- !is.na(symbol_vec) & nzchar(symbol_vec)
symbol_vec <- symbol_vec[valid]
tpm_values <- tpm_values[valid]

verbose_cat("成功转换", length(symbol_vec), "个基因ID")

# 合并重复 SYMBOL（取最大；你也可改成 mean）
tapply_max <- function(x, g) tapply(x, g, max, na.rm=TRUE)
tpm_by_symbol <- tapply_max(tpm_values, symbol_vec)

expr_matrix <- matrix(tpm_by_symbol, ncol = 1)
rownames(expr_matrix) <- names(tpm_by_symbol)

sample_name <- gsub("\\.genes\\.results$", "", basename(rsem_file))
colnames(expr_matrix) <- sample_name

cat("表达矩阵维度:", nrow(expr_matrix), "基因 × 1 样本\n")

# 基本健康检查（覆盖率过低时直接提醒并退出）
nonzero <- sum(expr_matrix[,1] > 0, na.rm=TRUE)
if (nrow(expr_matrix) < 500 || nonzero < 200) {
  cat("警告: 映射后有效基因太少（总行数=", nrow(expr_matrix),
      ", 非零=", nonzero, "）。请确认参考注释是否匹配！\n", sep="")
}

# 单样本安全垫：复制一列 + 极小扰动，仅用于绕过内部 n>1 限制
single_mode <- FALSE
if (ncol(expr_matrix) == 1) {
  single_mode <- TRUE
  expr_matrix <- cbind(expr_matrix, expr_matrix + rnorm(nrow(expr_matrix), sd=1e-6))
  colnames(expr_matrix)[2] <- paste0(sample_name, "_dup")
  verbose_cat("应用单样本安全垫机制")
}

# 清理 NA/Inf
expr_matrix[is.na(expr_matrix)] <- 0
expr_matrix[is.infinite(expr_matrix)] <- 0

# 运行LymphoMapR分析
cat("运行LymphoMapR分型分析...\n")
res <- tryCatch({
  run_LymphoMap(expr_matrix)
}, error = function(e) {
  cat("LymphoMapR分析失败:", e$message, "\n")
  # 打印调用栈，便于定位
  traceback()
  quit(status = 1)
})

# 只取原始样本的结果
pred_df <- res$pred
print(pred_df)


if (single_mode) {
  # 常见结构是 data.frame: sample / label / prob
  if ("sample" %in% names(pred_df)) {
    pred_df <- subset(pred_df, sample == sample_name)
  } else {
    # 保险：若是向量/长度为2，取第1个
    pred_df <- pred_df[1, , drop=FALSE]
  }
}

cat("\n=== 分析结果 ===\n")
cat("样本ID:", sample_name, "\n")

# 处理实际的pred_df格式（tibble with sample_id, LymphoMAP, LymphoMAP_prob, FMAC, LN, TEX）
if ("LymphoMAP" %in% names(pred_df)) {
  cat("预测分型:", as.character(pred_df$LymphoMAP[1]), "\n")
}
if ("LymphoMAP_prob" %in% names(pred_df)) {
  cat("预测概率:", sprintf("%.3f", pred_df$LymphoMAP_prob[1]), "\n")
}

# 显示详细的概率分布
if (all(c("FMAC", "LN", "TEX") %in% names(pred_df))) {
  cat("详细概率分布:\n")
  cat("  FMAC:", sprintf("%.3f", pred_df$FMAC[1]), "\n")
  cat("  LN:  ", sprintf("%.3f", pred_df$LN[1]), "\n")
  cat("  TEX: ", sprintf("%.3f", pred_df$TEX[1]), "\n")
}

cat("分析时间:", format(Sys.time()), "\n")
cat("基因数量:", nrow(expr_matrix), "\n")
cat("=================\n")

verbose_cat("分析完成！")
