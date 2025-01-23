# Fig3

## Fig3 A

```R
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)

plot_volcano_for_gene_call <- function(data, gene_call_column, output_filename) {
  object <- data
  object <- object %>%
    filter(!is.na(sample) & !is.na(!!sym(gene_call_column)) & str_starts(!!sym(gene_call_column), "IG")) %>%
    mutate(group = recode(group, "Healthy" = "HD"))
  
  object$sample.label <- paste0(object$cell.type, ",", object$sample_id, ",", object$group)
  
  sample.size <- table(object[,"sample.label"] %>% as.vector()) %>% as.data.frame()
  rownames(sample.size) <- sample.size$Var1
  sample.size$cell.type <- str_split(sample.size$Var1, ",", simplify = TRUE)[,1]
  sample.size$sample_id <- str_split(sample.size$Var1, ",", simplify = TRUE)[,2]
  sample.size$group <- str_split(sample.size$Var1, ",", simplify = TRUE)[,3]
  
  keep.samples <- sample.size[sample.size$Freq >= 5,]
  
  new <- object
  new$GeneCall <- NA
  for(i in 1:dim(new)[1]){
    tmp <- strsplit(as.vector(new[i, gene_call_column]), split = '\\*')
    new[i, "GeneCall"] <- tmp[[1]][1]
  }
  
  Usage <- prop.table(table(new$GeneCall, new[,"sample.label"]))
  for(j in 1:dim(Usage)[2]){
    Usage[,j] <- Usage[,j] / sum(Usage[,j])
  }
  
  Usage <- as.data.frame.array(t(Usage[, rownames(keep.samples)]))
  Usage$group <- keep.samples$group
  Usage[,-dim(Usage)[2]] <- Usage[,-dim(Usage)[2]] * 1000
  
  condition <- factor(Usage$group)
  colData <- data.frame(row.names = rownames(Usage), condition)
  Usage <- floor(Usage[,-dim(Usage)[2]])
  Usage <- t(Usage) %>% as.data.frame()
  Usage <- Usage + 1
  
  dds <- DESeqDataSetFromMatrix(countData = Usage, colData = colData, design = ~condition)
  dds$condition <- relevel(dds$condition, ref = "HD")
  dds2 <- DESeq(dds)
  res <- results(dds2, contrast = c("condition", "ITP", "HD")) %>% as.data.frame()
  res <- na.omit(res)
  
  prop_table <- prop.table(table(new$GeneCall, new$group), margin = 2) %>% as.data.frame()
  prop_table <- reshape2::dcast(prop_table, Var1 ~ Var2, value.var = "Freq", fill = 0)
  colnames(prop_table) <- c("GeneCall", "HD_Prop", "ITP_Prop")
  res <- merge(res, prop_table, by.x = "row.names", by.y = "GeneCall", all.x = TRUE)
  res <- na.omit(res)
  res$size <- ifelse(res$log2FoldChange > 0, res$ITP_Prop, res$HD_Prop)
  res$size <- res$size * 100
  res$significance <- "Not Significant"
  res$significance[res$padj < 0.05 & res$log2FoldChange > 0] <- "Upregulated in ITP"
  res$significance[res$padj < 0.05 & res$log2FoldChange < 0] <- "Upregulated in HD"
  
  volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance, size = size)) +
    geom_point(alpha = 1) +
    scale_color_manual(
      values = c(
        "Upregulated in ITP" = "#BFA3B2", 
        "Upregulated in HD" = "#6A92A9",
        "Not Significant" = "grey"
      ),
      name = "Group"
    ) +
    scale_size_continuous(
      name = "Usage (%)"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    labs(
      title = "",
      x = expression(Log[2]("FoldChange")),
      y = expression(-Log[10]("pvalue_adj"))
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      size = guide_legend(override.aes = list(color = "black"))
    ) +
    geom_text_repel(
      data = subset(res, padj < 0.05),
      aes(label = Row.names),
      size = 4,
      force = 2,
      fontface = "bold",
      min.segment.length = 0,
      nudge_y = 1,
      segment.size = 0.5
    )
  
  ggsave(output_filename, plot = volcano_plot, width = 8, height = 6)
}

plot_volcano_for_gene_call(data = db_obs, gene_call_column = "v_call", output_filename = "Volcano_Plot_v_call.pdf")
```

## Fig3 B

```R
library(ggplot2)
library(dplyr)

kv_upstream_genes <- c("IGKV1-NL1","IGKV3D-7","IGKV1D-8","IGKV1D-43","IGKV1D-42","IGKV2D-10","IGKV3D-11","IGKV1D-12","IGKV1D-13","IGKV2D-14","IGKV3D-15","IGKV1D-16","IGKV1D-17","IGKV6D-41","IGKV2D-18","IGKV2D-19","IGKV3D-20","IGKV6D-21","IGKV1D-22","IGKV2D-23","IGKV2D-24","IGKV3D-25","IGKV2D-26","IGKV1D-27","IGKV2D-28","IGKV2D-29","IGKV2D-30","IGKV3D-31","IGKV1D-32","IGKV1D-33","IGKV3D-34","IGKV1D-35","IGKV2D-36","IGKV1D-37","IGKV2D-38","IGKV1D-39","IGKV2D-40","IGKV2-40","IGKV1-39","IGKV2-38","IGKV1-37","IGKV2-36","IGKV1-35","IGKV3-34","IGKV1-33","IGKV1-32","IGKV3-31","IGKV2-30","IGKV2-29","IGKV2-28","IGKV1-27","IGKV2-26","IGKV3-25","IGKV2-24","IGKV2-23","IGKV1-22","IGKV6-21"
)
kv_downstream_genes <- c("IGKV3-20","IGKV2-19","IGKV2-18","IGKV1-17","IGKV1-16","IGKV3-15","IGKV2-14","IGKV1-13","IGKV1-12","IGKV3-11","IGKV2-10","IGKV1-9","IGKV1-8","IGKV3-7","IGKV1-6","IGKV1-5","IGKV2-4","IGKV7-3","IGKV5-2","IGKV4-1")
kj_upstream_genes <- c('IGKJ1', 'IGKJ2')
kj_downstream_genes <- c('IGKJ3', 'IGKJ4', 'IGKJ5')

hv_upstream_genes <- c('IGHV1-NL1', 'IGHV3-NL1', 'IGHV7-NL1', 'IGHV(III)-82', 'IGHV7-81', 'IGHV4-80', 'IGHV3-79', 'IGHV(II)-78-1', 'IGHV5-78', 'IGHV7-77', 'IGHV(III)-76-1', 'IGHV3-76', 'IGHV3-75', 'IGHV(II)-74-1', 'IGHV3-74', 'IGHV3-73', 'IGHV3-72', 'IGHV3-71', 'IGHV2-70', 'IGHV1-69D', 'IGHV1-68D', 'IGHV(III)-67-4D', 'IGHV(III)-67-3D', 'IGHV1-69-2', 'IGHV3-69-1', 'IGHV2-70D', 'IGHV1-69', 'IGHV1-68', 'IGHV(III)-67-4', 'IGHV(III)-67-3', 'IGHV(III)-67-2', 'IGHV(II)-67-1', 'IGHV1-67', 'IGHV3-66', 'IGHV(II)-65-1', 'IGHV3-65', 'IGHV3-64', 'IGHV3-63', 'IGHV(II)-62-1', 'IGHV3-62', 'IGHV4-61', 'IGHV(II)-60-1', 'IGHV3-60', 'IGHV4-59', 'IGHV1-58', 'IGHV3-57', 'IGHV7-56', 'IGHV4-55', 'IGHV3-54', 'IGHV(II)-53-1', 'IGHV3-53', 'IGHV3-52', 'IGHV(II)-51-2', 'IGHV8-51-1', 'IGHV5-51', 'IGHV3-50', 'IGHV(II)-49-1', 'IGHV3-49', 'IGHV3-48', 'IGHV(III)-47-1', 'IGHV3-47', 'IGHV(II)-46-1', 'IGHV1-46', 'IGHV1-45', 'IGHV(II)-44-2', 'IGHV(IV)-44-1', 'IGHV(III)-44', 'IGHV(II)-43-1', 'IGHV3-43', 'IGHV3-42', 'IGHV3-41', 'IGHV(II)-40-1', 'IGHV7-40', 'IGHV4-39', 'IGHV1-38-4', 'IGHV(III)-38-1D', 'IGHV3-38-3', 'IGHV(III)-44D', 'IGHV(II)-43-1D', 'IGHV3-43D', 'IGHV3-42D', 'IGHV7-40D', 'IGHV4-38-2', 'IGHV(III)-38-1', 'IGHV3-38', 'IGHV3-37', 'IGHV3-36', 'IGHV3-35', 'IGHV7-34-1', 'IGHV4-34', 'IGHV3-33-2', 'IGHV(II)-33-1', 'IGHV3-33', 'IGHV3-32', 'IGHV(II)-31-1', 'IGHV4-31', 'IGHV3-30-52', 'IGHV(II)-30-51', 'IGHV3-30-5', 'IGHV3-30-42', 'IGHV(II)-30-41', 'IGHV4-30-4', 'IGHV3-30-33', 'IGHV(II)-30-32', 'IGHV3-30-3', 'IGHV3-30-22', 'IGHV(II)-30-21', 'IGHV4-30-2', 'IGHV4-30-1', 'IGHV3-30-2', 'IGHV(II)-30-1')
hv_downstream_genes <- c(
 'IGHV3-30', 'IGHV3-29', 'IGHV(II)-28-1', 'IGHV4-28', 'IGHV7-27', 'IGHV(II)-26-2', 'IGHV(III)-26-1', 'IGHV2-26', 'IGHV(III)-25-1', 'IGHV3-25', 'IGHV1-24', 'IGHV3-23D', 'IGHV(III)-22-2D', 'IGHV(II)-22-1D', 'IGHV3-23', 'IGHV(III)-22-2', 'IGHV(II)-22-1', 'IGHV3-22', 'IGHV3-21', 'IGHV(II)-20-1', 'IGHV3-20', 'IGHV3-19', 'IGHV1-18', 'IGHV1-17', 'IGHV(III)-16-1', 'IGHV3-16', 'IGHV(II)-15-1', 'IGHV3-15', 'IGHV1-14', 'IGHV(III)-13-1', 'IGHV3-13', 'IGHV(II)-12-1', 'IGHV1-12', 'IGHV(III)-11-1', 'IGHV3-11', 'IGHV2-10', 'IGHV3-9', 'IGHV1-8', 'IGHV5-10-1', 'IGHV(III)-10-2', 'IGHV3-64D', 'IGHV3-7-2', 'IGHV(II)-7-1', 'IGHV3-7', 'IGHV3-6', 'IGHV(III)-5-2', 'IGHV(III)-5-1', 'IGHV2-5', 'IGHV7-4-1', 'IGHV(II)-4-4', 'IGHV4-4', 'IGHV1-3', 'IGHV(III)-2-1', 'IGHV1-2', 'IGHV(II)-1-1')
hj_upstream_genes <- c('IGHJ1P', 'IGHJ1', 'IGHJ2', 'IGHJ2P', 'IGHJ3') 
hj_downstream_genes <- c('IGHJ4', 'IGHJ5', 'IGHJ3P', 'IGHJ6')

lv_upstream_genes <- c(
'IGLV(I)-70', 'IGLV4-69', 'IGLV(I)-68', 'IGLV10-67', 'IGLV(IV)-66-1', 'IGLV(V)-66', 'IGLV(IV)-65', 'IGLV(IV)-64', 'IGLV(I)-63', 'IGLV1-62', 'IGLV8-61', 'IGLV4-60', 'IGLV(IV)-59', 'IGLV(V)-58', 'IGLV6-57', 'IGLV(I)-56', 'IGLV11-55', 'IGLV10-54', 'IGLV(IV)-53', 'IGLV5-52', 'IGLV1-51', 'IGLV1-50', 'IGLV9-49', 'IGLV5-48', 'IGLV1-47', 'IGLV7-46', 'IGLV5-45', 'IGLV1-44', 'IGLV7-43', 'IGLV(I)-42', 'IGLV(VII)-41-1', 'IGLV1-41', 'IGLV1-40', 'IGLV5-39', 'IGLV(I)-38', 'IGLV5-37', 'IGLV1-36', 'IGLV7-35', 'IGLV2-34', 'IGLV2-33', 'IGLV3-32', 'IGLV3-31', 'IGLV3-30', 'IGLV3-29', 'IGLV2-28', 'IGLV3-27', 'IGLV3-26', 'IGLV(VI)-25-1'
)

lv_downstream_genes <- c(
'IGLV3-25', 'IGLV3-24', 'IGLV2-23', 'IGLV(VI)-22-1', 'IGLV3-22', 'IGLV3-21', 'IGLV(I)-20', 'IGLV3-19', 'IGLV2-18', 'IGLV3-17', 'IGLV3-16', 'IGLV3-15', 'IGLV2-14', 'IGLV3-13', 'IGLV3-12', 'IGLV2-11', 'IGLV3-10', 'IGLV3-9', 'IGLV2-8', 'IGLV3-7', 'IGLV3-6', 'IGLV2-5', 'IGLV3-4', 'IGLV4-3', 'IGLV3-2', 'IGLV3-1'
)

lj_upstream_genes <- c(
'IGLJ1', 'IGLJ2'
)

lj_downstream_genes <- c(
'IGLJ3', 'IGLJ4', 'IGLJ5', 'IGLJ6', 'IGLJ7'
)

palette_dict <- c("HD" = "#6A92A9", "ITP" = "#BFA3B2")

filter_and_calculate <- function(data, igtype, upstream_genes, downstream_genes) {
  if (igtype %in% c("IGKV", "IGKJ")) {
    filter_col <- "locus_VJ"
    data <- data[!is.na(data[[filter_col]]) & grepl('IGK', data[[filter_col]]), ]
    gene_col <- if (igtype == "IGKV") "v_call_VJ_main" else "j_call_VJ_main"
  } else if (igtype %in% c("IGLV", "IGLJ")) {
    filter_col <- "locus_VJ"
    data <- data[!is.na(data[[filter_col]]) & grepl('IGL', data[[filter_col]]), ]
    gene_col <- if (igtype == "IGLV") "v_call_VJ_main" else "j_call_VJ_main"
  } else if (igtype %in% c("IGHV", "IGHJ")) {
    filter_col <- "locus_VDJ"
    data <- data[!is.na(data[[filter_col]]) & grepl('IGH', data[[filter_col]]), ]
    gene_col <- if (igtype == "IGHV") "v_call_VDJ_main" else "j_call_VDJ_main"
  } else {
    stop("Invalid igtype provided")
  }

  data <- data[!(data[[gene_col]] %in% c('None', 'No_contig')), ]

  extract_first_gene <- function(gene_list) {
    genes <- unlist(strsplit(gene_list, ','))
    for (gene in genes) {
      if (gene %in% c(upstream_genes, downstream_genes)) {
        return(gene)
      }
    }
    return(NA)
  }

  data[[gene_col]] <- sapply(data[[gene_col]], extract_first_gene)

  return(data)
}

calculate_percentage <- function(data, genes, group_col, igtype) {
  if (igtype %in% c("IGKV", "IGKJ")) {
    gene_col <- if (igtype == "IGKV") "v_call_VJ_main" else "j_call_VJ_main"
  } else if (igtype %in% c("IGLV", "IGLJ")) {
    gene_col <- if (igtype == "IGLV") "v_call_VJ_main" else "j_call_VJ_main"
  } else if (igtype %in% c("IGHV", "IGHJ")) {
    gene_col <- if (igtype == "IGHV") "v_call_VDJ_main" else "j_call_VDJ_main"
  } else {
    stop("Invalid igtype provided")
  }
  
  percentages <- c()
  for (ident in unique(data[[group_col]])) {
    sample_data <- data[data[[group_col]] == ident, ]
    gene_counts <- table(sample_data[[gene_col]])
    relevant_gene_counts <- sum(gene_counts[names(gene_counts) %in% genes])
    percentages <- c(percentages, (relevant_gene_counts / sum(gene_counts)) * 100)
  }
  return(percentages)
}

total_b_data_list <- list(
  IGKV = filter_and_calculate(sce@meta.data, "IGKV", kv_upstream_genes, kv_downstream_genes),
  IGKJ = filter_and_calculate(sce@meta.data, "IGKJ", kj_upstream_genes, kj_downstream_genes),
  IGLV = filter_and_calculate(sce@meta.data, "IGLV", lv_upstream_genes, lv_downstream_genes),
  IGLJ = filter_and_calculate(sce@meta.data, "IGLJ", lj_upstream_genes, lj_downstream_genes),
  IGHV = filter_and_calculate(sce@meta.data, "IGHV", hv_upstream_genes, hv_downstream_genes),
  IGHJ = filter_and_calculate(sce@meta.data, "IGHJ", hj_upstream_genes, hj_downstream_genes)
)

save_dir <- "/home/data/tmp_data/zyh-BCR/1222plot/up_down/"
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

for (igtype in names(total_b_data_list)) {
  total_b_data <- total_b_data_list[[igtype]]
  
  if (igtype == "IGLV") {
    upstream_genes <- lv_upstream_genes
    downstream_genes <- lv_downstream_genes
  } else if (igtype == "IGLJ") {
    upstream_genes <- lj_upstream_genes
    downstream_genes <- lj_downstream_genes
  } else if (igtype == "IGKV") {
    upstream_genes <- kv_upstream_genes
    downstream_genes <- kv_downstream_genes
  } else if (igtype == "IGKJ") {
    upstream_genes <- kj_upstream_genes
    downstream_genes <- kj_downstream_genes
  } else if (igtype == "IGHV") {
    upstream_genes <- hv_upstream_genes
    downstream_genes <- hv_downstream_genes
  } else if (igtype == "IGHJ") {
    upstream_genes <- hj_upstream_genes
    downstream_genes <- hj_downstream_genes
  } else {
    stop("Invalid igtype provided")
  }
  
  total_b_data$Group <- ifelse(grepl('HD_', total_b_data$orig.ident), 'HD', 'ITP')
  
  upstream_percentages <- data.frame(Cell_Type = character(), Group = character(), Percentage = numeric())
  downstream_percentages <- data.frame(Cell_Type = character(), Group = character(), Percentage = numeric())
  
  for (cell_type in unique(total_b_data$cell.type)) {
    cell_data <- total_b_data[total_b_data$cell.type == cell_type, ]
    for (group in c('HD', 'ITP')) {
      group_data <- cell_data[cell_data$Group == group, ]
      upstream_percentages <- rbind(
        upstream_percentages,
        data.frame(Cell_Type = cell_type, Group = group,
                   Percentage = calculate_percentage(group_data, upstream_genes, 'orig.ident', igtype))
      )
      downstream_percentages <- rbind(
        downstream_percentages,
        data.frame(Cell_Type = cell_type, Group = group,
                   Percentage = calculate_percentage(group_data, downstream_genes, 'orig.ident', igtype))
      )
    }
  }
  
  for (group in c('HD', 'ITP')) {
    group_data <- total_b_data[total_b_data$Group == group, ]
    upstream_percentages <- rbind(
      upstream_percentages,
      data.frame(Cell_Type = 'Total B', Group = group,
                 Percentage = calculate_percentage(group_data, upstream_genes, 'orig.ident', igtype))
    )
    downstream_percentages <- rbind(
      downstream_percentages,
      data.frame(Cell_Type = 'Total B', Group = group,
                 Percentage = calculate_percentage(group_data, downstream_genes, 'orig.ident', igtype))
    )
  }
  
  write.csv(upstream_percentages, file.path(save_dir, paste0("df_upstream_", igtype, ".csv")), row.names = FALSE)
  write.csv(downstream_percentages, file.path(save_dir, paste0("df_downstream_", igtype, ".csv")), row.names = FALSE)
}


gene_types <- c("IGKV", "IGKJ", "IGLV", "IGLJ", "IGHV", "IGHJ")

save_dir <- "/home/data/tmp_data/zyh-BCR/1222plot/up_down/"

cell_type_order <- c("Immature B", "Transitional B", "Naïve B", "CD27+ Memory B", "Plasma")

add_significance_annotations <- function(plot, data, x, y) {
  cell_types <- unique(data[[x]])
  for (cell_type in cell_types) {
    healthy_data <- data[data[[x]] == cell_type & data$Group == 'HD', ][[y]]
    itp_data <- data[data[[x]] == cell_type & data$Group == 'ITP', ][[y]]
    if (length(healthy_data) == 0 || length(itp_data) == 0) {
      next  
    }
    test_result <- t.test(healthy_data, itp_data)
    p_value <- test_result$p.value
    p_label <- ifelse(p_value <= 1e-4, '****', 
                      ifelse(p_value <= 1e-3, '***', 
                             ifelse(p_value <= 1e-2, '**', 
                                    ifelse(p_value <= 0.05, '*', 'ns'))))
    max_y <- max(data[data[[x]] == cell_type, ][[y]], na.rm = TRUE)
    text_size <- ifelse(p_label == 'ns', 3, 5)
    plot <- plot + annotate("text", x = cell_type, y = max_y, label = p_label, vjust = -0.5, size = text_size)
  }
  return(plot)
}

plot_gene_usage_with_error_bars <- function(data, title, filename) {
  summary_data <- data %>%
    group_by(Cell_Type, Group) %>%
    summarise(
      Mean = mean(Percentage, na.rm = TRUE),
      SE = sd(Percentage, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  p <- ggplot(summary_data, aes(x = Cell_Type, y = Mean, color = Group, group = Group)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    scale_color_manual(values = palette_dict) +
    ggtitle(title) +
    labs(x = '', y = 'Mean Percentage ± SE') +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  p <- add_significance_annotations(p, data, 'Cell_Type', 'Percentage')
  ggsave(filename, plot = p)
}

for (gene_type in gene_types) {
  df_upstream <- read.csv(file.path(save_dir, paste0("df_upstream_", gene_type, ".csv")))
  df_downstream <- read.csv(file.path(save_dir, paste0("df_downstream_", gene_type, ".csv")))
  
  df_upstream <- df_upstream %>%
    filter(Cell_Type %in% cell_type_order) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = cell_type_order))
  df_downstream <- df_downstream %>%
    filter(Cell_Type %in% cell_type_order) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = cell_type_order))
  
  plot_gene_usage_with_error_bars(df_upstream, paste('Upstream', gene_type, 'Gene Usage'), paste0('Upstream_', gene_type, '_Gene_Usage_Lineplot.pdf'))
  plot_gene_usage_with_error_bars(df_downstream, paste('Downstream', gene_type, 'Gene Usage'), paste0('Downstream_', gene_type, '_Gene_Usage_Lineplot.pdf'))
  
  for (cell_type in unique(df_upstream$Cell_Type)) {
    df_upstream_subset <- df_upstream[df_upstream$Cell_Type == cell_type, ]
    df_downstream_subset <- df_downstream[df_downstream$Cell_Type == cell_type, ]
    
    plot_gene_usage_with_error_bars(df_upstream_subset, paste('Upstream', gene_type, 'Gene Usage for', cell_type), paste0('Upstream_', gene_type, '_Gene_Usage_', cell_type, '_Lineplot.pdf'))
    plot_gene_usage_with_error_bars(df_downstream_subset, paste('Downstream', gene_type, 'Gene Usage for', cell_type), paste0('Downstream_', gene_type, '_Gene_Usage_', cell_type, '_Lineplot.pdf'))
  }
}

print("All plot saved.")
```

## Fig3 C

### Fig3 C Pre-process

```R
library(dplyr)
library(readr)

read_gene_positions <- function(file_path) {
  df <- read_csv(file_path)
  gene_positions <- setNames(df[[2]], df[[1]])
  return(gene_positions)
}

JH_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/JH_final.csv')
JK_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/JK_final.csv')
JL_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/JL_final.csv')
VH_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/VH_final.csv')
VK_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/VK_final.csv')
VL_positions <- read_gene_positions('/home/data/tmp_data/zyh-BCR/04_formal_analysis/01_dandelion/distance/location/VL_final.csv')

calculate_distance <- function(v_call, j_call, v_positions, j_positions) {
  v_genes <- if (!is.na(v_call)) unlist(strsplit(v_call, ",")) else character()
  j_genes <- if (!is.na(j_call)) unlist(strsplit(j_call, ",")) else character()
  
  v_value <- if (length(v_genes) > 0) v_positions[v_genes[1]] else NA
  j_value <- if (length(j_genes) > 0) j_positions[j_genes[1]] else NA
  
  if (!is.na(v_value) && !is.na(j_value)) {
    return(v_value - j_value)
  }
  return(NA)
}

resultsK <- list()
resultsL <- list()
resultsH <- list()

sce@meta.data$CellID <- rownames(sce@meta.data)
for (i in 1:nrow(sce@meta.data)) {
  row <- sce@meta.data[i, ]
  v_call_B_VJ <- row[['v_call_B_VJ_main']]
  j_call_B_VJ <- row[['j_call_B_VJ_main']]
  v_call_B_VDJ <- row[['v_call_B_VDJ_main']]
  j_call_B_VDJ <- row[['j_call_B_VDJ_main']]
  
  distanceK <- calculate_distance(v_call_B_VJ, j_call_B_VJ, VK_positions, JK_positions)
  distanceL <- calculate_distance(v_call_B_VJ, j_call_B_VJ, VL_positions, JL_positions)
  distanceH <- calculate_distance(v_call_B_VDJ, j_call_B_VDJ, VH_positions, JH_positions)
  
  if (!is.na(distanceK)) {
    resultsK <- append(resultsK, list(c(row['CellID'], row['cell.type'], row['group'], distanceK)))
  }
  
  if (!is.na(distanceL)) {
    resultsL <- append(resultsL, list(c(row['CellID'], row['cell.type'], row['group'], distanceL)))
  }
  
  if (!is.na(distanceH)) {
    resultsH <- append(resultsH, list(c(row['CellID'], row['cell.type'], row['group'], distanceH)))
  }
}

df_resultsK <- as.data.frame(do.call(rbind, resultsK))
names(df_resultsK) <- c('cell_id', 'cell.type', 'group', 'distanceK')

df_resultsL <- as.data.frame(do.call(rbind, resultsL))
names(df_resultsL) <- c('cell_id', 'cell.type', 'group', 'distanceL')

df_resultsH <- as.data.frame(do.call(rbind, resultsH))
names(df_resultsH) <- c('cell_id', 'cell.type', 'group', 'distanceH')

df_resultsK$cell_id <- unlist(df_resultsK$cell_id)
df_resultsK$cell.type <- unlist(df_resultsK$cell.type)
df_resultsK$group <- unlist(df_resultsK$group)
df_resultsK$distanceK <- as.numeric(unlist(df_resultsK$distanceK))

df_resultsL$cell_id <- unlist(df_resultsL$cell_id)
df_resultsL$cell.type <- unlist(df_resultsL$cell.type)
df_resultsL$group <- unlist(df_resultsL$group)
df_resultsL$distanceL <- as.numeric(unlist(df_resultsL$distanceL))

df_resultsH$cell_id <- unlist(df_resultsH$cell_id)
df_resultsH$cell.type <- unlist(df_resultsH$cell.type)
df_resultsH$group <- unlist(df_resultsH$group)
df_resultsH$distanceH <- as.numeric(unlist(df_resultsH$distanceH))

write.csv(df_resultsK, '/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceK_result.csv', row.names = FALSE)
write.csv(df_resultsL, '/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceL_result.csv', row.names = FALSE)
write.csv(df_resultsH, '/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceH_result.csv', row.names = FALSE)
```

### Fig3 C Plot

```R
library(ggplot2)
library(dplyr)
library(readr)
library(ggsignif)

sample_type.colors <- c("HD" = "#6A92A9", "ITP" = "#BFA3B2")

cell_type_order <- c("Immature B", "Transitional B", "Naïve B", "CD27+ Memory B", "Plasma")

distanceK_data <- read.csv("/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceK_result.csv")
distanceL_data <- read.csv("/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceL_result.csv")
distanceH_data <- read.csv("/home/data/tmp_data/zyh-BCR/CCA/distance/distance_result/distanceH_result.csv")

distanceL_data$distanceL <- as.numeric(distanceL_data$distanceL)
distanceK_data$distanceK <- as.numeric(distanceK_data$distanceK)
distanceH_data$distanceH <- as.numeric(distanceH_data$distanceH)

distanceK_data <- distanceK_data %>%
  mutate(group = recode(group, "Healthy" = "HD"))

distanceL_data <- distanceL_data %>%
  mutate(group = recode(group, "Healthy" = "HD"))

distanceH_data <- distanceH_data %>%
  mutate(group = recode(group, "Healthy" = "HD"))

filtered_distanceK_data <- distanceK_data %>%
  filter(cell.type %in% cell_type_order) %>%
  mutate(cell.type = factor(cell.type, levels = cell_type_order))

filtered_distanceL_data <- distanceL_data %>%
  filter(cell.type %in% cell_type_order) %>%
  mutate(cell.type = factor(cell.type, levels = cell_type_order))

filtered_distanceH_data <- distanceH_data %>%
  filter(cell.type %in% cell_type_order) %>%
  mutate(cell.type = factor(cell.type, levels = cell_type_order))

summary_distanceK <- filtered_distanceK_data %>%
  group_by(cell.type, group) %>%
  summarise(
    mean_distance = mean(distanceK, na.rm = TRUE),
    se_distance = sd(distanceK, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

summary_distanceL <- filtered_distanceL_data %>%
  group_by(cell.type, group) %>%
  summarise(
    mean_distance = mean(distanceL, na.rm = TRUE),
    se_distance = sd(distanceL, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

summary_distanceH <- filtered_distanceH_data %>%
  group_by(cell.type, group) %>%
  summarise(
    mean_distance = mean(distanceH, na.rm = TRUE),
    se_distance = sd(distanceH, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p_valuesK <- filtered_distanceK_data %>%
  group_by(cell.type) %>%
  summarise(
    group1 = "HD",
    group2 = "ITP",
    p.value = t.test(distanceK ~ group)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(signif = case_when(
    p.value <= 0.0001 ~ "****",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

p_valuesL <- filtered_distanceL_data %>%
  group_by(cell.type) %>%
  summarise(
    group1 = "HD",
    group2 = "ITP",
    p.value = t.test(distanceL ~ group)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(signif = case_when(
    p.value <= 0.0001 ~ "****",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))
  
p_valuesH <- filtered_distanceH_data %>%
  group_by(cell.type) %>%
  summarise(
    group1 = "HD",
    group2 = "ITP",
    p.value = t.test(distanceH ~ group)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(signif = case_when(
    p.value <= 0.0001 ~ "****",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

lineplot_distanceK <- ggplot(summary_distanceK, aes(x = cell.type, y = mean_distance, group = group, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance),
                width = 0.2) +
  scale_color_manual(values = sample_type.colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(-5, 0, -10, 5)
  ) +
  labs(
    x = "",
    y = "Mean Distance K",
    color = "Group"
  ) +
  ggtitle("") +
  stat_pvalue_manual(p_valuesK, x = "cell.type", label = "signif", y.position = max(summary_distanceK$mean_distance) * 1.1,
                     inherit.aes = FALSE)
lineplot_distanceL <- ggplot(summary_distanceL, aes(x = cell.type, y = mean_distance, group = group, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance),
                width = 0.2) +
  scale_color_manual(values = sample_type.colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(-5, 0, -10, 5)
  ) +
  labs(
    x = "",
    y = "Mean Distance L",
    color = "Group"
  ) +
  ggtitle("") +
  stat_pvalue_manual(p_valuesL, x = "cell.type", label = "signif", y.position = max(summary_distanceL$mean_distance) * 1.1,
                     inherit.aes = FALSE)

lineplot_distanceH <- ggplot(summary_distanceH, aes(x = cell.type, y = mean_distance, group = group, color = group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_distance - se_distance, ymax = mean_distance + se_distance),
                width = 0.2) +
  scale_color_manual(values = sample_type.colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(-5, 0, -10, 5)
  ) +
  labs(
    x = "",
    y = "Mean Distance H",
    color = "Group"
  ) +
  ggtitle("") +
  stat_pvalue_manual(p_valuesH, x = "cell.type", label = "signif", y.position = max(summary_distanceH$mean_distance) * 1.1,
                     inherit.aes = FALSE)

ggsave("distanceK_lineplot_with_signif.pdf", plot = lineplot_distanceK, width = 8, height = 5)
ggsave("distanceL_lineplot_with_signif.pdf", plot = lineplot_distanceL, width = 8, height = 5)
ggsave("distanceH_lineplot_with_signif.pdf", plot = lineplot_distanceH, width = 8, height = 5)
```

