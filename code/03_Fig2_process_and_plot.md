# Fig2

## Fig2 A

### Fig2 A Process

```R
library(Seurat)
library(scRepertoire)
library(tidyverse)
library(patchwork)

contig.output <- c("/home/data/tmp_data/zyh-BCR/CCA/01_dandelion/10x_vdj/")
contig.list <- loadContigs(input = contig.output, 
                           format = "AIRR")

samples <- c("H01", "H02", "H03", "H04", "H05", "P01", "P02", "P03", "P04", "P05")
base_path <- "/home/data/tmp_data/zyh-BCR/CCA/01_dandelion/10x_vdj/"

contig_list <- lapply(samples, function(sample) {
  file_path <- file.path(base_path, sample, "dandelion", "filtered_contig_dandelion.tsv")
  read.table(file_path, sep = "\t", header = TRUE)
})
head(contig_list[[1]])

contig.list <- loadContigs(input = contig_list, format = "AIRR")
head(contig.list[[1]])

combined.BCR <- combineBCR(contig.list, 
                           samples = c("H01", "H02", "H03", "H04", "H05", "P01", "P02", "P03", "P04", "P05"), 
                           ID = c("HD", "HD", "HD", "HD", "HD", "ITP", "ITP", "ITP", "ITP", "ITP"),
                           threshold = 0.85)

combined.BCR <- lapply(combined.BCR, function(df) {
  df$barcode <- sub("^.*?_.*?_.*?_(.*)", "\\1", df$barcode)
  df$barcode <- paste0(df$barcode, "-1") 
  return(df)
})

head(combined.BCR[[1]])

combined.BCR <- lapply(combined.BCR, function(df) {
  # 处理 barcode 列
  df$barcode <- sub("^.*?_.*?_.*?_(.*)", "\\1", df$barcode)
  df$barcode <- paste0(df$barcode, "-1")
  return(df)
})

head(combined.BCR[[1]])

cell_types <- c("Immature B", "Transitional B", "S100A8/A9 high B", "Naïve B", "CD27+ Memory B", "Plasma")
sub_objects <- list()

for (cell_type in cell_types) {
  barcodes <- colnames(scRNA)[scRNA$cell.type == cell_type]
  
  sub_objects[[cell_type]] <- lapply(combined.BCR, function(df) {
    df[df$barcode %in% barcodes, ]
  })
}

names(sub_objects) <- paste0("combined.BCR_", tolower(gsub(" ", "_", cell_types)))
sub_objects

cell_types <- c("Immature B", "Transitional B", "S100A8/A9 high B", "Naïve B", "CD27+ Memory B", "Plasma")
sub_objects <- list()

for (cell_type in cell_types) {
  barcodes <- colnames(scRNA)[scRNA$cell.type == cell_type]
  sub_objects[[cell_type]] <- lapply(combined.BCR, function(df) {
    df[df$barcode %in% barcodes, ]
  })
}

names(sub_objects) <- paste0("combined.BCR_", tolower(gsub(" ", "_", cell_types)))
sub_objects

seurat <- combineExpression(combined.BCR, scRNA, 
            cloneCall="CTgene", group.by = "sample", proportion = FALSE, 
            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=Inf))
table(seurat@meta.data$cloneSize)
```

### Fig2 A Plot

```R
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)

meta_data <- seurat@meta.data

celltype_levels <- c(
  "Immature B", "Transitional B", "Naïve B", 
  "CD27+ Memory B", "Plasma", "Total B"
)

meta_data <- meta_data %>%
  filter(cell.type %in% celltype_levels)

plot_data <- meta_data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(cell.type, group, cloneSize) %>%
  summarise(count = n(), .groups = 'drop')

total_counts <- plot_data %>%
  group_by(cell.type, group) %>%
  summarise(total_count = sum(count), .groups = 'drop')

plot_data <- plot_data %>%
  left_join(total_counts, by = c("cell.type", "group")) %>%
  mutate(percent = (count / total_count) * 100) %>%
  select(-count, -total_count)

total_b_data <- meta_data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(group, cloneSize) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(total_count = sum(count), .groups = 'drop') %>%
  mutate(percent = (count / total_count) * 100,
         cell.type = "Total B") %>%
  select(-count, -total_count)

plot_data <- bind_rows(plot_data, total_b_data)

plot_data <- plot_data %>%
  mutate(cloneSize = as.character(cloneSize)) %>% 
  mutate(cloneSize = replace_na(cloneSize, "No BCR information"),
         cloneSize = factor(cloneSize, levels = rev(c("Single (0 < X <= 1)", 
                                                      "Small (1 < X <= 5)", 
                                                      "Medium (5 < X <= 20)", 
                                                      "Large (20 < X <= 100)", 
                                                      "Hyperexpanded (100 < X <= 500)"))),
         group = ifelse(group == "Healthy", "HD", group)) %>%
  mutate(cloneSize_group = paste(cloneSize, group),
         cloneSize_group = factor(cloneSize_group, levels = rev(c("Single (0 < X <= 1) HD", 
                                                                  "Small (1 < X <= 5) HD", 
                                                                  "Medium (5 < X <= 20) HD", 
                                                                  "Large (20 < X <= 100) HD", 
                                                                  "Hyperexpanded (100 < X <= 500) HD",
                                                                  "Single (0 < X <= 1) ITP", 
                                                                  "Small (1 < X <= 5) ITP", 
                                                                  "Medium (5 < X <= 20) ITP", 
                                                                  "Large (20 < X <= 100) ITP", 
                                                                  "Hyperexpanded (100 < X <= 500) ITP"))))

plot_data$cell.type <- factor(plot_data$cell.type, levels = celltype_levels)

barwidth <- 0.4
gap <- 0.1 
group_gap <- 1 

cell_type_midpoints <- (1:length(celltype_levels)) * (group_gap + barwidth)

plot_data <- plot_data %>%
  mutate(cell_type_factor = as.numeric(factor(cell.type, levels = celltype_levels)),
         cell_type_adj = ifelse(group == "HD",
                                cell_type_factor * (group_gap + barwidth) - barwidth / 2 - gap / 2,
                                cell_type_factor * (group_gap + barwidth) + barwidth / 2 + gap / 2))

colors_hd <- c("Single (0 < X <= 1) HD" = "#92c5de", 
               "Small (1 < X <= 5) HD" = "#4393c3", 
               "Medium (5 < X <= 20) HD" = "#2166ac", 
               "Large (20 < X <= 100) HD" = "#053061", 
               "Hyperexpanded (100 < X <= 500) HD" = "#000033")

colors_itp <- c("Single (0 < X <= 1) ITP" = "#f4a582", 
                "Small (1 < X <= 5) ITP" = "#d6604d", 
                "Medium (5 < X <= 20) ITP" = "#b2182b", 
                "Large (20 < X <= 100) ITP" = "#67001f", 
                "Hyperexpanded (100 < X <= 500) ITP" = "#33000c")

clean_labels <- function(labels) {
  sapply(labels, function(label) {
    gsub(" HD| ITP", "", label)
  })
}

p <- ggplot() +
  geom_bar(data = filter(plot_data, group == "HD"),
           aes(x = cell_type_adj, y = percent, fill = cloneSize_group),
           stat = "identity",
           position = position_stack(),
           width = barwidth) +
  scale_fill_manual(name = "Clone Size (HD)", values = colors_hd, labels = clean_labels(levels(plot_data$cloneSize_group)[levels(plot_data$cloneSize_group) %in% names(colors_hd)])) +
  new_scale_fill() +
  geom_bar(data = filter(plot_data, group == "ITP"),
           aes(x = cell_type_adj, y = percent, fill = cloneSize_group),
           stat = "identity",
           position = position_stack(),
           width = barwidth) +
  scale_fill_manual(name = "Clone Size (ITP)", values = colors_itp, labels = clean_labels(levels(plot_data$cloneSize_group)[levels(plot_data$cloneSize_group) %in% names(colors_itp)])) +
  scale_x_continuous(breaks = cell_type_midpoints,
                     labels = celltype_levels) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

output_dir <- "/home/data/tmp_data/zyh-BCR/CCA/plot/p1_umap/BCR_diversity/"
output_file <- file.path(output_dir, "BCR_diversity_stacked_barplot.pdf")

ggsave(filename = output_file, plot = p, width = 10, height = 6)
```

## Fig2 B

```R
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(gridExtra)

data_1 <- db_obs
PhylipLineage_results <- list()

unique_clone_ids <- distinct(data_1, clone_id)$clone_id

for (unique_clone_id in unique_clone_ids) {
  clone_data <- subset(data_1, clone_id == unique_clone_id & locus == "IGH")
  if(nrow(clone_data) > 0) {
    max_length <- max(nchar(clone_data$sequence))
    clone_data <- clone_data[nchar(clone_data$sequence) == max_length, ]
    max_length1 <- max(nchar(clone_data$germline_alignment))
    clone_data <- clone_data[nchar(clone_data$germline_alignment) == max_length1, ]

    clone <- makeChangeoClone(clone_data, text_fields=c("c_call","cell.type"), 
                              num_fields="consensus_count", pad_end=TRUE)
    
    tryCatch({
      phylip_exec <- "/home/data/tmp_data/zyh-BCR/Phylip/phylip-3.697/exe/dnapars"
      PhylipLineage <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
      
      if (!is.null(PhylipLineage)) {
        print(paste("PhylipLineage is not null for clone_id:", unique_clone_id))
        
        num_edges <- length(E(PhylipLineage))
        
        PhylipLineage_results[[unique_clone_id]] <- num_edges
      }
    }, error = function(e) {
      message(paste("Error processing clone_id:", unique_clone_id, "-", e$message))
      next 
    })
  }
}

PhylipLineage_df <- data.frame(clone_id = names(PhylipLineage_results), num_edges = unlist(PhylipLineage_results))

top_50_PhylipLineage <- PhylipLineage_df %>%
  arrange(desc(num_edges)) %>%
  head(207)

print(top_50_PhylipLineage)

top_50_clone_ids <- top_50_PhylipLineage$clone_id
data <- db_obs

data$group <- recode(data$group, "Healthy" = "HD")

data$sample_id <- str_replace_all(data$sample_id, name_mapping)

for (id in top_50_clone_ids) {
  data_1 <- subset(data, clone_id == id)
  data_1 <- subset(data_1, locus == "IGH")
  
  clone <- makeChangeoClone(data_1, text_fields = c("c_call", "cell.type", "group", "sample_id", "cdr3_aa"), 
                            num_fields = "consensus_count", pad_end = TRUE)
  
  phylip_exec <- "/home/data/tmp_data/zyh-BCR/Phylip/phylip-3.697/exe/dnapars"
  PhylipLineage <- buildPhylipLineage(clone, phylip_exec, rm_temp = TRUE)

  V(PhylipLineage)$fill_color <- "steelblue"
  V(PhylipLineage)$fill_color[V(PhylipLineage)$name == "Germline"] <- "black"
  V(PhylipLineage)$fill_color[grepl("Inferred", V(PhylipLineage)$name)] <- "white"
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Resting Pro-B"] <- '#7E6148'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Cycling Pro-B"] <- '#91D1C2'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Cycling Pre-B"] <- '#3C5488'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Resting Pre-B"] <- '#8491B4'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Immature B"] <- '#F39B7F'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Transitional B"] <- '#4DBBD5'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "S100A8/A9 high B"] <- '#FFDC91'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Naïve B"] <- '#00A087'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "CD27+ Memory B"] <- '#EE4C97'
  V(PhylipLineage)$fill_color[V(PhylipLineage)$cell.type == "Plasma"] <- '#B09C85'
  V(PhylipLineage)$label <- V(PhylipLineage)$cell.type
  
  g <- ggraph(PhylipLineage, layout = 'tree') +
    geom_edge_link(arrow = arrow(length = unit(1, 'mm'), type = "closed"), end_cap = circle(1, 'mm')) +
    geom_node_point(aes(fill = fill_color), color = "black", size = 10, shape = 21, show.legend = TRUE) +
    geom_node_text(aes(label = ifelse(!is.na(group), paste0("[", sample_id, "] ", cdr3_aa), cdr3_aa)), 
                   nudge_y = -0.2, nudge_x = 0.1, size = 4) +
    theme_void() +
    ggtitle(paste0("Clone ID: ", id)) +
    theme(
      legend.position = c(1, 1),
      legend.justification = c("right", "top"),
      plot.margin = unit(c(1, 3, 1, 3), "cm")
    )

  g <- g + scale_fill_manual(
    values = c(
      "black" = "black", 
      "white" = "white", 
      "#7E6148" = "#7E6148", 
      "#91D1C2" = "#91D1C2", 
      "#3C5488" = "#3C5488", 
      "#8491B4" = "#8491B4", 
      "#F39B7F" = "#F39B7F", 
      "#4DBBD5" = "#4DBBD5", 
      "#FFDC91" = "#FFDC91", 
      "#00A087" = "#00A087", 
      "#EE4C97" = "#EE4C97", 
      "#B09C85" = "#B09C85",
      "steelblue" = "steelblue"
    ),
    name = "Node Color",
    labels = c(
      "black" = "Germline", 
      "white" = "Inferred", 
      "#7E6148" = "Resting Pro-B", 
      "#91D1C2" = "Cycling Pro-B", 
      "#3C5488" = "Cycling Pre-B", 
      "#8491B4" = "Resting Pre-B", 
      "#F39B7F" = "Immature B", 
      "#4DBBD5" = "Transitional B", 
      "#FFDC91" = "S100A8/A9 high B", 
      "#00A087" = "Naïve B", 
      "#EE4C97" = "CD27+ Memory B", 
      "#B09C85" = "Plasma",
      "steelblue" = "Sample"
    ),
    breaks = c("black", "white", "#7E6148", "#91D1C2", "#3C5488", "#8491B4", "#F39B7F", "#4DBBD5", "#FFDC91", "#00A087", "#EE4C97", "#B09C85", "steelblue")
  )
  
  g <- g + coord_cartesian(clip = 'off')
  
  current_sample <- unique(data_1$sample_id)

  ggsave(paste0(current_sample, "_lineage_", id, ".pdf"), g, height = 8, width = 8)
}
```

## Fig2 C

```R
library(Seurat)
library(dplyr)
library(ggplot2)

meta_data <- sce@meta.data

filtered_meta_data <- meta_data %>%
  filter(!is.na(clone_id) & clone_id != "" & clone_id != "No_contig")

ctgene_counts <- filtered_meta_data %>%
  group_by(clone_id) %>%
  summarise(cell_count = n(), .groups = 'drop')

top_ctgenes <- ctgene_counts %>%
  arrange(desc(cell_count)) %>%
  slice_head(n = 150) %>%
  pull(clone_id)
top_ctgene_data <- filtered_meta_data %>%
  filter(clone_id %in% top_ctgenes)
top_ctgene_data$clone_id <- factor(top_ctgene_data$clone_id, levels = top_ctgenes)
top_ctgene_data <- top_ctgene_data %>%
  mutate(Rank = as.numeric(factor(clone_id, levels = top_ctgenes)))
top_ctgene_data$group <- recode(top_ctgene_data$group, "Healthy" = "HD")

clonotype_group_proportion <- top_ctgene_data %>%
  group_by(Rank, group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / sum(count[group %in% c("HD", "ITP")])) %>%
  filter(group %in% c("HD", "ITP"))

hd_proportions <- clonotype_group_proportion %>%
  filter(group == "HD") %>%
  pull(proportion)

itp_proportions <- clonotype_group_proportion %>%
  filter(group == "ITP") %>%
  pull(proportion)

ks_test_result <- ks.test(hd_proportions, itp_proportions)

p_value <- ks_test_result$p.value
p_text <- paste0("p = ", format(p_value, digits = 3, scientific = TRUE))

p2 <- ggplot(clonotype_group_proportion, aes(x = Rank, y = proportion, color = group, fill = group, group = group)) +
  geom_smooth(se = TRUE, method = "loess") +
  theme_minimal() +
  labs(x = "Top 150 clonotype", y = "Proportion", title = "") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.85, 0.85),
    legend.title = element_text(face = "bold", size = 16),
    legend.key.size = unit(1.5, "lines"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_x_continuous(breaks = c(1, 50, 100, 150)) +
  scale_color_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
  scale_fill_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
  guides(fill = "none", color = guide_legend(title = "Group")) +
  annotate("text", x = 80, y = max(clonotype_group_proportion$proportion)/2.8, label = p_text, hjust = 0, size = 5)

ggsave("p2_smoothed_line_with_p_value.pdf", p2, height = 4, width = 4)
```

## Fig2 E

```R
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

cell_type_order <- c("CD27+ Memory B", "Plasma")

object <- db_obs
object$group <- recode(object$group, "Healthy" = "HD")
object <- subset(object, cell.type %in% cell_type_order)

object <- object %>%
  left_join(
    seurat@meta.data %>%
      rownames_to_column(var = "sample") %>%
      select(sample, cloneSize),
    by = c("sample" = "sample")
  )

object <- object %>%
  mutate(
    cloneSizeCategory = case_when(
      str_detect(cloneSize, "Large \\(20 < X <= 100\\)") | 
      str_detect(cloneSize, "Medium \\(5 < X <= 20\\)") ~ "Large",
      TRUE ~ "Other"
    )
  )

filtered_data <- object %>%
  filter(group == "ITP")

cdr3aa <- aminoAcidProperties(filtered_data, seq = "cdr3", trim = TRUE, label = "cdr3")

props <- cdr3aa %>%
  gather('cdr_prop', 'value', starts_with("cdr3_aa_")) %>%
  mutate(cdr_prop_short = str_remove(cdr_prop, 'cdr3_nt_aa_'))

props_summary <- props %>%
  group_by(cloneSizeCategory, cell.type, cdr_prop_short, locus) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()), 
    .groups = 'drop'
  )

calculate_ttest <- function(data) {
  ttest_results <- data %>%
    group_by(cell.type, cdr_prop_short) %>%
    filter(n_distinct(cloneSizeCategory) == 2) %>%
    summarise(p_value = t.test(value ~ cloneSizeCategory)$p.value, .groups = 'drop') %>%
    mutate(significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  return(ttest_results)
}

plot_with_significance <- function(props, original_data, locus) {
  # Add significance markers
  ttest_results <- calculate_ttest(original_data)
  props <- props %>%
    left_join(ttest_results, by = c("cell.type", "cdr_prop_short"))

  props$cloneSizeCategory <- factor(props$cloneSizeCategory, levels = c("Other", "Large"))

  dodge_width <- 0.7

  plot <- ggplot(props, aes(x = cell.type, y = mean_value, fill = cloneSizeCategory)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 0.7, preserve = "single"), width = 0.7) +
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                  position = position_dodge(width = dodge_width), width = 0.25) +
    geom_text(data = props %>% filter(cloneSizeCategory == "Large"), aes(label = significance), 
              position = position_dodge(width = dodge_width), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Other" = "#90B28D", "Large" = "#DDAEAA"),
                      breaks = c("Other", "Large")) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    facet_wrap(~cdr_prop_short, scales = 'free_y') +
    labs(title = paste(locus, 'CDR3 Amino Acid Properties (ITP Large vs Other)'), 
         x = "",
         y = "Mean value ± SE") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(plot)
}

props_summary$cell.type <- factor(props_summary$cell.type, levels = cell_type_order)

AA_plot_igh <- plot_with_significance(props_summary %>% filter(locus == "IGH"), props %>% filter(locus == "IGH"), "IGH")
AA_plot_igl <- plot_with_significance(props_summary %>% filter(locus == "IGL"), props %>% filter(locus == "IGL"), "IGL")
AA_plot_igk <- plot_with_significance(props_summary %>% filter(locus == "IGK"), props %>% filter(locus == "IGK"), "IGK")

ggsave("property_ITP_Large_vs_Other_IGH.pdf", AA_plot_igh, width = 8, height = 8)
ggsave("property_ITP_Large_vs_Other_IGL.pdf", AA_plot_igl, width = 8, height = 8)
ggsave("property_ITP_Large_vs_Other_IGK.pdf", AA_plot_igk, width = 8, height = 8)
```

## Fig2 J

### Fig2 J Pre-process

```R
scRNA_im <- subset(scRNA, subset = cell.type == "Immature B")
scRNA_tr <- subset(scRNA, subset = cell.type == "Transitional B")
scRNA_na <- subset(scRNA, subset = cell.type == "Naïve B")
scRNA_cd27 <- subset(scRNA, subset = cell.type == "CD27+ Memory B")
scRNA_pla <- subset(scRNA, subset = cell.type == "Plasma")

process_subset <- function(seurat_obj) {
  all.genes <- rownames(scRNA)
  
  seurat_obj <- NormalizeData(seurat_obj) %>% 
    FindVariableFeatures(nfeatures = 3000) %>% 
    ScaleData(features = all.genes)
  
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, dims = 1:20)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.cca", dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, cluster.name = "cca_clusters")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.cca", dims = 1:10, reduction.name = "umap.cca")
  
  return(seurat_obj)
}

scRNA_im <- process_subset(scRNA_im)
scRNA_tr <- process_subset(scRNA_tr)
scRNA_na <- process_subset(scRNA_na)
scRNA_cd27 <- process_subset(scRNA_cd27)
scRNA_pla <- process_subset(scRNA_pla)
```

### Fig2 J Plot

```R
library(Seurat)
library(patchwork)

cluster_colors <- c(
  "0" = "#8987B7",
  "1" = "#6EA6A5",
  "2" = "#97BCA3",
  "3" = "#51859D",
  "4" = "#D35D6C",
  "5" = "#5F7174",
  "6" = "#EED0D8",
  "7" = "#C7ABB4"
)

p1 <- DimPlot(scRNA_im, reduction = "umap.cca", label = TRUE) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Immature B")
p2 <- DimPlot(scRNA_tr, reduction = "umap.cca", label = TRUE) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Transitional B")
p3 <- DimPlot(scRNA_na, reduction = "umap.cca", label = TRUE) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Naïve B")
p4 <- DimPlot(scRNA_cd27, reduction = "umap.cca", label = TRUE) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("CD27+ Memory B")
p5 <- DimPlot(scRNA_pla, reduction = "umap.cca", label = TRUE) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Plasma")

combined_plot <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2)

ggsave("subset.pdf", combined_plot, width = 10, height = 10)
```

## Fig2 L

```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

cluster_colors <- c(
  "0" = "#8987B7",
  "1" = "#6EA6A5",
  "2" = "#97BCA3",
  "3" = "#51859D",
  "4" = "#D35D6C",
  "5" = "#5F7174",
  "6" = "#EED0D8",
  "7" = "#C7ABB4"
)

get_cluster_proportion <- function(seurat_obj) {
  cluster_distribution <- seurat_obj@meta.data %>%
    group_by(cca_clusters, group) %>%
    summarise(count = n()) %>%
    ungroup()
  
  cluster_distribution <- cluster_distribution %>%
    group_by(group) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  return(cluster_distribution)
}

plot_cluster_proportion <- function(cluster_distribution, title) {
  ggplot(cluster_distribution, aes(x = group, y = proportion, fill = as.factor(cca_clusters))) +
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    scale_fill_manual(values = cluster_colors) +
    ggtitle(title) +
    xlab("Group (HD vs ITP)") +
    ylab("Proportion of Cells") +
    labs(fill = "Clusters") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

prop_im <- get_cluster_proportion(scRNA_im)
prop_tr <- get_cluster_proportion(scRNA_tr)
prop_na <- get_cluster_proportion(scRNA_na)
prop_cd27 <- get_cluster_proportion(scRNA_cd27)
prop_pla <- get_cluster_proportion(scRNA_pla)

p1 <- plot_cluster_proportion(prop_im, "Immature B")
p2 <- plot_cluster_proportion(prop_tr, "Transitional B")
p3 <- plot_cluster_proportion(prop_na, "Naïve B")
p4 <- plot_cluster_proportion(prop_cd27, "CD27+ Memory B")
p5 <- plot_cluster_proportion(prop_pla, "Plasma")

combined_barplot <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2)

ggsave("subset_cluster_proportion.pdf", combined_barplot, width = 8, height = 10)
```

## Fig2 M

```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)

plot_difference_bar_chart <- function(seurat_obj, title) {
  meta_data <- seurat_obj@meta.data
  cell_type_order <- sort(unique(meta_data$cca_clusters))
  meta_data$cell.type <- factor(meta_data$cca_clusters, levels = cell_type_order)
  
  sample_cell_counts <- meta_data %>%
    group_by(orig.ident, cca_clusters, group) %>%
    summarise(count = n(), .groups = 'drop')
  
  sample_total_counts <- meta_data %>%
    group_by(orig.ident) %>%
    summarise(total_count = n())
  
  sample_cell_proportion <- sample_cell_counts %>%
    left_join(sample_total_counts, by = "orig.ident") %>%
    mutate(proportion = count / total_count * 100) %>%
    mutate(group = if_else(group == "Healthy", "HD", group))
  
  p_values <- sample_cell_proportion %>%
    group_by(cca_clusters) %>%
    summarise(
      group1 = "HD",
      group2 = "ITP",
      p.value = wilcox.test(proportion ~ group)$p.value
    ) %>%
    mutate(signif = case_when(
      p.value <= 0.0001 ~ "****",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  p <- ggplot(sample_cell_proportion, aes(x = cca_clusters, y = proportion, fill = group)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), width = 0.7, color = "black") +
    geom_errorbar(stat = "summary", fun.data = mean_cl_normal, position = position_dodge(width = 0.9), width = 0.2) +
    geom_beeswarm(dodge.width = 0.8, size = 1, color = "black") +
    scale_fill_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(0, 80),
                       breaks = seq(0, 80, 20)) +
    labs(x = "", y = "Fraction of cells (%)", fill = "Group", title = title) +
    theme_classic() +
    theme(axis.line = element_line(size = 0.8),
          text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.ticks = element_line(colour = "black"),
          strip.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 12),
          legend.position = "right",
          strip.background = element_blank()) +
    stat_pvalue_manual(p_values, x = "cca_clusters", label = "signif", y.position = 70,
                       position = position_dodge(0.9), inherit.aes = FALSE)
  
  return(p)
}

p1 <- plot_difference_bar_chart(scRNA_im, "Immature B")
p2 <- plot_difference_bar_chart(scRNA_tr, "Transitional B")
p3 <- plot_difference_bar_chart(scRNA_na, "Naïve B")
p4 <- plot_difference_bar_chart(scRNA_cd27, "CD27+ Memory B")
p5 <- plot_difference_bar_chart(scRNA_pla, "Plasma")

ggsave("Immature_B_cell_type_proportion_difference_bar_plot.pdf", plot = p1, height = 4, width = 6)
ggsave("Transitional_B_cell_type_proportion_difference_bar_plot.pdf", plot = p2, height = 4, width = 6)
ggsave("Naive_B_cell_type_proportion_difference_bar_plot.pdf", plot = p3, height = 4, width = 6)
ggsave("CD27+_Memory_B_cell_type_proportion_difference_bar_plot.pdf", plot = p4, height = 4, width = 6)
ggsave("Plasma_cell_type_proportion_difference_bar_plot.pdf", plot = p5, height = 4, width = 6)

combined_plot <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2)

ggsave("combined_cell_type_proportion_difference_bar_plot.pdf", plot = combined_plot, height = 12, width = 16)
```

## Fig2 P

```R
library(ggplot2)
library(patchwork)

cluster_colors <- c(
  "0" = "#8987B7",
  "1" = "#6EA6A5",
  "2" = "#97BCA3",
  "3" = "#51859D",
  "4" = "#D35D6C",
  "5" = "#5F7174",
  "6" = "#EED0D8",
  "7" = "#C7ABB4"
)

max_clusters <- max(length(unique(proportion_cd27$cca_clusters)), length(unique(proportion_pla$cca_clusters)))
max_y_value <- max(c(proportion_cd27$proportion, proportion_pla$proportion))

plot_large_proportion <- function(proportion_data, title, max_clusters, max_y_value) {
  n_clusters <- length(unique(proportion_data$cca_clusters))
  width_ratio <- n_clusters / max_clusters
  
  ggplot(proportion_data, aes(x = factor(cca_clusters), y = proportion, fill = factor(cca_clusters))) +
    geom_bar(stat = "identity", width = 0.7 * width_ratio) +
    scale_fill_manual(values = cluster_colors) +
    theme_minimal() +
    labs(
      title = title,
      x = "Cluster",
      y = "Proportion of Large Cells (%)",
      fill = "Cluster"
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank() 
    ) +
    coord_cartesian(ylim = c(0, max_y_value))
}

p_cd27 <- plot_large_proportion(proportion_cd27, "CD27+ Memory B", max_clusters, max_y_value)
p_pla <- plot_large_proportion(proportion_pla, "Plasma", max_clusters, max_y_value)

combined_barplot <- p_cd27 + p_pla + plot_layout(ncol = 2)
ggsave("large_cell_proportion_barplot_with_overall_border.pdf", combined_barplot, width = 10, height = 6)
```

## Fig2 R

```R
library(scRNAtoolVis)
library(ggplot2)

Idents(scRNA_pla) <- "seurat_clusters"
other_clusters <- setdiff(levels(Idents(scRNA_pla)), "2")

all_diff_results <- data.frame()

for (cluster in other_clusters) {
  diff <- FindMarkers(scRNA_pla, ident.1 = "2", ident.2 = cluster,
                      logfc.threshold = 0.25, min.pct = 0.25)
  
  diff <- diff %>% filter(p_val < 0.05)
  
  diff$gene <- rownames(diff)
  diff$group <- ifelse(diff$avg_log2FC > 0, "up", "down")
  diff$cluster <- paste0("Cluster 2 vs Cluster ", cluster)
  
  rownames(diff) <- paste0(rownames(diff), "_", cluster)
  
  all_diff_results <- rbind(all_diff_results, diff)
}

cluster_colors <- c(
  "Cluster 2 vs Cluster 0" = "#8987B7", 
  "Cluster 2 vs Cluster 1" = "#6EA6A5",
  "Cluster 2 vs Cluster 3" = "#51859D",
  "Cluster 2 vs Cluster 4" = "#D35D6C",
  "Cluster 2 vs Cluster 5" = "#5F7174",
  "Cluster 2 vs Cluster 6" = "#EED0D8",
  "Cluster 2 vs Cluster 7" = "#C7ABB4"
)

volcano_plot <- jjVolcano(diffData = all_diff_results,
                          tile.col = cluster_colors,
                          size = 3,
                          topGeneN = 5,
                          fontface = 'bold.italic',
                          polar = TRUE) +
  ylim(-8, 10)

ggsave("volcano_plot.pdf", plot = volcano_plot, width = 14, height = 14)
```
