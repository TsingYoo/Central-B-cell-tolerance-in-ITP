

# Fig1

## Fig1 B

```R
#install.packages('tidydr')
#install.packages('ggsci')
library(tidydr)
library(ggplot2)
library(ggrepel)
library(ggsci)
packageVersion("ggplot2")

umap = scRNA@reductions$umap.cca@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cellType = scRNA@meta.data$cell.type)

celltypepos <- umap %>%
  group_by(cellType) %>%
  summarise(
    umap_1 = median(umapcca_1),
    umap_2 = median(umapcca_2))

celltypepos$cellType <- factor(celltypepos$cellType, levels = c(
  "Resting Pro-B", "Cycling Pro-B", "Cycling Pre-B", "Resting Pre-B", "Immature B", 
  "Transitional B", "S100A8/A9 high B", "Naïve B", "CD27+ Memory B", 
  "Plasma"
))

p <-  DimPlot(scRNA, reduction = "umap.cca", label = FALSE, pt.size = 1.2) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = cellType,color = cellType), 
                   fontface = "bold",
                   data = celltypepos,
                   box.padding = 0.5) +
  labs(x = "umap_1", y = "umap_2")

ggsci_db_custom <- vector("list")
ggsci_db_custom$"npg"$"nrc" <- c(
  "FrenchRose" = "#EE4C97", "Shakespeare" = "#4DBBD5",
  "PersianGreen" = "#00A087", "Chambray" = "#3C5488",
  "Apricot" = "#F39B7F", "WildBlueYonder" = "#8491B4",
  "MonteCarlo" = "#91D1C2", "Sandrift" = "#B09C85",
  "RomanCoffee" = "#7E6148", "Salomie" = "#FFDC91"
 )

pal_npg_custom <- function(palette = c("nrc"), alpha = 1) {
  palette <- match.arg(palette)

  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")

  raw_cols <- ggsci_db_custom$"npg"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  scales::manual_pal(unname(alpha_cols))
}

is_installed <- function(pkg, version = NULL) {
  installed <- isNamespaceLoaded(pkg) || nzchar(system_file_cached(package = pkg))

  if (is.null(version)) {
    return(installed)
  }

  if (!is.character(version) && !inherits(version, "numeric_version")) {
    # Avoid https://bugs.r-project.org/show_bug.cgi?id=18548
    alert <- if (identical(Sys.getenv("TESTTHAT"), "true")) stop else warning
    alert("`version` must be a character string or a `package_version` or `numeric_version` object.")

    version <- numeric_version(sprintf("%0.9g", version))
  }

  installed && isTRUE(get_package_version(pkg) >= version)
}

get_package_version <- function(pkg) {
  # `utils::packageVersion()` can be slow, so first try the fast path of
  # checking if the package is already loaded.
  ns <- .getNamespace(pkg)
  if (is.null(ns)) {
    utils::packageVersion(pkg)
  } else {
    as.package_version(ns$.__NAMESPACE__.$spec[["version"]])
  }
}
is_ggplot2_350 <- function() {
  is_installed("ggplot2", version = "3.5.0")
}

scale_color_npg_custom <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_npg_custom(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "npg", palette = pal_npg_custom(palette, alpha), ...)
  }
}

p_final <- p + scale_color_npg_custom()

ggsave("umap_final.pdf", p_final, width = 12, height = 8)

p_final <- p_final + theme(aspect.ratio = 1)

############################################################################################

library(Seurat)
library(cowplot)
library(ggplot2)
library(grid)

calculate_mean_expression <- function(scRNA, gene_family) {
  all_genes <- rownames(scRNA)
  gene_names <- grep(paste0("^", gene_family), all_genes, value = TRUE)
  
  if (length(gene_names) == 0) {
    stop(paste("No genes found for gene family:", gene_family))
  }

  expression_data <- FetchData(scRNA, vars = gene_names)
  
  gene_names_no_na <- gene_names[apply(expression_data, 2, function(col) all(!is.na(col)))]

  if (length(gene_names_no_na) == 0) {
    stop(paste("No non-NA genes found for gene family:", gene_family))
  }

  mean_expression <- rowMeans(expression_data[, gene_names_no_na, drop = FALSE])
  return(mean_expression)
}

add_mean_expression_to_seurat <- function(scRNA, gene_family, new_col_name) {
  mean_expression <- calculate_mean_expression(scRNA, gene_family)
  scRNA <- AddMetaData(scRNA, metadata = mean_expression, col.name = new_col_name)
  return(scRNA)
}

create_feature_plot <- function(scRNA, feature, color, title) {
  p <- FeaturePlot(scRNA, features = feature, min.cutoff = "q10", max.cutoff = "q99", reduction = "umap.cca") + 
    scale_color_gradient(low = "white", high = color) + 
    theme_void() + 
    theme(legend.position = "none", plot.title = element_blank())
  p <- ggdraw(p) + draw_label(title, color = color, size = 16, fontface = "bold", 
                              hjust = 1, vjust = 0, x = 0.95, y = 0.05)
  return(p)
}

scRNA <- add_mean_expression_to_seurat(scRNA, "IGHG", "IGHG_mean")
scRNA <- add_mean_expression_to_seurat(scRNA, "IGHA", "IGHA_mean")
scRNA <- add_mean_expression_to_seurat(scRNA, "IGKV", "IGKV_mean")
scRNA <- add_mean_expression_to_seurat(scRNA, "IGLV", "IGLV_mean")
scRNA <- add_mean_expression_to_seurat(scRNA, "IGLC", "IGLC_mean")

p1 <- create_feature_plot(scRNA, "IGHM", "#F197D4", "IGHM")
p2 <- create_feature_plot(scRNA, "IGHD", "#85A6A1", "IGHD")
p3 <- create_feature_plot(scRNA, "IGHG_mean", "#556EA1", "IGHG")
p4 <- create_feature_plot(scRNA, "IGHA_mean", "#722955", "IGHA")
p5 <- create_feature_plot(scRNA, "IGKV_mean", "#5AD3AA", "IGKV")
p6 <- create_feature_plot(scRNA, "IGKC", "#AB448B", "IGKC")
p7 <- create_feature_plot(scRNA, "IGLV_mean", "#A8919D", "IGLV")
p8 <- create_feature_plot(scRNA, "IGLC_mean", "#5E7DE3", "IGLC")

combined_plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, rel_widths = c(1, 1, 1, 1))

combined_plot_with_border <- ggdraw(combined_plot) +
  annotation_custom(
    grob = roundrectGrob(
      x = 0.5, y = 0.5, width = 0.99, height = 0.99, r = unit(0.02, "npc"),
      gp = gpar(fill = NA, col = "grey", lty = "99", lwd = 1) 
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

for (x_pos in c(0.25, 0.5, 0.75)) {
  combined_plot_with_border <- combined_plot_with_border +
    draw_line(x = c(x_pos, x_pos), y = c(0.02, 0.98), color = "grey", lty = "77")
}

combined_plot_with_border <- combined_plot_with_border +
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", lty = "77") 

final_combined_plot <- plot_grid(p_final, combined_plot_with_border, ncol = 1, align = 'v', axis = 'tb', rel_heights = c(1, 0.5))

#ggsave("final_combined_plot.png", final_combined_plot, width = 10, height = 15, dpi = 300)
ggsave("p1_umap_final_combined_plot.pdf", final_combined_plot, width = 10, height = 15, dpi = 300)
```

## Fig1 C

```R
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

genes_list <- list(
  c("CD34", "IGLL1", "DNTT", "MKI67"),
  c("RAG1", "RAG2", "MS4A1"),
  c("MME", "CD72", "TCL1A"),
  c("PLD4", "IL4R", "IGHD", "IGHM"),
  c("IFITM1", "MX1", "S100A8", "S100A9", "ITGAX"),
  c("CD27", "AIM2"),
  c("MZB1", "SDC1", "JCHAIN")
)

celltype_levels <- c(
  "Resting Pro-B", "Cycling Pro-B", "Cycling Pre-B", "Resting Pre-B", 
  "Immature B", "S100A8/A9 high B", "Transitional B", "Naïve B", 
  "CD27+ Memory B", "Plasma"
)

all_genes <- unique(unlist(genes_list))
plot_data_all <- DotPlot(scRNA, features = all_genes)$data

global_min_expr <- min(plot_data_all$avg.exp.scaled)
global_max_expr <- max(plot_data_all$avg.exp.scaled)
global_min_pct <- min(plot_data_all$pct.exp)
global_max_pct <- max(plot_data_all$pct.exp)

dotplots <- lapply(genes_list, function(genes) {
  plot_data <- DotPlot(scRNA, features = genes)$data
  plot_data$features.plot <- factor(plot_data$features.plot, levels = genes)
  plot_data$id <- factor(plot_data$id, levels = celltype_levels)
  ggplot(plot_data, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(global_min_expr, global_max_expr),
                          guide = guide_colorbar(title.position = "left", title.theme = element_text(hjust = 0.5),
                                                 barwidth = 10, barheight = 0.5)) +
    scale_size(range = c(0.1, 6), limits = c(global_min_pct, global_max_pct), 
               guide = guide_legend(title.position = "left", title.hjust = 0.5)) +
    theme(
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
      axis.ticks.x = element_line(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "lightgrey"),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(x = "", y = "", color = "Average Expression", size = "Percent Expressed")
})

dotplots <- lapply(seq_along(dotplots), function(i) {
  if (i == 1) {
    dotplots[[i]] + theme(axis.title.y = element_blank())
  } else {
    dotplots[[i]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
  }
})

dotplots[[1]] <- dotplots[[1]] + theme(legend.position = "bottom", legend.direction = "horizontal")

combined_dotplot <- wrap_plots(dotplots, ncol = length(dotplots)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave("p2_marker_combined_dotplot.pdf", plot = combined_dotplot, width = 12, height = 6)
```

## Fig1 D

```R
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(SeuratWrappers)
library(viridis)

setwd("/home/data/tmp_data/zyh-BCR/CCA/plot/p1_umap/psedotime/")
scRNA <- readRDS("/home/data/tmp_data/zyh-BCR/CCA/CCA_annotated_0826.rds")
umap_cca_embeddings <- Embeddings(scRNA, reduction = "umap.cca")
scRNA[["umap"]] <- CreateDimReducObject(embeddings = umap_cca_embeddings, key = "UMAP_")
scRNA[["umap.cca"]] <- NULL
scRNA
 
cds <- as.cell_data_set(scRNA)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell.type") + ggtitle('cds.umap')

get_earliest_principal_node <- function(cds, time_bin=c('DEGC')){
  cell_ids <- which(colData(cds)[, "cell.type"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
    }

nodes_vec <- c(get_earliest_principal_node(cds,"Cycling Pre-B"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")

p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

ggsave(plot=p,file=paste0('monocle3.pdf'),width=14,height=8)
```

## Fig1 E

```R
setwd("/home/data/tmp_data/zyh-BCR/CCA/plot/p1_umap/propotion/")
library(Seurat)
library(ggplot2)
library(dplyr)

meta_data <- sce@meta.data

cell_type_order <- c("Resting Pro-B", "Cycling Pro-B", "Cycling Pre-B", "Resting Pre-B", 
                     "Immature B", "Transitional B", "S100A8/A9 high B", "Naïve B", 
                     "CD27+ Memory B", "Plasma")

meta_data$cell.type <- factor(meta_data$cell.type, levels = cell_type_order)

total_counts_group <- meta_data %>%
  group_by(group) %>%
  summarise(total_count = n())

total_counts_sample <- meta_data %>%
  group_by(orig.ident, group) %>%
  summarise(total_count = n())

cell_type_group <- meta_data %>%
  group_by(cell.type, group) %>%
  summarise(count = n()) %>%
  left_join(total_counts_group, by = "group") %>%
  mutate(proportion = count / total_count * 100) %>%
  mutate(orig.ident = paste(group)) %>%
  mutate(orig.ident = if_else(orig.ident == "Healthy", "HD", orig.ident)) %>%
  mutate(group = "Total")
  
cell_type_sample <- meta_data %>%
  group_by(cell.type, orig.ident, group) %>%
  summarise(count = n()) %>%
  left_join(total_counts_sample, by = c("orig.ident", "group")) %>%
  mutate(group = if_else(group == "Healthy", "HD", group)) %>%
  mutate(proportion = count / total_count * 100)

combined_data <- bind_rows(cell_type_group, cell_type_sample)

p <- ggplot(combined_data, aes(x = orig.ident, y = proportion, fill = cell.type)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.6) + 
  scale_fill_manual(values = c(
    "Resting Pro-B" = "#7E6148", "Cycling Pro-B" = "#91D1C2", "Cycling Pre-B" = "#3C5488",
    "Resting Pre-B" = "#8491B4", "Immature B" = "#F39B7F", "Transitional B" = "#4DBBD5",
    "S100A8/A9 high B" = "#FFDC91", "Naïve B" = "#00A087", "CD27+ Memory B" = "#EE4C97",
    "Plasma" = "#B09C85")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Fraction of cells (%)",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.1, "lines"),
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 12, hjust = 0.5, vjust = 1),
        strip.placement = "outside",
        panel.border = element_rect(fill = NA, color = "black")) +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x") 

p <- p + theme(
  strip.text.x = element_text(face = "bold", size = 12, hjust = 0.5, vjust = -1.5, margin = margin(b = 10)),
  strip.background = element_rect(fill = "white", color = "black", size = 1)
)

ggsave("proportion_stacked_combined.pdf", plot = p, height = 6, width = 8)
```

## Fig1 F

```R
# Import required packages
library(alakazam)
library(dplyr)
library(ggplot2)
library(shazam)
library(cowplot)
library(tidyr)
library(ggnewscale)
library(tibble)

db <- vdj_data

db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)

# Show new mutation frequency columns
db_obs %>% 
    select(sequence_id, starts_with("mu_freq")) %>%
    head(n=4)
    
db_obs$sample <- sub("^(.*)_[^_]+_[^_]+$", "\\1-1", db_obs$sequence_id)

head(db_obs$sample)

db_obs %>% 
    select(sequence_id, starts_with("mu_count_")) %>%
    head(n=4)
db_obs %>% 
    select(sequence_id, starts_with("mu_freq_")) %>%
    head(n=4)

db_obs <- db_obs %>%
  mutate(SHM_levels = case_when(
    mu_freq * 100 <= 1 ~ "Low",
    mu_freq * 100 > 1 & mu_freq * 100 <= 5 ~ "Median",
    mu_freq * 100 > 5 ~ "High"
  ))

head(db_obs$SHM_levels)

meta_data <- scRNA@meta.data

head(meta_data)

meta_data <- rownames_to_column(meta_data, var = "cell_name")
head(meta_data$cell_name)

db_obs <- db_obs %>%
  left_join(meta_data %>% select(cell_name, group), by = c("sample" = "cell_name"))

head(db_obs)

selected_cell_types <- c("Immature B", "Transitional B", "Naïve B", 
                         "CD27+ Memory B", "Plasma")

plot_df <- db_obs %>%
  filter(cell.type %in% selected_cell_types) %>%
  group_by(cell.type, SHM_levels, group) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cell.type, group) %>%
  mutate(sum_n = sum(n)) %>%
  ungroup() %>%
  mutate(percent = (n / sum_n) * 100)

plot_df <- plot_df %>%
  mutate(group = ifelse(group == "Healthy", "HD", group))

total_b_df <- db_obs %>%
  filter(cell.type %in% selected_cell_types) %>%
  group_by(SHM_levels, group) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(sum_n = sum(n)) %>%
  ungroup() %>%
  mutate(percent = (n / sum_n) * 100) %>%
  mutate(group = ifelse(group == "Healthy", "HD", group)) %>%
  mutate(cell.type = "Total B")

plot_df <- bind_rows(plot_df, total_b_df)

plot_df$SHM_levels <- factor(plot_df$SHM_levels, levels = c("Low", "Median", "High"))
plot_df$cell.type <- factor(plot_df$cell.type, levels = c(selected_cell_types, "Total B")) 

plot_df <- plot_df %>%
  mutate(SHM_levels_group = factor(paste(SHM_levels, group), levels = c("Low HD", "Low ITP", 
                                                                        "Median HD", "Median ITP", 
                                                                        "High HD", "High ITP")))
plot_df <- plot_df %>%
  mutate(cell_type_factor = as.numeric(factor(cell.type, levels = c(selected_cell_types, "Total B"))),
         cell_type_adj = ifelse(group == "HD",
                                cell_type_factor - 0.2, 
                                cell_type_factor + 0.2))

plot_df_median <- db_obs %>%
  filter(cell.type %in% selected_cell_types) %>%
  group_by(cell.type, group) %>% 
  summarise(median_SHM = mean(mu_freq)) %>%
  ungroup() %>%
  mutate(group = ifelse(group == "Healthy", "HD", group))

total_b_median <- db_obs %>%
  filter(cell.type %in% selected_cell_types) %>%
  group_by(group) %>%
  summarise(median_SHM = mean(mu_freq)) %>%
  mutate(cell.type = "Total B") %>%
  mutate(group = ifelse(group == "Healthy", "HD", group))

plot_df_median <- bind_rows(plot_df_median, total_b_median)

plot_df_median$cell.type <- factor(plot_df_median$cell.type, levels = c(selected_cell_types, "Total B"))
plot_df_median <- plot_df_median %>%
  mutate(cell_type_factor = as.numeric(factor(cell.type, levels = c(selected_cell_types, "Total B"))),
         cell_type_adj = ifelse(group == "HD",
                                cell_type_factor - 0.2,
                                cell_type_factor + 0.2))

colors_hd <- c("Low HD" = "#d1e5f0", "Median HD" = "#6BABD2", "High HD" = "#5987BD")
colors_itp <- c("Low ITP" = "#fde0dd", "Median ITP" = "#D67E6F", "High ITP" = "#BF305A")

x_breaks <- 1:length(c(selected_cell_types, "Total B"))
x_labels <- c(selected_cell_types, "Total B")

common_limits <- range(plot_df$cell_type_adj, plot_df_median$cell_type_adj)
common_limits[1] <- common_limits[1] - 0.3 
common_limits[2] <- common_limits[2] + 0.3 

p2 <- ggplot() +
  geom_bar(data = filter(plot_df, group == "HD"),
           aes(x = cell_type_adj, y = percent, fill = SHM_levels_group),
           stat = "identity", position = position_stack(), width = 0.4) +
  scale_fill_manual(name = "SHM Levels (HD)", values = colors_hd[levels(plot_df$SHM_levels_group)[levels(plot_df$SHM_levels_group) %in% names(colors_hd)]]) +
  new_scale_fill() +
  geom_bar(data = filter(plot_df, group == "ITP"),
           aes(x = cell_type_adj, y = percent, fill = SHM_levels_group), 
           stat = "identity", position = position_stack(), width = 0.4) +
  scale_fill_manual(name = "SHM Levels (ITP)", values = colors_itp[levels(plot_df$SHM_levels_group)[levels(plot_df$SHM_levels_group) %in% names(colors_itp)]]) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = common_limits) +
  theme_minimal() +
  labs(title = "", x = "Cell Type", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

p3 <- ggplot(plot_df_median, aes(x = cell_type_adj, y = median_SHM, color = group)) +
  geom_point(shape = 16, stroke = 0, size = 3) +
  scale_color_manual(values = c("HD" = "#2166ac", "ITP" = "#b2182b")) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = common_limits) +
  scale_y_continuous(breaks = seq(0, 0.08, by = 0.02), limits = c(0, 0.08)) +
  labs(x = NULL, y = "SHM rate", color = "Group") +  
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "right"
  )

combined_plot <- plot_grid(p3, p2, ncol = 1, align = "v", rel_heights = c(0.3, 0.7))
ggsave(filename = "combined_BCR_diversity_plot_with_legend.pdf", plot = combined_plot, width = 8, height = 6)
```

