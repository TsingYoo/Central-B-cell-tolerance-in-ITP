# FigS2

## FigS2 F

```R
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(grid)

meta_data <- seurat@meta.data

cell_types <- c("Immature B", "Transitional B", "Naïve B", "CD27+ Memory B", "Plasma")

for (cell_type in cell_types) {
  filtered_meta_data <- meta_data %>%
    filter(!is.na(CTgene) & CTgene != "" & cell.type == cell_type) 
  sample_counts <- filtered_meta_data %>%
    group_by(orig.ident) %>%
    summarise(cell_count = n()) %>%
    deframe()

  overlap_matrix <- matrix(0, nrow = 10, ncol = 10, dimnames = list(samples_name_new, samples_name_new))

  sample_ctgenes <- filtered_meta_data %>%
    group_by(orig.ident) %>%
    summarise(CTgene_set = list(unique(CTgene))) %>%
    deframe()

  for (sample_1 in samples_name_new) {
    for (sample_2 in samples_name_new) {
      ctgenes_1 <- sample_ctgenes[[sample_1]]
      ctgenes_2 <- sample_ctgenes[[sample_2]]
      overlap_count <- length(intersect(ctgenes_1, ctgenes_2))
      overlap_percentage <- (overlap_count / length(ctgenes_1)) * 100
      overlap_matrix[sample_1, sample_2] <- overlap_percentage
    }
  }

  col_fun <- colorRamp2(c(min(overlap_matrix), median(overlap_matrix), max(overlap_matrix[overlap_matrix < 100])), 
                        c("#3399CC", "#FCA758", "#CB5A5A"))

  ht <- Heatmap(overlap_matrix,
                 name = "Overlap Percentage",
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_names = TRUE,
                 show_column_names = TRUE,
                 row_names_side = "left",
                 column_names_side = "bottom",
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if (i == j) {
                     grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = "white", col = NA))
                     grid.text(sprintf("n = %d", sample_counts[samples_name_new[i]]), x, y)
                   } else {
                     grid.text(sprintf("%.2f%%", overlap_matrix[i, j]), x, y)
                   }
                 },
                 heatmap_legend_param = list(title = "Percentage (%)",
                                             legend_direction = "vertical",
                                             legend_height = unit(4, "cm"),
                                             title_position = "topleft",
                                             at = c(0, 5, 10),
                                             border = NA),
                 width = unit(20, "cm"),
                 height = unit(20, "cm"))

  pdf(paste0("BCR_overlap_heatmap_complex_", gsub(" ", "_", cell_type), ".pdf"), width = 10.5, height = 9.5)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

  decorate_heatmap_body("Overlap Percentage", {
    grid.lines(x = c(0.5, 0.5), y = c(0, 1))
    grid.lines(x = c(0, 1), y = c(0.5, 0.5))
    grid.text("HD", x = unit(0.25, "npc"), y = unit(-0.1, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.text("ITP", x = unit(0.75, "npc"), y = unit(-0.1, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
    grid.rect(gp = gpar(col = "black", fill = NA, lwd = 2))
  })

  dev.off()
}
```

## FigS2 J

```R
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

diversity_data <- clonalDiversity(combined.BCR, 
                                  cloneCall = "CTgene", 
                                  group.by = "sample", 
                                  x.axis = "ID", 
                                  exportTable = TRUE)

head(diversity_data)

diversity_data_long <- melt(diversity_data, id.vars = c("sample", "ID"), 
                            measure.vars = c("shannon", "inv.simpson", "norm.entropy", 
                                             "gini.simpson", "chao1", "ACE"),
                            variable.name = "Method", value.name = "Value")

p_values <- diversity_data_long %>%
  group_by(Method) %>%
  summarize(p_value = wilcox.test(Value[ID == "HD"], Value[ID == "ITP"])$p.value)

print(p_values)

p_top <- ggplot(diversity_data_long %>% filter(Method %in% c("shannon", "inv.simpson", "norm.entropy")), 
                aes(x = ID, y = Value, fill = ID)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Method, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") 

p_bottom <- ggplot(diversity_data_long %>% filter(Method %in% c("gini.simpson", "chao1", "ACE")), 
                   aes(x = ID, y = Value, fill = ID)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Method, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

p_top <- p_top + stat_compare_means(method = "wilcox.test", 
                                    label = "p.signif", 
                                    comparisons = list(c("HD", "ITP")))

p_bottom <- p_bottom + stat_compare_means(method = "wilcox.test", 
                                          label = "p.signif", 
                                          comparisons = list(c("HD", "ITP")))

legend <- get_legend(
  ggplot(diversity_data_long, aes(x = ID, y = Value, fill = ID)) +
    geom_boxplot() +
    scale_fill_manual(values = c("HD" = "#6A92A9", "ITP" = "#BFA3B2")) +
    theme_minimal() +
    labs(fill = "Group") +
    theme(legend.position = "right")
)

combined_plot <- plot_grid(p_top, p_bottom, ncol = 1, rel_heights = c(1, 1))
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(0.85, 0.15))
final_plot <- ggdraw(final_plot) +
  draw_label("Group", x = 0.5, y = 0, vjust = 1, angle = 0, fontface = "plain", size = 10) +
  draw_label("Index Value", x = -0.02, y = 0.5, vjust = 1, angle = 90, fontface = "plain", size = 10)

final_plot <- final_plot + theme(plot.margin = margin(5, 0, 10, 10))

ggsave(filename = "clonalDiversity.pdf", plot = final_plot, width = 6, height = 6)
```
