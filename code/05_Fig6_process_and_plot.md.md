# Fig6

## Fig6 A

```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

gene_expression <- FetchData(sce, vars = c("RAG1", "RAG2"))
cell_type_info <- sce@meta.data$cell.type
group_info <- sce@meta.data$group
gene_expression$cell.type <- cell_type_info
gene_expression$group <- group_info
gene_expression$group <- gsub("Healthy", "HD", gene_expression$group)

batch1_cell_types <- c("Resting Pro-B", "Cycling Pro-B", "Cycling Pre-B", "Resting Pre-B", 
                       "Immature B", "Transitional B", "Naïve B", 
                       "CD27+ Memory B", "Plasma")

batch2_cell_types <- c("Immature B", "Transitional B", "Naïve B", 
                       "CD27+ Memory B", "Plasma")

filtered_gene_expression <- gene_expression %>% filter(group %in% c("ITP", "HD"))

plot_and_save <- function(cell_types, filename) {
  filtered_data <- filtered_gene_expression %>% filter(cell.type %in% cell_types)
  filtered_data$cell.type <- factor(filtered_data$cell.type, levels = cell_types)
  long_gene_expression <- filtered_data %>%
    pivot_longer(cols = c("RAG1", "RAG2"), names_to = "gene", values_to = "expression")

  mean_expression <- long_gene_expression %>%
    group_by(gene, cell.type, group) %>%
    summarise(
      mean_expression = mean(expression),
      se = sd(expression) / sqrt(n()),
      .groups = 'drop'
    )

  calculate_p_value <- function(data, gene) {
    tryCatch({
      t.test(expression ~ group, data = data)$p.value
    }, error = function(e) {
      NA
    })
  }

  p_values <- long_gene_expression %>%
    group_by(gene, cell.type) %>%
    summarise(
      p.value = calculate_p_value(cur_data(), gene),
      group1 = "ITP",
      group2 = "HD", 
      .groups = 'drop'
    )

  p_values <- p_values %>%
    mutate(signif = case_when(
      p.value <= 0.0001 ~ "****",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ))

  custom_colors <- c("HD" = "#6A92A9", "ITP" = "#BFA3B2")

  bar_plot <- ggplot(mean_expression, aes(x = cell.type, y = mean_expression, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(ymin = mean_expression - se, ymax = mean_expression + se),
                  position = position_dodge(width = 0.9), width = 0.25) +
    facet_wrap(~gene, scales = "free") +
    theme_minimal() +
    ggtitle("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(size = 14, face = "bold")) +
    labs(x = "", y = "Mean Expression", fill = "Group") +
    scale_fill_manual(values = custom_colors)

  bar_plot <- bar_plot +
    stat_pvalue_manual(
      data = p_values, 
      x = "cell.type", 
      y.position = max(mean_expression$mean_expression) * 1.05,
      label = "signif", 
      group1 = "group1",
      group2 = "group2",
      position = position_dodge(0.9),
      inherit.aes = FALSE
    )

  ggsave(filename, plot = bar_plot, width = 12, height = 8)
}

plot_and_save(batch1_cell_types, "RAG1_RAG2_expression_batch1.pdf")
plot_and_save(batch2_cell_types, "RAG1_RAG2_expression_batch2.pdf")

print(p_values)
```

