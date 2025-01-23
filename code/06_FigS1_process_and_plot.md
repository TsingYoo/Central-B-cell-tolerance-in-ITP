# FigS1

## FigS1 A

```R
library(Seurat)
library(dplyr)
library(ggplot2)

total.markers <- FindAllMarkers(scRNA,                               
                                 only.pos = TRUE,                               
                                 min.pct = 0.25,                               
                                 logfc.threshold = 0.25)

total.markers <- total.markers %>% filter(!startsWith(gene, "ENSG"))

cluster_colors <- c(
  "Resting Pro-B" = "#7E6148", 
  "Cycling Pro-B" = "#91D1C2", 
  "Cycling Pre-B" = "#3C5488",
  "Resting Pre-B" = "#8491B4", 
  "Immature B" = "#F39B7F", 
  "S100A8/A9 high B" = "#FFDC91", 
  "Transitional B" = "#4DBBD5",
  "Naïve B" = "#00A087", 
  "CD27+ Memory B" = "#EE4C97",
  "Plasma" = "#B09C85"
)

total.markers$cluster <- factor(total.markers$cluster, levels = names(cluster_colors))

top20 <- total.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster)

all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)

scRNA@meta.data$cell.type <- factor(scRNA@meta.data$cell.type, levels = names(cluster_colors))
Idents(scRNA) <- scRNA@meta.data$cell.type

p <- DoHeatmap(scRNA, features = top20$gene, group.colors = cluster_colors, label = FALSE) +
  scale_color_manual(values = cluster_colors, name = "Cell Type") +
  scale_fill_gradientn(colors = c("#509BD1", "white", "#E26E65"), 
                       limits = c(-2, 2),
                       oob = scales::squish,
                       name = "Expression") +
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1))
  )+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(3, "lines")
  )

ggsave("heatmap_marker.pdf", plot = p, width = 18, height = 20)
message("Heatmap saved as heatmap_marker.pdf")
```

## FigS1 B

```R
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(stringr)
library(DOSE)
library(forcats)

Idents(scRNA) <- "cell.type"
S100A8_A9_cells <- WhichCells(scRNA, idents = "S100A8/A9 high B")
other_cells <- WhichCells(scRNA, idents = setdiff(levels(Idents(scRNA)), "S100A8/A9 high B"))
scRNA$comparison_group <- "Other Clusters"
scRNA$comparison_group[S100A8_A9_cells] <- "S100A8/A9 high B"
Idents(scRNA) <- "comparison_group"
table(scRNA$comparison_group)

diff <- FindMarkers(scRNA, ident.1 = "S100A8/A9 high B", ident.2 = "Other Clusters",
                    logfc.threshold = 0.25, min.pct = 0.25)
diff <- diff %>% filter(p_val < 0.05)

if (nrow(diff) == 0) {
  message("No significant markers found for S100A8/A9 high B vs Other Clusters")
} else {
  diff$gene <- rownames(diff)
  diff$group <- ifelse(diff$avg_log2FC > 0, "up", "down")
  
  Gene_ID <- bitr(diff$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  data <- merge(Gene_ID, data.frame(gene = diff$gene, group = diff$group), by.x = 'SYMBOL', by.y = 'gene')
  
  diff_KEGG <- compareCluster(ENTREZID ~ group,
                              data = data,
                              fun = "enrichKEGG",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.01,
                              qvalueCutoff = 0.01,
                              organism = "hsa")
  
  if (is.null(diff_KEGG@compareClusterResult) || nrow(diff_KEGG@compareClusterResult) == 0) {
    message("No significant enrichment found for S100A8/A9 high B vs Other Clusters")
  } else {
    diff_KEGG <- setReadable(diff_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
    diff_KEGG <- diff_KEGG@compareClusterResult
    
    diff_KEGG <- diff_KEGG %>%
      mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
    
    diff_KEGG_up <- diff_KEGG %>% filter(group == "up") %>% top_n(n = 10, wt = -qvalue)
    diff_KEGG_down <- diff_KEGG %>% filter(group == "down") %>% top_n(n = 10, wt = -qvalue)
    
    diff_KEGG_up$Description <- str_wrap(diff_KEGG_up$Description, width = 30)
    diff_KEGG_down$Description <- str_wrap(diff_KEGG_down$Description, width = 30)
    
    if (nrow(diff_KEGG_up) > 0) {
      p_up <- ggplot(diff_KEGG_up, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment), size = Count, color = -log10(qvalue))) +
        geom_point(alpha = 0.7) +
        scale_size_continuous(range = c(2, 8)) +
        scale_color_gradient(low = "#6EA6A5", high = "#8987B7") +
        theme_classic() +
        theme(axis.title.y = element_text(colour = 'black', size = 12),
              axis.line = element_line(colour = 'black', linewidth = 0.5),
              axis.text.x = element_text(colour = 'black', size = 10),
              axis.ticks.x = element_line(colour = 'black'),
              axis.title.x = element_text(colour = 'black', size = 12),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.text.y = element_text(size = 8),
              legend.position = "right") +
        guides(size = guide_legend(order = 1),
               color = guide_colorbar(order = 2)) +
        labs(title = "KEGG Enrichment (Up): S100A8/A9 high B vs Other Clusters",
             y = "", x = "Fold Enrichment",
             size = "Gene Count", color = "-log10(qvalue)")
      
      ggsave("KEGG_S100A8_A9_vs_Other_Clusters_up_bubble.pdf", plot = p_up, width = 6, height = 5)
      message("KEGG enrichment bubble plot (Up) saved for S100A8/A9 high B vs Other Clusters")
    }
    
    if (nrow(diff_KEGG_down) > 0) {
      p_down <- ggplot(diff_KEGG_down, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment), size = Count, color = -log10(qvalue))) +
        geom_point(alpha = 0.7) +
        scale_size_continuous(range = c(2, 8)) +
        scale_color_gradient(low = "#6EA6A5", high = "#8987B7") +
        theme_classic() +
        theme(axis.title.y = element_text(colour = 'black', size = 12),
              axis.line = element_line(colour = 'black', linewidth = 0.5),
              axis.text.x = element_text(colour = 'black', size = 10),
              axis.ticks.x = element_line(colour = 'black'),
              axis.title.x = element_text(colour = 'black', size = 12),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.text.y = element_text(size = 8),
              legend.position = "right") +
        guides(size = guide_legend(order = 1),
               color = guide_colorbar(order = 2)) +
        labs(title = "KEGG Enrichment (Down): S100A8/A9 high B vs Other Clusters",
             y = "", x = "Fold Enrichment",
             size = "Gene Count", color = "-log10(qvalue)")
      
      ggsave("KEGG_S100A8_A9_vs_Other_Clusters_down_bubble.pdf", plot = p_down, width = 6, height = 5)
      message("KEGG enrichment bubble plot (Down) saved for S100A8/A9 high B vs Other Clusters")
    }
  }
}
```

## FigS1 C

```R
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(stringr)
library(DOSE)
library(forcats)

Idents(scRNA) <- "cell.type"
S100A8_A9_cells <- WhichCells(scRNA, idents = "S100A8/A9 high B")
other_cells <- WhichCells(scRNA, idents = setdiff(levels(Idents(scRNA)), "S100A8/A9 high B"))
scRNA$comparison_group <- "Other Clusters"
scRNA$comparison_group[S100A8_A9_cells] <- "S100A8/A9 high B"
Idents(scRNA) <- "comparison_group"
table(scRNA$comparison_group)

diff <- FindMarkers(scRNA, ident.1 = "S100A8/A9 high B", ident.2 = "Other Clusters",
                    logfc.threshold = 0.25, min.pct = 0.25)
diff <- diff %>% filter(p_val < 0.05)

if (nrow(diff) == 0) {
  message("No significant markers found for S100A8/A9 high B vs Other Clusters")
} else {
  diff$gene <- rownames(diff)
  diff$group <- ifelse(diff$avg_log2FC > 0, "up", "down")
  
  Gene_ID <- bitr(diff$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  data <- merge(Gene_ID, data.frame(gene = diff$gene, group = diff$group), by.x = 'SYMBOL', by.y = 'gene')
  
  diff_GO <- compareCluster(
    ENTREZID ~ group, 
    data = data, 
    fun = "enrichGO", 
    OrgDb = "org.Hs.eg.db", 
    ont = "BP",
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
  )
  
  if (!is.null(diff_GO@compareClusterResult) && nrow(diff_GO@compareClusterResult) > 0) {
    diff_GO <- setReadable(diff_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
    diff_GO <- diff_GO@compareClusterResult
    
    diff_GO$geneID <- sapply(strsplit(diff_GO$geneID, "/"), function(x) paste(x[1:min(5, length(x))], collapse = "/"))
    diff_GO_up <- diff_GO %>% filter(group == "up") %>% top_n(n = 10, wt = -qvalue)
    diff_GO_down <- diff_GO %>% filter(group == "down") %>% top_n(n = 10, wt = -qvalue)
    
    if (nrow(diff_GO_up) > 0) {
      diff_GO_up$Description <- factor(diff_GO_up$Description, levels = rev(diff_GO_up$Description))
      p_up <- ggplot(diff_GO_up, aes(x = -log10(qvalue), y = Description)) +
        geom_bar(stat = "identity", width = 0.4, fill = "#CB5640") +
        geom_text(aes(x = 0.1, y = Description, label = Description), size = 3, hjust = 0) +
        geom_text(aes(x = 0.1, y = Description, label = geneID), 
                  color = "#CB5640", size = 3, fontface = 'italic', hjust = 0, vjust = 2.3) +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_text(colour = 'black', size = 12),
              axis.line = element_line(colour = 'black', linewidth = 0.5),
              axis.text.x = element_text(colour = 'black', size = 10),
              axis.ticks.x = element_line(colour = 'black'),
              axis.title.x = element_text(colour = 'black', size = 12),
              plot.title = element_text(size = 10, hjust = 0.5),
              legend.position = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        labs(title = "GO Enrichment (Up): S100A8/A9 high B vs Other Clusters",
             y = "Enriched Pathways")
      
      ggsave("GO_S100A8_A9_vs_Other_Clusters_up_bar.pdf", plot = p_up, width = 4.1, height = 5)
      message("GO enrichment plot (Up) saved for S100A8/A9 high B vs Other Clusters")
    }
    
    if (nrow(diff_GO_down) > 0) {
      diff_GO_down$Description <- factor(diff_GO_down$Description, levels = rev(diff_GO_down$Description))
      p_down <- ggplot(diff_GO_down, aes(x = -log10(qvalue), y = Description)) +
        geom_bar(stat = "identity", width = 0.4, fill = "#65B0C6") +
        geom_text(aes(x = 0.1, y = Description, label = Description), size = 3, hjust = 0) +
        geom_text(aes(x = 0.1, y = Description, label = geneID), 
                  color = "#65B0C6", size = 3, fontface = 'italic', hjust = 0, vjust = 2.3) +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_text(colour = 'black', size = 12),
              axis.line = element_line(colour = 'black', linewidth = 0.5),
              axis.text.x = element_text(colour = 'black', size = 10),
              axis.ticks.x = element_line(colour = 'black'),
              axis.title.x = element_text(colour = 'black', size = 12),
              plot.title = element_text(size = 10, hjust = 0.5), 
              legend.position = "none") +
        scale_x_continuous(expand = c(0, 0)) +
        labs(title = "GO Enrichment (Down): S100A8/A9 high B vs Other Clusters",
             y = "Enriched Pathways")
      
      ggsave("GO_S100A8_A9_vs_Other_Clusters_down_bar.pdf", plot = p_down, width = 4.1, height = 5)
      message("GO enrichment plot (Down) saved for S100A8/A9 high B vs Other Clusters")
    }
  }
}
```

