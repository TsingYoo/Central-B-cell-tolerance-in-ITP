# FigS3

## FigS3 F

```R
library(ggplot2)
library(dplyr)
library(scales)
library(Cairo)

cell_type_order <- c("CD27+ Memory B", "Plasma")

object <- db_obs
object$group <- recode(object$group, "Healthy" = "HD")
object <- subset(object, cell.type %in% cell_type_order)

object <- object %>%
  left_join(
    seurat@meta.data %>% 
      rownames_to_column(var = "sample") %>% 
      dplyr::select(sample, cloneSize),
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

generate_locus_plot <- function(cdr3aa, locus_type = NULL, output_file) {
  locus_data <- if (!is.null(locus_type)) {
    cdr3aa %>% filter(locus == locus_type)
  } else {
    cdr3aa
  }
  median_value <- median(locus_data$cdr3_aa_charge, na.rm = TRUE)
  print(median_value)
  
  if (is.null(locus_type)) {
    locus_data <- locus_data %>%
      mutate(charge_category = factor(case_when(
        cdr3_aa_charge >= 0 & cdr3_aa_charge <= 1 ~ "0",
        cdr3_aa_charge >= 1  ~ "≥1"
      ), levels = c("0", "≥1")))
  } else if (locus_type == "IGH") {
    locus_data <- locus_data %>%
      mutate(charge_category = factor(case_when(
        cdr3_aa_charge >= 0 & cdr3_aa_charge < 1 ~ "0",
        cdr3_aa_charge >= 1 & cdr3_aa_charge < 2 ~ "1",
        cdr3_aa_charge >= 2 & cdr3_aa_charge < 3 ~ "2",
        cdr3_aa_charge >= 3 ~ "≥3"
      ), levels = c("0", "1", "2", "≥3")))
  } else if (locus_type == "IGK") {
    locus_data <- locus_data %>%
      mutate(charge_category = factor(case_when(
        cdr3_aa_charge >= 0 & cdr3_aa_charge < 1 ~ "0",
        cdr3_aa_charge >= 1 & cdr3_aa_charge < 2 ~ "1",
        cdr3_aa_charge >= 2  ~ "≥2"
      ), levels = c("0", "1", "≥2")))
  } else if (locus_type == "IGL") {
    locus_data <- locus_data %>%
      mutate(charge_category = factor(case_when(
        cdr3_aa_charge >= 0 & cdr3_aa_charge < 1 ~ "0",
        cdr3_aa_charge >= 1 & cdr3_aa_charge < 2 ~ "1",
        cdr3_aa_charge >= 2 ~ "≥2"
      ), levels = c("0", "1", "≥2")))
  }
  
  locus_data <- locus_data %>%
    mutate(cloneSizeCategory = factor(cloneSizeCategory, levels = c("Other", "Large")))
  
  charge_proportion <- locus_data %>%
    group_by(cloneSizeCategory, charge_category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(cloneSizeCategory) %>%
    mutate(proportion = count / sum(count))
  
  max_categories = 4
  n_categories <- length(unique(charge_proportion$charge_category))
  width_ratio <- n_categories / max_categories
  
  p <- ggplot(charge_proportion, aes(x = charge_category, y = proportion, fill = cloneSizeCategory)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85 * width_ratio, preserve = "single"), color = "black", width = 0.65 * width_ratio) +
    geom_text(aes(label = sprintf("%.2f%%", proportion * 100)), 
      position = position_dodge2(width = 0.85 * width_ratio, preserve = "single"), 
      vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c("Other" = "#90B28D", "Large" = "#DDAEAA")) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal(base_line_size = 0) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Charge Category",
      y = "Proportion (%)",
      fill = "Group",
      title = ifelse(is.null(locus_type), "All Loci CDR3 Positive Charges", paste(locus_type, "CDR3 Positive Charges"))
    )
  
  CairoPDF(file = output_file, width = 6, height = 4)
  print(p)
  dev.off()
  
  return(p)
}

p_IGH <- generate_locus_plot(cdr3aa, "IGH", "charge_pos_IGH.pdf")
p_IGL <- generate_locus_plot(cdr3aa, "IGL", "charge_pos_IGL.pdf")
p_IGK <- generate_locus_plot(cdr3aa, "IGK", "charge_pos_IGK.pdf")

p_all <- generate_locus_plot(cdr3aa, output_file = "charge_pos_all_loci.pdf")
```



