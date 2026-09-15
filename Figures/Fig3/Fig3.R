library(tidyverse)

# ---- Folder containing the files ----
folder <- "~/Desktop/DENV_EVO_2025/HyPhy_results/RELAX_results"  # <-- change if needed
files <- file.path(folder, paste0("DENV", 1:4, "_combined_summary.csv"))

# ---- Desired protein order ----
protein_order  <- c("C","PrM","E","NS1","NS2A","NS2B","NS3","NS4A","NS4B","NS5")
serotype_order <- paste0("DENV", 1:4)

# ---- Per-serotype comparison groups to KEEP ----
keep_groups <- list(
  DENV1 = c("1III", "1IV", "1V", "1VII"),
  DENV2 = c("2II", "2III", "2V"),
  DENV3 = c("3II", "3III"),
  DENV4 = c("4II")
)

# ---- Read + combine ----
df_all <- purrr::map_dfr(files, function(f) {
  readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(serotype = stringr::str_extract(basename(f), "DENV[1-4]"))
})

# ---- Clean / transform ----
df_plot <- df_all %>%
  mutate(
    serotype = factor(serotype, levels = serotype_order),
    comparison_group = as.character(comparison_group),
    
    # protein = last chunk after underscore (e.g., capsid_protein_C -> C)
    protein_raw = stringr::str_extract(gene, "[^_]+$"),
    protein = case_when(tolower(protein_raw) == "prm" ~ "PrM", TRUE ~ protein_raw),
    
    RELAX_K = as.numeric(RELAX_K),
    RELAX_overall_pval = as.numeric(RELAX_overall_pval),
    
    sig = RELAX_overall_pval < 0.01,
    direction = if_else(RELAX_K >= 1, "Intensified (k>1)", "Relaxed (k<1)"),
    
    # magnitude as fold-change away from 1 (symmetric for k<1 vs k>1)
    fold_change_from_1 = if_else(RELAX_K >= 1, RELAX_K, 1 / RELAX_K),
    
    # 4-shade class for fill
    sig_class = if_else(sig, "Significant (p < 0.01)", "Not significant (p ≥ 0.01)"),
    dir_sig = interaction(direction, sig_class, sep = " • ")
  ) %>%
  filter(
    protein %in% protein_order,
    is.finite(RELAX_K), RELAX_K > 0
  ) %>%
  # exclude per-serotype reference groups
  filter(
    !(
      (serotype == "DENV1" & comparison_group %in% c("1_I", "1I")) |
        (serotype == "DENV2" & comparison_group %in% c("2_I", "2I")) |
        (serotype == "DENV3" & comparison_group %in% c("3_I", "3I")) |
        (serotype == "DENV4" & comparison_group %in% c("4_I", "4I"))
    )
  ) %>%
  # keep only specified x-groups for each serotype
  filter(
    case_when(
      serotype == "DENV1" ~ comparison_group %in% keep_groups$DENV1,
      serotype == "DENV2" ~ comparison_group %in% keep_groups$DENV2,
      serotype == "DENV3" ~ comparison_group %in% keep_groups$DENV3,
      serotype == "DENV4" ~ comparison_group %in% keep_groups$DENV4,
      TRUE ~ TRUE
    )
  ) %>%
  mutate(
    protein = factor(protein, levels = rev(protein_order)),
    comparison_group = factor(comparison_group),
    
    # control legend order (significant first within each direction)
    dir_sig = factor(
      dir_sig,
      levels = c(
        "Intensified (k>1) • Significant (p < 0.01)",
        "Intensified (k>1) • Not significant (p ≥ 0.01)",
        "Relaxed (k<1) • Significant (p < 0.01)",
        "Relaxed (k<1) • Not significant (p ≥ 0.01)"
      )
    )
  )

# ---- Wrap legend labels so nothing gets clipped ----
legend_labels <- c(
  "Intensified (k>1) • Significant (p < 0.01)"     = "Intensified (k > 1)\nSignificant (p < 0.01)",
  "Intensified (k>1) • Not significant (p ≥ 0.01)" = "Intensified (k > 1)\nNot significant (p ≥ 0.01)",
  "Relaxed (k<1) • Significant (p < 0.01)"         = "Relaxed (k < 1)\nSignificant (p < 0.01)",
  "Relaxed (k<1) • Not significant (p ≥ 0.01)"     = "Relaxed (k < 1)\nNot significant (p ≥ 0.01)"
)

# ---- Plot (faceted left-to-right, free x per facet) ----
p <- ggplot(df_plot, aes(x = comparison_group, y = protein)) +
  
  # Non-significant: light fill + thin grey outline, still visible
  geom_point(
    data = df_plot %>% filter(!sig),
    aes(size = fold_change_from_1, fill = dir_sig),
    shape = 21,
    color = "grey45",
    stroke = 0.35,
    alpha = 0.55
  ) +
  
  # Significant: darker fill + thicker black outline
  geom_point(
    data = df_plot %>% filter(sig),
    aes(size = fold_change_from_1, fill = dir_sig),
    shape = 21,
    color = "black",
    stroke = 1.0,
    alpha = 1
  ) +
  
  facet_wrap(~ serotype, nrow = 1, scales = "free_x") +
  
  scale_fill_manual(
    name = "Selection direction (color)\nand significance (shade)",
    values = c(
      "Intensified (k>1) • Significant (p < 0.01)"      = "#B30000",
      "Intensified (k>1) • Not significant (p ≥ 0.01)"  = "#F8B6B6",
      "Relaxed (k<1) • Significant (p < 0.01)"          = "#084594",
      "Relaxed (k<1) • Not significant (p ≥ 0.01)"      = "#BFD9F2"
    ),
    labels = legend_labels
  ) +
  
  scale_size_continuous(
    name = "Fold-change from k = 1",
    range = c(3, 18),
    breaks = c(1, 2, 5, 10, 50),
    labels = c("1×", "2×", "5×", "10×", "50x")
  ) +
  
  scale_y_discrete(drop = FALSE) +
  labs(x = "Comparison group", y = "Protein") +
  
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    size = guide_legend(nrow = 1)
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 0, size = 16),
    axis.text.y = element_text(size = 16),
    
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    
    strip.text = element_text(size = 18, face = "bold"),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title.align = 0,
    legend.text.align = 0,
    
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.y = unit(0.25, "cm"),
    
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10),
    
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "black")
  )

print(p)

# ---- Minimum-info version for manuscript figure ----
legend_labels_min <- c(
  "Intensified (k>1) • Significant (p < 0.01)"     = "Intensified, significant",
  "Intensified (k>1) • Not significant (p ≥ 0.01)" = "Intensified, not significant",
  "Relaxed (k<1) • Significant (p < 0.01)"         = "Relaxed, significant",
  "Relaxed (k<1) • Not significant (p ≥ 0.01)"     = "Relaxed, not significant"
)

p_min <- p +
  scale_fill_manual(
    name = "Selection intensity",
    values = c(
      "Intensified (k>1) • Significant (p < 0.01)"      = "#B30000",
      "Intensified (k>1) • Not significant (p ≥ 0.01)"  = "#F8B6B6",
      "Relaxed (k<1) • Significant (p < 0.01)"          = "#084594",
      "Relaxed (k<1) • Not significant (p ≥ 0.01)"      = "#BFD9F2"
    ),
    labels = legend_labels_min
  )

# ---- Save detailed version ----
ggsave(
  "~/Desktop/DENV_EVO_2025/HyPhy_results/RELAX_results/RELAX_bubbleplot_detailed.svg",
  p,
  width = 12,
  height = 8,
  device = "svg",
  bg = "white"
)

ggsave(
  "~/Desktop/DENV_EVO_2025/HyPhy_results/RELAX_results/RELAX_bubbleplot_detailed.png",
  p,
  width = 12,
  height = 8,
  device = "png",
  bg = "white"
)

# ---- Save minimum-info version ----
ggsave(
  "~/Desktop/DENV_EVO_2025/HyPhy_results/RELAX_results/RELAX_bubbleplot_minimal.svg",
  p_min,
  width = 14,
  height = 8,
  device = "svg",
  bg = "white"
)

ggsave(
  "~/Desktop/DENV_EVO_2025/HyPhy_results/RELAX_results/RELAX_bubbleplot_minimal.png",
  p_min,
  width = 14,
  height = 8,
  device = "png",
  bg = "white"
)
