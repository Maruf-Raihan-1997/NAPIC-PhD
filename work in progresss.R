# =====================================================================
# Free amino acid release from Lemna minor duckweed after in vitro digestion
# Statistics (One-way ANOVA + Tukey HSD)
# Conference Poster Version (Portrait)
# =====================================================================

library(tidyverse)
library(broom)
library(multcompView)

# ---------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------

codes <- c("ND", "W", "Bn", "Br")

labels <- c(
  ND = "Undigested DWP",
  W = "10% DWP Water",
  Bn = "10% DWP  Banana",
  Br = "10% DWP  Bread"
)

cols <- c(
  ND = "grey40",
  W = "steelblue",
  Bn = "forestgreen",
  Br = "saddlebrown"
)

dat <- tribble(
  ~amino_acid, ~treatment, ~r1, ~r2, ~r3,
  
  "Glutamic acid", "ND", 0.44, 0.45, 0.48,
  "Glutamic acid", "W", 0.23, 0.21, 0.36,
  "Glutamic acid", "Bn", 0.05, 0.05, 0.08,
  "Glutamic acid", "Br", 0.04, 0.04, 0.05,
  
  "Aspartic acid", "ND", 0.96, 0.75, 0.79,
  "Aspartic acid", "W", 0.85, NA, 0.92,
  "Aspartic acid", "Bn", 0.86, 1.03, 0.79,
  "Aspartic acid", "Br", 0.99, NA, 1.27,
  
  "Leucine", "ND", 0.22, 0.16, 0.21,
  "Leucine", "W", 0.45, 0.45, 0.44,
  "Leucine", "Bn", 0.16, 0.20, 0.14,
  "Leucine", "Br", 0.17, 0.16, 0.17
) |>
  pivot_longer(
    r1:r3,
    names_to = "replicate",
    values_to = "conc"
  ) |>
  dplyr::filter(!is.na(conc)) |>
  mutate(
    amino_acid = factor(
      amino_acid,
      levels = c(
        "Glutamic acid",
        "Aspartic acid",
        "Leucine"
      )
    ),
    treatment = factor(
      treatment,
      levels = codes
    )
  )

# ---------------------------------------------------------------------
# 2. Summary statistics
# ---------------------------------------------------------------------

summ <- dat |>
  group_by(amino_acid, treatment) |>
  summarise(
    n = n(),
    mean = mean(conc),
    sd = sd(conc),
    ymax = mean + sd,
    .groups = "drop"
  )

print(summ, n = Inf)

# ---------------------------------------------------------------------
# 3. Levene Test
# ---------------------------------------------------------------------

levene_tbl <- dat |>
  group_by(amino_acid) |>
  summarise(
    levene_p = car::leveneTest(conc ~ treatment)$`Pr(>F)`[1],
    .groups = "drop"
  )

print(levene_tbl)

# ---------------------------------------------------------------------
# 4. ANOVA
# ---------------------------------------------------------------------

anova_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~{
    
    fit <- aov(conc ~ treatment, data = .x)
    
    broom::tidy(fit) |>
      dplyr::filter(term == "treatment") |>
      dplyr::select(df, statistic, p.value)
    
  }) |>
  ungroup()

print(anova_tbl)

# ---------------------------------------------------------------------
# 5. Tukey HSD
# ---------------------------------------------------------------------

tukey_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~{
    
    fit <- aov(conc ~ treatment, data = .x)
    
    as_tibble(
      TukeyHSD(fit)$treatment,
      rownames = "comparison"
    )
    
  }) |>
  ungroup()

print(tukey_tbl, n = Inf)

# ---------------------------------------------------------------------
# 6. Compact Letter Display
# ---------------------------------------------------------------------

letters_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~{
    
    fit <- aov(conc ~ treatment, data = .x)
    
    let <- multcompLetters4(
      fit,
      TukeyHSD(fit)
    )$treatment$Letters
    
    tibble(
      treatment = names(let),
      letter = unname(let)
    )
    
  }) |>
  ungroup() |>
  mutate(
    treatment = factor(
      treatment,
      levels = codes
    )
  )

plot_summ <- left_join(
  summ,
  letters_tbl,
  by = c("amino_acid", "treatment")
)

# ---------------------------------------------------------------------
# 7. Percent Change Relative to Water
# ---------------------------------------------------------------------

pct_tbl <- summ |>
  dplyr::select(amino_acid, treatment, mean) |>
  pivot_wider(
    names_from = treatment,
    values_from = mean
  ) |>
  mutate(
    across(
      c(ND, Bn, Br),
      ~ round(100 * (.x - W) / W, 1),
      .names = "pct_{.col}_vs_W"
    )
  )

print(pct_tbl)

# ---------------------------------------------------------------------
# 8. Plot
# ---------------------------------------------------------------------

dodge <- position_dodge(width = 0.8)

p <- ggplot(
  plot_summ,
  aes(
    x = amino_acid,
    y = mean,
    fill = treatment
  )
) +
  
  geom_col(
    position = dodge,
    width = 0.8,
    colour = "white"
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    position = dodge,
    width = 0.25,
    linewidth = 0.8
  ) +
  
  geom_text(
    data = dplyr::filter(
      plot_summ,
      amino_acid != "Aspartic acid"
    ),
    aes(
      y = ymax + 0.05,
      label = letter
    ),
    position = dodge,
    fontface = "bold",
    size = 6
  ) +
  
  annotate(
    "segment",
    x = 1.7,
    xend = 2.3,
    y = 1.35,
    yend = 1.35,
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = 2,
    y = 1.39,
    label = "ns",
    fontface = "bold",
    size = 6
  ) +
  
  scale_fill_manual(
    values = cols,
    labels = labels,
    name = "Sample Legend"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.6),
    breaks = seq(0, 1.6, 0.2),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Free Amino Acid",
    y = "Free Amino Acid Concentration (mg/g)"
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 120),
    
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 8,
      face = "bold",
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 1.2
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 1.2
    ),
    
    panel.grid.major.y = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank()
  )

print(p)

# ---------------------------------------------------------------------
# 9. Export (Portrait)
# ---------------------------------------------------------------------

ggsave(
  "Free_Amino_Acids_Portrait_Poster.png",
  p,
  width = 7,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Free_Amino_Acids_Portrait_Poster.pdf",
  p,
  width = 7,
  height = 9
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    "Free_Amino_Acids_Portrait_Poster.svg",
    p,
    width = 7,
    height = 9
  )
}

# ---------------------------------------------------------------------
# 10. Save statistics
# ---------------------------------------------------------------------

write_csv(anova_tbl, "anova_results.csv")
write_csv(tukey_tbl, "tukey_results.csv")
write_csv(plot_summ, "summary_means_sd_letters.csv")