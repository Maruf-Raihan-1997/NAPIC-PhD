# =====================================================================
# Free amino acid release from Lemna minor duckweed after in vitro digestion
# Statistics (TWO-way ANOVA + INTERACTION)
#A two-factor ANOVA always tests:
  
 # Main effect of Factor A
#Main effect of Factor B
#A × B interaction
# Conference Poster Version (Portrait)
# =====================================================================

library(tidyverse)

# ---------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------

codes <- c("ND", "W", "Bn", "Br")

labels <- c(
  ND = "Undigested DWP",
  W = "10% (w/v) DWP Water",
  Bn = "10% (w/w) DWP  Banana",
  Br = "10% (w/w) DWP  Bread"
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

# =====================================================================
# TWO-WAY ANOVA
# =====================================================================

library(tidyverse)
library(car)
library(broom)

# ---------------------------------------------------------------------
# Two-way ANOVA
# ---------------------------------------------------------------------

fit2 <- aov(conc ~ treatment * amino_acid, data = dat)

cat("\n=== TWO-WAY ANOVA ===\n")
print(summary(fit2))

# =====================================================================
# ASSUMPTION CHECKS
# =====================================================================

# ---------------------------------------------------------------------
# Normality of residuals
# ---------------------------------------------------------------------

res <- residuals(fit2)

cat("\n=== SHAPIRO-WILK TEST ===\n")
print(shapiro.test(res))

qqnorm(res)
qqline(res, col = "red", lwd = 2)

hist(
  res,
  main = "Residual Distribution",
  xlab = "Residuals"
)

# ---------------------------------------------------------------------
# Homogeneity of Variance
# ---------------------------------------------------------------------

cat("\n=== LEVENE TEST ===\n")

levene_result <- leveneTest(
  conc ~ interaction(treatment, amino_acid),
  data = dat
)

print(levene_result)

# =====================================================================
# SUMMARY STATISTICS
# =====================================================================

plot_data <- dat %>%
  group_by(amino_acid, treatment) %>%
  summarise(
    n = n(),
    mean = mean(conc),
    sd = sd(conc),
    .groups = "drop"
  )

cat("\n=== SUMMARY STATISTICS ===\n")
print(plot_data)

# =====================================================================
# BAR CHART (MEAN ± SD)
# =====================================================================

dodge <- position_dodge(width = 0.8)

p_twoway <- ggplot(
  plot_data,
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
  
  scale_fill_manual(
    values = cols,
    labels = labels,
    name = "Sample Legend"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
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
      size = 12,
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

print(p_twoway)

dev.new()
print(p_interaction)

# ---------------------------------------------------------------------
# Interaction Plot with 95% CI

#Think of 95% confidence interval (CI) as:

#"Given my sample data, the true population mean is likely to lie somewhere within this range."

#It tells you how precise your estimate is, whereas the mean tells you the best estimate itself.
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Interaction Plot with Mean ± SD
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Interaction Plot (Mean ± SD)
# ---------------------------------------------------------------------

interaction_summary <- dat %>%
  group_by(amino_acid, treatment) %>%
  summarise(
    n = n(),
    mean = mean(conc),
    sd = sd(conc),
    .groups = "drop"
  )

# Wrapped labels for poster readability
interaction_labels <- c(
  ND = "Undigested\nDWP",
  W = "10% (w/v)\nDWP Water",
  Bn = "10% (w/w)\nDWP Banana",
  Br = "10% (w/w)\nDWP Bread"
)

p_interaction <- ggplot(
  interaction_summary,
  aes(
    x = treatment,
    y = mean,
    colour = amino_acid,
    group = amino_acid
  )
) +
  
  geom_line(linewidth = 1.8) +
  
  geom_point(size = 5) +
  
  geom_errorbar(
    aes(
      ymin = pmax(mean - sd, 0),
      ymax = mean + sd
    ),
    width = 0.08,
    linewidth = 1
  ) +
  
  scale_x_discrete(
    labels = interaction_labels
  ) +
  
  scale_colour_manual(
    values = c(
      "Glutamic acid" = "#D55E00",
      "Aspartic acid" = "#0072B2",
      "Leucine" = "#009E73"
    )
  ) +
  
  labs(
    title = "Food Matrix × Amino Acid Interaction",
    x = "Treatment (Food Matrix )",
    y = "Free Amino Acid Concentration (mg/g)",
    colour = "Amino Acid"
  ) +
  
  theme_classic(base_size = 18) +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 16
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    axis.text.x = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      lineheight = 0.9
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    panel.grid.minor = element_blank()
  )

print(p_interaction)


# ---------------------------------------------------------------------
# Interaction Plot (Mean ± SD)
# Food Matrix × Amino Acid Interaction
# Thick x- and y-axes, no top/right border
# ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(grid)

# Summary statistics
interaction_summary <- dat %>%
  group_by(amino_acid, treatment) %>%
  summarise(
    n = n(),
    mean = mean(conc, na.rm = TRUE),
    sd = sd(conc, na.rm = TRUE),
    .groups = "drop"
  )

# Wrapped labels for poster readability
interaction_labels <- c(
  ND = "Undigested\nDWP",
  W  = "10% (w/v)\nDWP Water",
  Bn = "10% (w/w)\nDWP Banana",
  Br = "10% (w/w)\nDWP Bread"
)

# Plot
p_interaction <- ggplot(
  interaction_summary,
  aes(
    x = treatment,
    y = mean,
    colour = amino_acid,
    group = amino_acid
  )
) +
  
  geom_line(linewidth = 2.0) +
  
  geom_point(size = 5) +
  
  geom_errorbar(
    aes(
      ymin = pmax(mean - sd, 0),
      ymax = mean + sd
    ),
    width = 0.08,
    linewidth = 1.2
  ) +
  
  scale_x_discrete(
    labels = interaction_labels
  ) +
  
  scale_colour_manual(
    values = c(
      "Glutamic acid" = "#D55E00",
      "Aspartic acid" = "#0072B2",
      "Leucine"       = "#009E73"
    )
  ) +
  
  labs(
    title = "Food Matrix × Amino Acid Interaction",
    x = "Treatment (Food Matrix)",
    y = "Free Amino Acid Concentration (mg/g)",
    colour = "Amino Acid"
  ) +
  
  theme_classic(base_size = 18) +
  
  theme(
    
    # Thick x and y axes only
    axis.line = element_line(
      colour = "black",
      linewidth = 2
    ),
    
    # Thick ticks
    axis.ticks = element_line(
      colour = "black",
      linewidth = 2
    ),
    
    axis.ticks.length = unit(0.25, "cm"),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 16
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    axis.text.x = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      lineheight = 0.9
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_interaction)


# =====================================================================
# POST-HOC TESTS (Tukey)
# =====================================================================

library(emmeans)
library(multcomp)

# Estimated marginal means
emm <- emmeans(
  fit2,
  ~ treatment * amino_acid
)

# Compact letter display
letters_df <- cld(
  emm,
  by = "amino_acid",
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame()

# Remove spaces from letters
letters_df$.group <- gsub(" ", "", letters_df$.group)

print(letters_df)

# =====================================================================
# SUMMARY STATISTICS FOR PLOT
# =====================================================================

plot2 <- dat %>%
  group_by(amino_acid, treatment) %>%
  summarise(
    mean = mean(conc),
    sd = sd(conc),
    .groups = "drop"
  ) %>%
  left_join(
    letters_df %>%
      dplyr::select(
        amino_acid,
        treatment,
        .group
      ),
    by = c(
      "amino_acid",
      "treatment"
    )
  )

# =====================================================================
# BAR PLOT WITH TUKEY LETTERS
# =====================================================================

library(ggplot2)
library(grid)

dodge <- position_dodge(width = 0.8)

p_twoway <- ggplot(
  plot2,
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
      ymin = pmax(mean - sd, 0),
      ymax = mean + sd
    ),
    position = dodge,
    width = 0.25,
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = .group,
      y = mean + sd + 0.08
    ),
    position = dodge,
    size = 5,
    fontface = "bold",
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = cols,
    labels = labels,
    name = "Sample Legend"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  labs(
    title = "Free Amino Acid Release (mg/g)",
    x = "Amino Acid",
    y = "Free Amino Acid Concentration (mg/g)"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    # Remove legend
    legend.position = "none",
    
    # Center title
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5
    ),
    
    # Axis text
    axis.text.x = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    # Thick bottom and left axes only
    axis.line = element_line(
      colour = "black",
      linewidth = 2
    ),
    
    # Thick tick marks
    axis.ticks = element_line(
      colour = "black",
      linewidth = 2
    ),
    
    axis.ticks.length = unit(0.25, "cm"),
    
    # No top/right border, no grid
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_twoway)

# legend position



legend_plot <- ggplot(
  plot2,
  aes(x = amino_acid, y = mean, fill = treatment)
) +
  geom_col() +
  scale_fill_manual(
    values = cols,
    labels = labels,
    name = "Sample Legend"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_text(
      size = 20,
      face = "bold"
    ),
    legend.text = element_text(
      size = 18
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      keywidth = unit(1.5, "cm"),
      keyheight = unit(1.0, "cm")
    )
  )

legend <- cowplot::get_legend(legend_plot)

cowplot::ggdraw(legend)