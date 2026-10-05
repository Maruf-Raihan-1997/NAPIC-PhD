# =====================================================================
# Free amino acid release from Lemna minor duckweed after in vitro digestion
# Statistics (one-way ANOVA + Tukey HSD per amino acid) and figure
# (mean +/- SD, compact letters)
#
# Packages: install.packages(c("tidyverse", "broom", "car", "multcompView"))
# (car is only used via car::leveneTest)
# Optional for SVG output: install.packages("svglite")
# Requires R >= 4.1 (native pipe |>)
# =====================================================================

library(tidyverse)
library(broom)
# car is called as car::leveneTest() (not attached) because attaching it
# loads MASS, which masks dplyr::select().
library(multcompView)  # multcompLetters4

# ---------------------------------------------------------------------
# 1. Data (mg/g, n = 3 replicates per treatment)
#    Short codes are used internally because "Non-disruptive" contains a
#    hyphen, which breaks Tukey comparison names ("A-B").
# ---------------------------------------------------------------------
codes  <- c("ND", "W", "Bn", "Br")
labels <- c(ND = "Undigested duckweed powder",
            W  = "10% (w/v) DW Water",
            Bn = "10% (w/w) DW Banana",
            Br = "10% (w/w) DW Bread")
# Colours taken from the OPA script: Water = steelblue, Banana = forestgreen,
# Bread = saddlebrown. ND (not in the OPA script) = grey40; change as you like.
cols   <- c(ND = "grey40", W = "steelblue", Bn = "forestgreen", Br = "saddlebrown")

dat <- tribble(
  ~amino_acid,     ~treatment, ~r1,  ~r2,  ~r3,
  "Glutamic acid", "ND",       0.44, 0.45, 0.48,
  "Glutamic acid", "W",        0.23, 0.21, 0.36,
  "Glutamic acid", "Bn",       0.05, 0.05, 0.08,
  "Glutamic acid", "Br",       0.04, 0.04, 0.05,
  "Aspartic acid", "ND",       0.96, 0.75, 0.79,
  "Aspartic acid", "W",        0.85, 1.01, 0.92,
  "Aspartic acid", "Bn",       0.86, 1.03, 0.79,
  "Aspartic acid", "Br",       0.99, 0.67, 1.27,
  "Leucine",       "ND",       0.22, 0.16, 0.21,
  "Leucine",       "W",        0.45, 0.45, 0.44,
  "Leucine",       "Bn",       0.16, 0.20, 0.14,
  "Leucine",       "Br",       0.17, 0.16, 0.17
) |>
  pivot_longer(r1:r3, names_to = "replicate", values_to = "conc") |>
  mutate(
    amino_acid = factor(amino_acid,
                        levels = c("Glutamic acid", "Aspartic acid", "Leucine")),
    treatment  = factor(treatment, levels = codes)
  )

# ---------------------------------------------------------------------
# 2. Descriptive statistics
# ---------------------------------------------------------------------
summ <- dat |>
  group_by(amino_acid, treatment) |>
  summarise(n = n(),
            mean = mean(conc),
            sd = sd(conc),
            ymax = mean + sd,                    # top of error bar
            .groups = "drop")

print(summ, n = Inf)

# ---------------------------------------------------------------------
# 3. Assumption check: Levene's test (median-centred) per amino acid
# ---------------------------------------------------------------------
levene_tbl <- dat |>
  group_by(amino_acid) |>
  summarise(levene_p = car::leveneTest(conc ~ treatment)$`Pr(>F)`[1],
            .groups = "drop")
print(levene_tbl)

# ---------------------------------------------------------------------
# 4. One-way ANOVA per amino acid
# ---------------------------------------------------------------------
anova_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~ tidy(aov(conc ~ treatment, data = .x)) |>
                 dplyr::filter(term == "treatment") |>
                 dplyr::select(df, statistic, p.value)) |>
  ungroup()
print(anova_tbl)

# ---------------------------------------------------------------------
# 5. Tukey HSD pairwise comparisons per amino acid
# ---------------------------------------------------------------------
tukey_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~ {
    fit <- aov(conc ~ treatment, data = .x)
    as_tibble(TukeyHSD(fit)$treatment, rownames = "comparison")
  }) |>
  ungroup()
print(tukey_tbl, n = Inf)

# ---------------------------------------------------------------------
# 6. Compact letter display (a = highest mean within each amino acid)
# ---------------------------------------------------------------------
letters_tbl <- dat |>
  group_by(amino_acid) |>
  group_modify(~ {
    fit <- aov(conc ~ treatment, data = .x)
    let <- multcompLetters4(fit, TukeyHSD(fit))$treatment$Letters
    tibble(treatment = names(let), letter = unname(let))
  }) |>
  ungroup() |>
  mutate(treatment = factor(treatment, levels = codes))

plot_summ <- left_join(summ, letters_tbl, by = c("amino_acid", "treatment"))

# ---------------------------------------------------------------------
# 7. Percent change relative to water-digested duckweed
# ---------------------------------------------------------------------
pct_tbl <- summ |>
  dplyr::select(amino_acid, treatment, mean) |>
  pivot_wider(names_from = treatment, values_from = mean) |>
  mutate(across(c(ND, Bn, Br), ~ round(100 * (.x - W) / W, 1),
                .names = "pct_{.col}_vs_W"))
print(pct_tbl)

# ---------------------------------------------------------------------
# 8. Figure
# ---------------------------------------------------------------------
dodge <- position_dodge(width = 0.8)

p <- ggplot(plot_summ, aes(x = amino_acid, y = mean, fill = treatment)) +
  geom_col(position = dodge, width = 0.8, colour = "white") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = dodge, width = 0.25, linewidth = 0.6) +
  # significance letters (aspartic acid is n.s., handled by the bracket below)
  geom_text(data = dplyr::filter(plot_summ, amino_acid != "Aspartic acid"),
            aes(y = ymax + 0.03, label = letter),
            position = dodge, vjust = 0, fontface = "bold", size = 5) +
  annotate("segment", x = 1.7, xend = 2.3, y = 1.35, yend = 1.35,
           linewidth = 0.5) +
  annotate("text", x = 2, y = 1.38, label = "ns",
           fontface = "bold", size = 5, vjust = 0) +
  scale_fill_manual(values = cols, labels = labels, name = NULL) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.4, 0.2),
                     expand = c(0, 0)) +
  labs(x = "Free Amino Acid",
       y = "Concentration (mg/g)",
       title = "Effect of food matrix on free amino acid release from\nLemna minor duckweed following in vitro digestion",
       caption = paste("Mean \u00b1 SD, n = 3 replicates.",
                       "Within each amino acid, bars sharing a letter do not differ",
                       "(one-way ANOVA with Tukey HSD, p < 0.05); ns = no significant differences.",
                       sep = "\n")) +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.caption = element_text(face = "italic", colour = "grey30",
                                    size = 10, hjust = 0),
        panel.grid.major.y = element_line(colour = "grey90"),
        axis.text = element_text(colour = "black"))

print(p)


