#This Code is made by Maruf Raihan (PhD,NICHE,Ulster)for Folin data statistical analysis and visual presentation


library(dplyr)
library(ggplot2)
library(multcompView)

# -----------------------------
# Raw data
# -----------------------------
TPC <- c(
  22.1, 21.0, 23.1,        # Undigested DW
  80.0, 65.4, 64.1,        # 10%(w/v)\nDW Water
  52.4, 55.1, 53.2,        # 10%(w/w)\nDW Banana
  105.13, 89.83, 94.30     # 10%(w/w)\nDW Bread
)

group <- factor(rep(c("Undigested DW", "10%(w/v)\nDW Water", "10%(w/w)\nDW Banana", "10%(w/w)\nDW Bread"), each = 3))
df <- data.frame(TPC, group)

# -----------------------------
# ANOVA + Tukey
# -----------------------------
anova_result <- aov(TPC ~ group, data = df)
tukey <- TukeyHSD(anova_result)

# Extract p-values
pvals <- tukey$group[, "p adj"]
names(pvals) <- rownames(tukey$group)

# -----------------------------
# Tukey LETTERS (NEW)
# -----------------------------
letters <- multcompLetters(pvals)$Letters
letters_df <- data.frame(
  Sample = names(letters),
  Letter = letters
)

# -----------------------------
# Summary table (means + SD)
# -----------------------------
df_summary <- df %>%
  group_by(group) %>%
  summarise(
    Mean = mean(TPC),
    SD   = sd(TPC),
    .groups = "drop"
  )

df_summary$Sample <- factor(df_summary$group,
                            levels = c("Undigested DW", "10%(w/v)\nDW Water", "10%(w/w)\nDW Banana", "10%(w/w)\nDW Bread"))

# Merge Tukey letters into summary
df_summary <- merge(df_summary, letters_df, by.x = "Sample", by.y = "Sample")

# -----------------------------
# Function: convert p-value → asterisks
# -----------------------------
p_to_star <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  return("ns")
}

# -----------------------------
# Build automatic bracket table
# -----------------------------
comparisons <- data.frame(
  comp = names(pvals),
  p = pvals,
  stars = sapply(pvals, p_to_star)
)

# Split "A-B" into two groups
comparisons <- comparisons %>%
  mutate(
    g1 = sub("-.*", "", comp),
    g2 = sub(".*-", "", comp)
  )

# -----------------------------
# REMOVE ONLY:
# Undigested DW – 10%(w/w)\nDW Banana
# Undigested DW – 10%(w/w) \nDW Bread
# 10%(w/v)\nDW Water – 10%(w/w) \nDW Bread
# -----------------------------
remove_comps <- c(
  "Undigested DW-10%(w/w)\nDW Banana",
  "10%(w/w)\nDW Banana-Undigested DW",
  "10%(w/v)\nDW Water-10%(w/w)\nDW Bread",
  "10%(w/w)\nDW Bread-10%(w/v)\nDW Water",
  "Undigested DW-10%(w/w)\nDW Bread"
)

comparisons <- comparisons %>%
  filter(!(comp %in% remove_comps))

# -----------------------------
# Assign bracket height automatically
# -----------------------------
max_y <- max(df_summary$Mean + df_summary$SD)
comparisons$y <- seq(max_y + 10, max_y + 10 + 15*(nrow(comparisons)-1), by = 15)

# -----------------------------
# ggplot
# -----------------------------
p <- ggplot(df_summary, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  
  # -----------------------------
# ADD TUKEY LETTERS ABOVE BARS
# -----------------------------
geom_text(
  aes(label = Letter, y = Mean + SD + 5),
  size = 7,
  fontface = "bold"
) +
  
  scale_fill_manual(
    values = c("Undigested DW" = "darkgreen",
               "10%(w/v)\nDW Water"   = "steelblue",
               "10%(w/w)\nDW Banana"  = "forestgreen",
               "10%(w/w)\nDW Bread"   = "saddlebrown"),
    labels = c(
      "Undigested DW",
      "10% (w/v) DW Water",
      "10% (w/w) DW Banana",
      "10% (w/w) DW Bread"
    ),
    name = "Sample Legend"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Total Phenolic Content (Folin Assay)",
    x = "Food Vehicle",
    y = "TPC (mg GAE/g dry weight)"
  )

# -----------------------------
# Add automatic brackets + asterisks
# -----------------------------
for (i in 1:nrow(comparisons)) {
  g1 <- comparisons$g1[i]
  g2 <- comparisons$g2[i]
  y  <- comparisons$y[i]
  stars <- comparisons$stars[i]
  
  x1 <- which(levels(df_summary$Sample) == g1)
  x2 <- which(levels(df_summary$Sample) == g2)
  
  p <- p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.8) +
    annotate("text", x = (x1 + x2)/2, y = y + 3, label = stars, size = 6)
}

print(p)

# -----------------------------
# Console output
# -----------------------------
cat("\n=== DESCRIPTIVE STATS + TUKEY LETTERS ===\n")
print(df_summary)

cat("\n=== ANOVA SUMMARY ===\n")
print(summary(anova_result))

cat("\n=== TUKEY POST-HOC ===\n")
print(tukey)

cat("\n=== SIGNIFICANCE STARS USED IN PLOT ===\n")
print(comparisons[, c("comp", "p", "stars")])
