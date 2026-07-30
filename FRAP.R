#This Code is made by Maruf Raihan (PhD,NICHE,Ulster)for FRAP data statistical analysis and visual presentation

library(dplyr)
library(ggplot2)
library(multcompView)

# -----------------------------
# Raw FRAP data
# -----------------------------
FRAP <- c(
  33.80, 32.20, 24.80,      # Undigested
  68.80, 68.10, 70.20,      # DW-Water
  41.60, 45.30, 49.00,      # DW-Banana
  81.40, 65.61, 64.97       # DW-Bread
)

group <- factor(rep(c("Undigested", "DW_Water", "DW_Banana", "DW_Bread"), each = 3))
df <- data.frame(FRAP, group)

# -----------------------------
# ANOVA + Tukey
# -----------------------------
anova_result <- aov(FRAP ~ group, data = df)
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
    Mean = mean(FRAP),
    SD   = sd(FRAP),
    .groups = "drop"
  )

df_summary$Sample <- factor(df_summary$group,
                            levels = c("Undigested", "DW_Water", "DW_Banana", "DW_Bread"))

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
# Undigested – DW_Banana
# Undigested – DW_Bread
# DW_Water – DW_Bread
# -----------------------------
remove_comps <- c(
  "Undigested-DW_Banana",
  "DW_Banana-Undigested",
  "DW_Water-DW_Bread",
  "DW_Bread-DW_Water",
  "Undigested-DW_Bread"
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
    values = c("Undigested" = "darkgreen",
               "DW_Water"   = "steelblue",
               "DW_Banana"  = "forestgreen",
               "DW_Bread"   = "saddlebrown"),
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
    title = "Antioxidant Capacity (FRAP assay)",
    x = "Food Vehicle",
    y = "FRAP (µmol TE/g dry weight)"
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
