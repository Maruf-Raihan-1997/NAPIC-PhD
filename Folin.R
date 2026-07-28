library(dplyr)
library(ggplot2)
library(multcompView)

# -----------------------------
# Raw data
# -----------------------------
TPC <- c(
  22.1, 21.0, 23.1,        # Undigested
  80.0, 65.4, 64.1,        # DW-Water
  52.4, 55.1, 53.2,        # DW-Banana
  105.13, 89.83, 94.30     # DW-Bread
)

group <- factor(rep(c("Undigested", "DW_Water", "DW_Banana", "DW_Bread"), each = 3))
df <- data.frame(TPC, group)

# -----------------------------
# ANOVA + Tukey
# -----------------------------
anova_result <- aov(TPC ~ group, data = df)
tukey <- TukeyHSD(anova_result)

# Extract p-values
pvals <- tukey$group[, "p adj"]
names(pvals) <- rownames(tukey$group)

# Tukey letters
letters <- multcompLetters(pvals)$Letters

# -----------------------------
# Summary table (means + SD + letters)
# -----------------------------
df_summary <- df %>%
  group_by(group) %>%
  summarise(
    Mean = mean(TPC),
    SD   = sd(TPC),
    .groups = "drop"
  ) %>%
  mutate(letter = letters[group])

df_summary$Sample <- factor(df_summary$group,
                            levels = c("Undigested", "DW_Water", "DW_Banana", "DW_Bread"))

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
# REMOVE ONLY the two pairs you specified
# -----------------------------
comparisons <- comparisons %>%
  filter(!(g1 == "Undigested" & g2 == "DW_Banana")) %>%   # remove Undigested–DW_Banana
  filter(!(g1 == "DW_Water" & g2 == "DW_Bread")) %>%      # remove DW_Water–DW_Bread
  filter(!(g1 == "DW_Bread" & g2 == "DW_Water"))          # remove reversed order too

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
  geom_text(aes(label = letter, y = Mean + SD + 10), size = 7) +
  
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
    title = "In Vitro Digestion Increased the Total Phenolic Content of Duckweed",
    x = "Food Vehicle",
    y = "TPC (mg GAE/g dry weight)"
  )

# -----------------------------
# Add automatic brackets + asterisks (except removed pairs)
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
