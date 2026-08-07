#This Code is made by Maruf Raihan (PhD,NICHE,Ulster)for OPA data statistical analysis and visual presentation

library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)

# ---------------------------------------------------------
# 1. RAW OPA DATA
# ---------------------------------------------------------
df_raw <- data.frame(
  Sample = factor(rep(c("10%(w/v)\nDW Water", "10%(w/w)\nDW Banana", "10%(w/w) \nDW Bread"), each = 3),
                  levels = c("10%(w/v)\nDW Water", "10%(w/w)\nDW Banana", "10%(w/w) \nDW Bread")),
  Value = c(
    28.36, 24.83, 33.51,   # 10%(w/v)\nDW Water
    32.78, 43.67, 50.13,   # 10%(w/w)\nDW Banana
    76.13, 84.81, 81.92    # 10%(w/w) \nDW Bread
  )
)

# ---------------------------------------------------------
# 2. DESCRIPTIVE STATS (auto)
# ---------------------------------------------------------
df_stats <- aggregate(Value ~ Sample, df_raw, function(x) {
  c(Mean = mean(x), SD = sd(x), CV = sd(x)/mean(x)*100)
})

df_stats <- do.call(data.frame, df_stats)
names(df_stats) <- c("Sample", "Mean", "SD", "CV")

# ---------------------------------------------------------
# 3. ANOVA + TUKEY
# ---------------------------------------------------------
model <- aov(Value ~ Sample, data = df_raw)
emm <- emmeans(model, ~ Sample)
tukey <- summary(pairs(emm, adjust = "tukey"))

# ---------------------------------------------------------
# 4. Tukey LETTERS
# ---------------------------------------------------------
cld_out <- suppressMessages(
  cld(emm, Letters = letters, adjust = "tukey")
)

letters_df <- as.data.frame(cld_out)

# Merge letters into stats
df_stats2 <- merge(df_stats, letters_df[, c("Sample", ".group")], by = "Sample")

# Extract p-values
p_WB  <- tukey$p.value[1]   # Water vs Banana
p_WBr <- tukey$p.value[2]   # Water vs Bread
p_BBr <- tukey$p.value[3]   # Banana vs Bread

# ---------------------------------------------------------
# 5. Convert p-values → asterisks
# ---------------------------------------------------------
p_to_star <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.01)    return("**")
  if (p < 0.05)    return("*")
  return("ns")
}

star_WB  <- p_to_star(p_WB)
star_WBr <- p_to_star(p_WBr)
star_BBr <- p_to_star(p_BBr)

# ---------------------------------------------------------
# 6. BASE PLOT
# ---------------------------------------------------------
p_opa <- ggplot(df_stats2, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  
  # Tukey letters ABOVE bars
  geom_text(
    aes(x = Sample, y = Mean + SD + 5, label = .group),
    size = 7,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c("10%(w/v)\nDW Water" = "steelblue",
               "10%(w/w)\nDW Banana" = "forestgreen",
               "10%(w/w) \nDW Bread" = "saddlebrown"),
    labels = c(
      "10% (w/v) DW Water",
      "10% (w/w) DW Banana",
      "10% (w/w) DW Bread"
    ),
    
    name = "Sample Legend"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 170),      # push legend downward
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 1.2)
  ) +
  labs(
    title = "Apparent Protein Digestibility % (OPA assay)",
    x = "Food Vehicle",
    y = "Degree of protein hydrolysis (%)"
  )

# ---------------------------------------------------------
# 7. AUTOMATED BRACKETS BASED ON ANOVA OUTPUT
# ---------------------------------------------------------

# Water vs Banana
p_opa <- p_opa +
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 100, label = star_WB, size = 7)

# Banana vs Bread
p_opa <- p_opa +
  annotate("segment", x = 2, xend = 3, y = 105, yend = 105, linewidth = 0.8) +
  annotate("text", x = 2.5, y = 110, label = star_BBr, size = 7)

# Water vs Bread
p_opa <- p_opa +
  annotate("segment", x = 1, xend = 3, y = 120, yend = 120, linewidth = 0.8) +
  annotate("text", x = 2, y = 125, label = star_WBr, size = 7)

print(p_opa)

cat("\n=== DESCRIPTIVE STATS ===\n")
print(df_stats2)

cat("\n=== ANOVA SUMMARY ===\n")
print(summary(model))

cat("\n=== TUKEY POST-HOC ===\n")
print(tukey)
