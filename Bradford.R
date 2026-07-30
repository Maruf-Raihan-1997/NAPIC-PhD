#This Code is made by Maruf Raihan (PhD,NICHE,Ulster)for Bradford data statistical analysis and visual presentation

library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)

# ---------------------------------------------------------
# 1. RAW DATA
# ---------------------------------------------------------
df <- data.frame(
  Sample = factor(rep(c("10%(w/v) DW-Water", "10%(w/w) DW-Banana", "10%(w/w) DW-Bread"), each = 3),
                  levels = c("10%(w/v) DW-Water", "10%(w/w) DW-Banana", "10%(w/w) DW-Bread")),
  Protein = c(
    639.00, 472.04, 701.52,   # 10%(w/v) DW-Water
    400.72, 441.69, 447.00,   # 10%(w/w) DW-Banana
    448.05, 480.52, 450.14    # 10%(w/w) DW-Bread"
  )
)

# ---------------------------------------------------------
# 2. DESCRIPTIVE STATS
# ---------------------------------------------------------
df_stats <- aggregate(Protein ~ Sample, df, function(x) {
  c(Mean = mean(x), SD = sd(x), CV = sd(x)/mean(x)*100)
})

df_stats <- do.call(data.frame, df_stats)
names(df_stats) <- c("Sample", "Mean", "SD", "CV")

# ---------------------------------------------------------
# 3. ANOVA + TUKEY
# ---------------------------------------------------------
model <- aov(Protein ~ Sample, data = df)
emm <- emmeans(model, ~ Sample)
tukey <- summary(pairs(emm, adjust = "tukey"))

# ---------------------------------------------------------
# 4. Tukey LETTERS
# ---------------------------------------------------------
cld_out <- suppressMessages(
  cld(emm, Letters = letters, adjust = "tukey")
)

df_plot <- merge(df_stats, cld_out[, c("Sample", ".group")], by = "Sample")

# ---------------------------------------------------------
# 5. BASE PLOT
# ---------------------------------------------------------
p_bradford <- ggplot(df_plot, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = .group, y = Mean + SD + 30), size = 7) +
  
  scale_fill_manual(
    name= "Sample Legend",
    values = c("10%(w/v) DW-Water" = "steelblue",
               "10%(w/w) DW-Banana" = "forestgreen",
               "10%(w/w) DW-Bread" = "saddlebrown")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 1.2)
  ) +
  labs(
    title = "Protein Concentration (Bradford Assay)",
    x = "Food Vehicle",
    y = "Protein Concentration (µg/mL)"
  ) 

# ---------------------------------------------------------
# 6. SINGLE ANOVA BRACKET ACROSS ALL THREE BARS
# ---------------------------------------------------------
#p_bradford <- p_bradford +
#  annotate("segment", x = 1, xend = 3, y = 1050, yend = 1050, linewidth = 1.2) +
# annotate("text", x = 2, y = 1080, label = "*", size = 8)

# ---------------------------------------------------------
# 7. PRINT
# ---------------------------------------------------------
print(p_bradford)

cat("\n=== DESCRIPTIVE STATS ===\n")
print(df_stats)

cat("\n=== ANOVA SUMMARY ===\n")
print(summary(model))

cat("\n=== TUKEY POST-HOC ===\n")
print(tukey)
