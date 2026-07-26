# -----------------------------
# Apparent Protein Digestibility Data
# -----------------------------
df_protein <- data.frame(
  Sample = c("DW-Water", "DW-Banana", "DW-Bread"),
  Mean   = c(28.90, 42.19, 80.95),
  SD     = c(4.36, 8.77, 4.42),
  CV     = c(15.10, 20.78, 5.46),
  Tukey  = c("b", "c", "d")   # <-- YOUR REQUEST
)

df_protein$Sample <- factor(df_protein$Sample,
                            levels = c("DW-Water", "DW-Banana", "DW-Bread"))

library(ggplot2)

# -----------------------------
# Base Plot
# -----------------------------
p_protein <- ggplot(df_protein, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = Tukey, y = Mean + SD + 5), size = 7) +
  
  scale_fill_manual(
    values = c("DW-Water" = "steelblue",
               "DW-Banana" = "forestgreen",
               "DW-Bread" = "saddlebrown"),
    name = "Sample Legend"
  ) +
  
  geom_point(aes(shape = Tukey), alpha = 0) +
  scale_shape_manual(
    name = "Tukey Groups",
    values = c("b" = 16, "c" = 16, "d" = 16),
    labels = c("b = DW-Water",
               "c = DW-Banana",
               "d = DW-Bread")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    
    # ⭐ Title ABOVE panel (journal style)
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 1.2)
  ) +
  labs(
    title = "Apparent Protein Digestibility of Duckweed in Different Food Matrices",
    x = "Food Vehicle",
    y = "Apparent Protein Digestibility (%)"
  )

# -----------------------------
# Asterisks + comparison brackets
# -----------------------------
p_protein <- p_protein +
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 100, label = "**", size = 7) +
  annotate("segment", x = 2, xend = 3, y = 105, yend = 105, linewidth = 0.8) +
  annotate("text", x = 2.5, y = 110, label = "****", size = 7)

# -----------------------------
# Significance legend inside plot
# -----------------------------
p_protein <- p_protein +
  annotate("text", x = 1, y = 120,
           label = "** p < 0.01   |   **** p < 0.0001",
           hjust = 0, size = 5)

print(p_protein)
