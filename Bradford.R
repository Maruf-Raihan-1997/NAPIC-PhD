# -----------------------------
# Bradford Protein Concentration Data
# -----------------------------
df_bradford <- data.frame(
  Sample = c("DW-Water", "DW-Banana", "DW-Bread"),
  Mean   = c(604.19, 429.80, 459.57),
  SD     = c(118.63, 25.33, 18.17),
  CV     = c(19.64, 5.89, 3.95),
  Tukey  = c("b", "c", "d")   # <-- YOUR REQUEST
)

df_bradford$Sample <- factor(df_bradford$Sample,
                             levels = c("DW-Water", "DW-Banana", "DW-Bread"))

library(ggplot2)

# -----------------------------
# Base Bradford Plot
# -----------------------------
p_bradford <- ggplot(df_bradford, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = Tukey, y = Mean + SD + 30), size = 7) +
  
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
    title = "Protein Concentration of Duckweed After In Vitro Digestion in Different Food Matrices",
    x = "Food Vehicle",
    y = "Protein Concentration (µg/mL)"
  )

# -----------------------------
# Asterisks + comparison brackets
# -----------------------------
p_bradford <- p_bradford +
  annotate("segment", x = 1, xend = 2, y = 800, yend = 800, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 830, label = "**", size = 7) +
  annotate("segment", x = 2, xend = 3, y = 880, yend = 880, linewidth = 0.8) +
  annotate("text", x = 2.5, y = 910, label = "****", size = 7)

# -----------------------------
# Significance legend inside plot
# -----------------------------
p_bradford <- p_bradford +
  annotate("text", x = 1, y = 950,
           label = "** p < 0.01   |   **** p < 0.0001",
           hjust = 0, size = 5)

print(p_bradford)
