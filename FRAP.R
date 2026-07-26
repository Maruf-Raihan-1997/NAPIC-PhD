# -----------------------------
# FRAP Data
# -----------------------------
df_frap <- data.frame(
  Sample = c("Undigested", "DW-Water", "DW-Banana", "DW-Bread"),
  Mean   = c(30.20, 69.80, 45.30, 70.66),
  SD     = c(4.80, 1.00, 3.70, 9.30),
  CV     = c(15.80, 1.50, 8.20, 13.17),
  Tukey  = c("a", "b", "c", "d")   # DW-Bread = d
)

df_frap$Sample <- factor(df_frap$Sample,
                         levels = c("Undigested", "DW-Water", "DW-Banana", "DW-Bread"))

library(ggplot2)

# -----------------------------
# Base FRAP Plot
# -----------------------------
p_frap <- ggplot(df_frap, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = Tukey, y = Mean + SD + 8), size = 7) +
  
  scale_fill_manual(
    values = c("Undigested" = "darkgreen",
               "DW-Water" = "steelblue",
               "DW-Banana" = "forestgreen",
               "DW-Bread" = "saddlebrown"),
    name = "Sample Legend"
  ) +
  
  geom_point(aes(shape = Tukey), alpha = 0) +
  scale_shape_manual(
    name = "Tukey Groups",
    values = c("a" = 16, "b" = 16, "c" = 16, "d" = 16),
    labels = c("a = Undigested",
               "b = DW-Water",
               "c = DW-Banana",
               "d = DW-Bread")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    
    # ⭐ TITLE ABOVE PANEL (journal style)
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 1.2)
  ) +
  labs(
    title = "In Vitro Digestion Increased the Antioxidant Capacity (FRAP) of Duckweed",
    x = "Food Vehicle",
    y = "FRAP (µmol TE/g dry weight)"
  )

# -----------------------------
# Asterisks + comparison brackets
# -----------------------------
p_frap <- p_frap +
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 100, label = "**", size = 7) +
  annotate("segment", x = 2, xend = 4, y = 105, yend = 105, linewidth = 0.8) +
  annotate("text", x = 3, y = 110, label = "****", size = 7)

# -----------------------------
# Significance legend inside plot
# -----------------------------
p_frap <- p_frap +
  annotate("text", x = 1, y = 125,
           label = "**** p < 0.0001   |   ** p < 0.01",
           hjust = 0, size = 5)

print(p_frap)
