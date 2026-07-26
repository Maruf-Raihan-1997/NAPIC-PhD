#This Code is made by Maruf Raihan for FOLIN data statistical analysis and visual presentation

# -----------------------------
# Raw data
# -----------------------------
undig <- c(22.1, 21.0, 23.1)
water <- c(80.0, 65.4, 64.1)
banana <- c(64.4, 68.4, 65.3)
banana_ctrl <- c(12.0, 13.3, 12.1)
dw_banana <- c(52.4, 55.1, 53.2)
dw_bread <- c(114.31, 99.36, 104.06)
bread_ctrl <- c(9.18, 9.53, 9.76)

value <- c(undig, water, banana, banana_ctrl, dw_banana, dw_bread, bread_ctrl)
group <- factor(rep(c("Undigested", "Water", "Banana", "BananaCtrl",
                      "DW_Banana", "DW_Bread", "BreadCtrl"), each = 3))

df <- data.frame(value, group)

# -----------------------------
# ANOVA + Tukey
# -----------------------------
anova_result <- aov(value ~ group, data = df)
summary(anova_result)
TukeyHSD(anova_result)

# -----------------------------
# ggplot2
# -----------------------------
library(ggplot2)

df_bar <- data.frame(
  Sample = c("Undigested", "DW-Water", "DW-Banana", "DW-Bread"),
  Mean   = c(22.00, 69.80, 53.50, 105.90),
  SD     = c(1.00, 8.40, 1.30, 7.65),
  Tukey  = c("a", "b", "c", "d")
)

df_bar$Sample <- factor(df_bar$Sample,
                        levels = c("Undigested", "DW-Water", "DW-Banana", "DW-Bread"))

p <- ggplot(df_bar, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = Tukey, y = Mean + SD + 10), size = 7) +
  
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
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),   # <-- TITLE FIXED
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                               size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 1.2)
  ) +
  labs(
    title = "In Vitro Digestion Increased the Total Phenolic Content of Duckweed",
    x = "Sample Type",
    y = "TPC (mg GAE/g dry weight)"
  )

# -----------------------------
# Asterisks + comparison brackets
# -----------------------------
p <- p +
  annotate("segment", x = 2, xend = 4, y = 130, yend = 130, linewidth = 0.8) +
  annotate("text", x = 3, y = 135, label = "****", size = 7) +
  annotate("segment", x = 3, xend = 4, y = 145, yend = 145, linewidth = 0.8) +
  annotate("text", x = 3.5, y = 150, label = "****", size = 7) +
  annotate("segment", x = 1, xend = 2, y = 120, yend = 120, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 125, label = "**", size = 7)

# -----------------------------
# Significance legend inside plot
# -----------------------------
p <- p +
  annotate("text", x = 1, y = 180,
           label = "**** p < 0.0001   |   ** p = 0.0031",
           hjust = 0, size = 5)

print(p)
