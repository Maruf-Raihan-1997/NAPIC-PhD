# ---------------------------------------------------------
# DL-SERINE STANDARD CURVE 
# ---------------------------------------------------------

library(ggplot2)

# 1. Standard curve data
std_df <- data.frame(
  Serine = c(0, 50, 100, 150, 200),
  Abs    = c(0.00, 0.09, 0.18, 0.26, 0.34)
)

pos_df <- data.frame(
  Serine = 186.19,
  Abs    = 0.32,
  Type   = "Positive control (BSA)"
)

neg_df <- data.frame(
  Serine = 0,
  Abs    = 0.00,
  Type   = "Negative control (DDH2O)"
)

# 2. Regression model (from excel master sheet)
eq <- "y = 0.0017x + 0.0063"
r2 <- "R² = 0.9989"

# 4. Inset plot with crisp fonts
p_std_small <- ggplot(std_df, aes(x = Serine, y = Abs)) +
  geom_point(size = 3.5, colour = "blue") +
  
  # Solid regression line
  geom_abline(
    intercept = 0.0063,
    slope     = 0.0017,
    colour    = "black",
    linewidth = 1.5
  ) +
  
  geom_point(data = pos_df, aes(colour = Type), size = 4.5) +
  geom_point(data = neg_df, aes(colour = Type), size = 4.5) +
  
  scale_colour_manual(
    values = c("Negative control (DDH2O)" = "red",
               "Positive control (BSA)" = "green"),
    labels = c(
      expression("Negative control (DDH"["2"]*"O)"),
      "Positive control (BSA)"
    )
  ) +
  
  labs(
    title = "DL-Serine Standard Curve",
    x = "Serine concentration (mg/L)",
    y = "Blank corrected Absorbance at 340nm",
    colour = NULL
  ) +
  
  # CRISP, LARGE EQUATION + R²
  annotate("text", x = 20, y = 0.34, label = eq, size = 7, fontface = "bold", hjust = 0) +
  annotate("text", x = 20, y = 0.30, label = r2, size = 7, fontface = "bold", hjust = 0) +
  
  theme_bw(base_size = 16) +
  theme(
    # Centered title
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    
    # Axis titles (Y-axis smaller)
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    # Axis text
    axis.text = element_text(size = 14),
    
    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.width  = unit(0.15, "cm"),
    legend.key.height = unit(0.15, "cm"),
    legend.spacing.x  = unit(0.15, "cm"),
    legend.spacing.y  = unit(0.15, "cm"),
    legend.margin     = margin(t = -4)
  )

# ---------------------------------------------------------
# 5. PRINT PLOT
# ---------------------------------------------------------

print(p_std_small)
