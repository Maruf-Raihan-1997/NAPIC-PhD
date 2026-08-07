## Trolox Standard Curve – FRAP Assay
## ---------------------------------------------------------
library(ggplot2)

# 1. Standard curve data
std_df <- data.frame(
  Trolox = c(100, 200, 400, 600, 800),
  Abs    = c(0.229, 0.507, 0.706, 0.909, 1.215)
)

# Positive control: Quercetin (40 mg/mL)
pos_df <- data.frame(
  Trolox = 1670.769231,   # TEAC measured (µM)
  Abs    = 2.339,
  Type   = "Positive control (Quercetin)"
)

# Negative control: 95% Methanol
neg_df <- data.frame(
  Trolox = 0,
  Abs    = 0.00,
  Type   = "Negative control (95% MeOH)"
)

# 2. Linear regression model
model <- lm(Abs ~ Trolox, data = std_df)

intercept <- coef(model)[1]
slope     <- coef(model)[2]

# regression
eq <- sprintf("y = %.4fx + %.4f", slope, intercept)
r2 <- sprintf("R² = %.4f", summary(model)$r.squared)

# 3. Plot (matching OPA style exactly)
p_trolox <- ggplot(std_df, aes(x = Trolox, y = Abs)) +
  geom_point(size = 3.5, colour = "blue") +
  
  geom_abline(
    intercept = intercept,
    slope     = slope,
    colour    = "black",
    linewidth = 1.5
  ) +
  
  geom_point(data = pos_df, aes(colour = Type), size = 4.5) +
  geom_point(data = neg_df, aes(colour = Type), size = 4.5) +
  
  # ⭐ Prevent auto-rescaling (positive control is above standards)
  scale_y_continuous(limits = c(0, 2.5)) +
  
  scale_colour_manual(
    values = c(
      "Negative control (95% MeOH)" = "red",
      "Positive control (Quercetin)" = "green"
    ),
    labels = c(
      "Negative control (95% MeOH)",
      "Positive control (Quercetin)"
    )
  ) +
  
  labs(
    title = "Trolox Standard Curve",
    x = "Trolox concentration (µM)",
    y = "Blank corrected Absorbance at 593 nm",
    colour = NULL
  ) +
  
  # ⭐ Annotation positioned relative to top of standards (same as OPA)
  annotate("text", x = 120, y = 2.44, label = eq,
           size = 7, fontface = "bold", hjust = 0) +
  annotate("text", x = 120, y = 2.20, label = r2,
           size = 7, fontface = "bold", hjust = 0) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 14),
    
    axis.line  = element_line(colour = "black", linewidth = 1.5),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    legend.key      = element_rect(fill = NA, colour = NA),
    legend.key.width  = unit(0.15, "cm"),
    legend.key.height = unit(0.15, "cm"),
    legend.spacing.x  = unit(0.15, "cm"),
    legend.spacing.y  = unit(0.15, "cm"),
    legend.margin     = margin(t = -4)
  )

print(p_trolox)
