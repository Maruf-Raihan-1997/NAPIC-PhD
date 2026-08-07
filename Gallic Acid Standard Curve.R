## Gallic Acid Standard Curve – Folin–Ciocalteu Assay
## ---------------------------------------------------------
library(ggplot2)

# 1. Standard curve data
std_df <- data.frame(
  GA  = c(10, 20, 30, 40, 50),
  Abs = c(0.095, 0.182, 0.259, 0.361, 0.465)
)

# Positive control: Quercetin
pos_df <- data.frame(
  GA   = 70.16,
  Abs  = 0.6104,
  Type = "Positive control (Quercetin)"
)

# Negative control: 95% MeOH
neg_df <- data.frame(
  GA   = 0,
  Abs  = 0.00,
  Type = "Negative control (95% MeOH)"
)

# 2. Linear regression model
model <- lm(Abs ~ GA, data = std_df)

intercept <- coef(model)[1]
slope     <- coef(model)[2]

# Forced minus sign (your requirement)
eq <- sprintf("y = %.4fx - %.4f", slope, abs(intercept))
r2 <- sprintf("R² = %.4f", summary(model)$r.squared)

# 3. Plot (matching OPA style exactly)
p_ga <- ggplot(std_df, aes(x = GA, y = Abs)) +
  geom_point(size = 3.5, colour = "blue") +
  
  geom_abline(
    intercept = intercept,
    slope     = slope,
    colour    = "black",
    linewidth = 1.5
  ) +
  
  geom_point(data = pos_df, aes(colour = Type), size = 4.5) +
  geom_point(data = neg_df, aes(colour = Type), size = 4.5) +
  
  # ⭐ FIX: Prevent auto-rescaling so fonts match OPA
  scale_y_continuous(limits = c(0, 0.65)) +
  
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
    title = "Gallic Acid Standard Curve",
    x = "Gallic acid concentration (mg/L)",
    y = "Blank corrected Absorbance at 760 nm",
    colour = NULL
  ) +
  
  # ⭐ MATCH OPA EXACTLY: annotation relative to top of standards
  annotate("text", x = 8,  y = 0.62, label = eq,
           size = 7, fontface = "bold", hjust = 0) +
  annotate("text", x = 8,  y = 0.56, label = r2,
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

print(p_ga)
