#This Code is made by Maruf Raihan (PhD,NICHE,Ulster)for ## Bradford assay standard curve – quadratic regression statistical analysis and visual presentation
library(ggplot2)

# 1. Standard curve data
std_df <- data.frame(
  BSA = c(200, 400, 600, 800),
  Abs = c(0.0657, 0.0947, 0.1677, 0.2970)
)

# Positive control at actual protein concentration
pos_df <- data.frame(
  BSA = 495.30,
  Abs = 0.12,
  Type = "Positive control (10% (w/v) BSA)"
)

neg_df <- data.frame(
  BSA = 0,
  Abs = 0.00,
  Type = "Negative control (DDH2O)"
)

# 2. Quadratic regression model
model <- lm(Abs ~ poly(BSA, 2, raw = TRUE), data = std_df)

c <- coef(model)[1]   # intercept
b <- coef(model)[2]   # linear term
a <- coef(model)[3]   # quadratic term

# CLEAN equation with correct + / – handling
eq <- sprintf(
  "y = %.3e x^2 %s %.3e x %s %.4f",
  a,
  ifelse(b < 0, " -", " +"),
  abs(b),
  ifelse(c < 0, " -", " +"),
  abs(c)
)

r2 <- sprintf("R² = %.4f", summary(model)$r.squared)

# 3. Plot (matching DL‑Serine style exactly)
p_bradford_quad <- ggplot(std_df, aes(x = BSA, y = Abs)) +
  geom_point(size = 3.5, colour = "blue") +
  
  # Quadratic regression curve
  stat_function(
    fun = function(x) a*x^2 + b*x + c,
    colour = "black",
    linewidth = 1.5
  ) +
  
  # Controls
  geom_point(data = pos_df, aes(colour = Type), size = 4.5) +
  geom_point(data = neg_df, aes(colour = Type), size = 4.5) +
  
  scale_colour_manual(
    values = c(
      "Negative control (DDH2O)" = "red",
      "Positive control (10% (w/v) BSA)" = "green"
    ),
    labels = c(
      expression("Negative control (DDH"["2"]*"O)"),
      "Positive control (10% (w/v) BSA)"
    )
  ) +
  
  labs(
    title = "BSA Standard Curve (Quadratic Model)",
    x = "BSA concentration (µg/mL)",
    y = "Blank corrected Absorbance at 595 nm",
    colour = NULL
  ) +
  
  # CRISP, LARGE EQUATION + R²
  annotate("text", x = 50, y = 0.30, label = eq,
           size = 7, fontface = "bold", hjust = 0) +
  annotate("text", x = 50, y = 0.26, label = r2,
           size = 7, fontface = "bold", hjust = 0) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    
    # Journal‑grade axis lines
    axis.line = element_line(colour = "black", linewidth = 1.5),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.width = unit(0.15, "cm"),
    legend.key.height = unit(0.15, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.margin = margin(t = -4)
  )

# ---------------------------------------------------------
# 5. PRINT PLOT
# ---------------------------------------------------------
print(p_bradford_quad)
