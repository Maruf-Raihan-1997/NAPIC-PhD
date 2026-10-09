#script that:
  
#Creates the dataset
#Runs the two-way ANOVA
#Checks assumptions
#Performs Tukey post hoc
#Generates Tukey letters
#Creates the thesis-standard grouped bar plot
#Creates the thesis-standard interaction plot
# =====================================================================
# FREE AMINO ACID RELEASE
# TWO-WAY ANOVA + TUKEY LETTERS
# THESIS STANDARD FORMAT
# =====================================================================

library(tidyverse)
library(car)
library(emmeans)
library(multcompView)

# =====================================================================
# DATA
# =====================================================================

codes <- c("ND", "W", "Bn", "Br")

labels <- c(
  ND = "Undigested DWP",
  W  = "10% (w/v) DWP Water",
  Bn = "10% (w/w) DWP Banana",
  Br = "10% (w/w) DWP Bread"
)

cols <- c(
  ND = "grey40",
  W  = "steelblue",
  Bn = "forestgreen",
  Br = "saddlebrown"
)

dat <- tribble(
  ~amino_acid, ~treatment, ~r1, ~r2, ~r3,
  
  "Glutamic acid","ND",0.44,0.45,0.48,
  "Glutamic acid","W",0.23,0.21,0.36,
  "Glutamic acid","Bn",0.05,0.05,0.08,
  "Glutamic acid","Br",0.04,0.04,0.05,
  
  "Aspartic acid","ND",0.96,0.75,0.79,
  "Aspartic acid","W",0.85,NA,0.92,
  "Aspartic acid","Bn",0.86,1.03,0.79,
  "Aspartic acid","Br",0.99,NA,1.27,
  
  "Leucine","ND",0.22,0.16,0.21,
  "Leucine","W",0.45,0.45,0.44,
  "Leucine","Bn",0.16,0.20,0.14,
  "Leucine","Br",0.17,0.16,0.17
) %>%
  pivot_longer(
    r1:r3,
    names_to = "rep",
    values_to = "conc"
  ) %>%
  filter(!is.na(conc)) %>%
  mutate(
    amino_acid = factor(
      amino_acid,
      levels = c(
        "Glutamic acid",
        "Aspartic acid",
        "Leucine"
      )
    ),
    treatment = factor(
      treatment,
      levels = codes
    )
  )

# =====================================================================
# TWO-WAY ANOVA
# =====================================================================

fit2 <- aov(
  conc ~ treatment * amino_acid,
  data = dat
)

cat("\n=== TWO-WAY ANOVA ===\n")
print(summary(fit2))

# =====================================================================
# ASSUMPTIONS
# =====================================================================

cat("\n=== SHAPIRO TEST ===\n")
print(shapiro.test(residuals(fit2)))

cat("\n=== LEVENE TEST ===\n")
print(
  leveneTest(
    conc ~ interaction(treatment, amino_acid),
    data = dat
  )
)

# =====================================================================
# TUKEY LETTERS
# =====================================================================

emm <- emmeans(
  fit2,
  ~ treatment | amino_acid
)

letters_df <- multcomp::cld(
  emm,
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame()

letters_df$.group <- gsub(
  " ",
  "",
  letters_df$.group
)

# =====================================================================
# SUMMARY DATA
# =====================================================================

plot2 <- dat %>%
  group_by(
    amino_acid,
    treatment
  ) %>%
  summarise(
    mean = mean(conc),
    sd = sd(conc),
    .groups = "drop"
  ) %>%
  left_join(
    letters_df %>%
      dplyr::select(
        amino_acid,
        treatment,
        .group
      ),
    by = c(
      "amino_acid",
      "treatment"
    )
  )

# =====================================================================
# BAR PLOT
# =====================================================================

dodge <- position_dodge(0.7)

p_twoway <- ggplot(
  plot2,
  aes(
    x = amino_acid,
    y = mean,
    fill = treatment
  )
) +
  
  geom_col(
    position = dodge,
    width = 0.7,
    colour = "white"
  ) +
  
  geom_errorbar(
    aes(
      ymin = pmax(mean - sd, 0),
      ymax = mean + sd
    ),
    position = dodge,
    width = 0.2,
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = .group,
      y = mean + sd + 0.08
    ),
    position = dodge,
    size = 10,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = cols
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20))
  ) +
  
  labs(
    title = "Free Amino Acid Release (mg/g)",
    x = "Amino Acid",
    y = "Free Amino Acid Concentration (mg/g)"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 5),
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 25,
      face = "bold",
      colour = "black",
    ),
    
    axis.title.x = element_text(
      size = 25,
      face = "bold",
      margin = margin(t = 20)
    ),
    
    axis.title.y = element_text(
      size = 22,
      face = "bold",
      margin = margin(r = 20)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 2.5
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 2.5
    ),
    
    axis.ticks.length = unit(
      0.25,
      "cm"
    )
  )

print(p_twoway)

## =====================================================================
# INTERACTION PLOT
# =====================================================================

interaction_labels <- c(
  ND = "Undigested\nDWP",
  W  = "10% (w/v)\nDWP Water",
  Bn = "10% (w/w)\nDWP Banana",
  Br = "10% (w/w)\nDWP Bread"
)

interaction_summary <- dat %>%
  group_by(
    amino_acid,
    treatment
  ) %>%
  summarise(
    mean = mean(conc),
    sd = sd(conc),
    .groups = "drop"
  )

p_interaction <- ggplot(
  interaction_summary,
  aes(
    x = treatment,
    y = mean,
    colour = amino_acid,
    group = amino_acid
  )
) +
  
  geom_line(
    linewidth = 2
  ) +
  
  geom_point(
    size = 6
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    width = 0.08,
    linewidth = 1
  ) +
  
  coord_cartesian(
    ylim = c(0, NA)
  ) +
  
  scale_x_discrete(
    labels = interaction_labels
  ) +
  
  scale_colour_manual(
    values = c(
      "Glutamic acid" = "#D55E00",
      "Aspartic acid" = "#0072B2",
      "Leucine"       = "#009E73"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  labs(
    title = "Food Matrix × Free Amino Acid Interaction",
    x = "Food Matrix",
    y = "Free Amino Acid Concentration (mg/g)",
    colour = "Amino Acid"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 5),
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 25,
      face = "bold",
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 25,
      face = "bold",
      margin = margin(t = 20)
    ),
    
    axis.title.y = element_text(
      size = 22,
      face = "bold" , 
      margin = margin(r = 20)
    ),
    
    legend.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 18
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 2.5
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 2.5
    ),
    
    axis.ticks.length = unit(
      0.25,
      "cm"
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_interaction)