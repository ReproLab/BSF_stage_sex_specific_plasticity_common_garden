# Load necessary libraries ------------------------------------------------
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(cowplot)
library(grid)
library(gridExtra)
library(rlang)
library(tidyr)
library(lmodel2)
library(purrr)
library(patchwork)
library(lmerTest)
library(agricolae)

### plotting setup -------------------------------------------------

custom_colors <- c(
  "A"    = "red4",
  "B"    = "midnightblue",
  "C"    = "darkorange2",
  "D"    = "goldenrod",
  "E"    = "darkgreen",
  "WT_1" = "grey70",
  "WT_2" = "grey34",
  "SW"   = "black"
)

custom_shapes <- c(
  "A"    = 1,
  "B"    = 2,
  "C"    = 3,
  "D"    = 4,
  "E"    = 5,
  "WT_1" = 6,
  "WT_2" = 7,
  "SW"   = 8
)

### Data import and basic formatting --------------------------------------

# Larval–pupal dataset
larval_pupal_data <- read.csv(
  "D:/My research/Common garden project/data and analysis/larvae.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    Genetic_line = factor(Genetic_line),
    Treatment    = factor(Treatment),
    Replicate    = factor(Replicate)
  )

# Adult survival data
adult_survival_data <- read.csv(
  "D:/My research/Common garden project/data and analysis/survival.csv"
) %>%
  filter(!is.na(Date.of.Death), !is.na(Sex), Date.of.Death != "") %>%
  mutate(
    Date.of.Death    = as.Date(Date.of.Death, "%d/%m/%Y"),
    Date.of.Eclosion = as.Date(Date.of.Eclosion, "%d/%m/%Y"),
    lifespan         = as.numeric(Date.of.Death - Date.of.Eclosion),
    Group            = factor(Group),
    Genetic_line     = factor(Genetic_line),
    Treatment        = factor(Treatment),
    Sex              = factor(Sex)
  ) %>%
  filter(lifespan != 0)

# Adult morphology / reproductive data
adult_repro_data <- read.csv(
  "D:/My research/Common garden project/data and analysis/repro_morphs.csv"
) %>%
  mutate(
    Genetic_line = factor(Genetic_line),
    Treatment    = factor(Treatment),
    Sex          = factor(Sex)
  )

# Ensure numeric for organ weights
adult_repro_data$ow <- as.numeric(as.character(adult_repro_data$ow))
adult_repro_data$tw <- as.numeric(as.character(adult_repro_data$tw))
adult_repro_data$ag <- as.numeric(as.character(adult_repro_data$ag))

### Generic summariser for multi-trait reaction norms --------------------

summarize_traits <- function(df, trait_cols) {
  df %>%
    pivot_longer(
      cols      = all_of(trait_cols),
      names_to  = "Trait",
      values_to = "Value"
    ) %>%
    group_by(Genetic_line, Treatment, Trait) %>%
    summarize(
      Mean     = mean(Value, na.rm = TRUE),
      SE       = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
      CI_lower = Mean - qt(0.975, df = sum(!is.na(Value)) - 1) * SE,
      CI_upper = Mean + qt(0.975, df = sum(!is.na(Value)) - 1) * SE,
      .groups  = "drop"
    )
}

# Larval/pupal plotting helper (no facet) -------------------------------

plot_reaction <- function(df, trait_name, y_lab) {
  ggplot(dplyr::filter(df, Trait == trait_name),
         aes(x = Treatment,
             y = Mean,
             colour = Genetic_line,
             shape  = Genetic_line,
             group  = Genetic_line)) +
    geom_point(
      position = position_dodge(0.2),
      size     = 2,
      stroke   = 1
    ) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper),
      position = position_dodge(0.2),
      width    = 0.8,
      linewidth = 0.6
    ) +
    geom_line(
      linetype  = "dashed",
      linewidth = 0.9
    ) +
    scale_color_manual(name = "Genetic line", values = custom_colors) +
    scale_shape_manual(name = "Genetic line", values = custom_shapes) +
    labs(
      x = "Diet",
      y = y_lab
    ) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid    = element_blank(),
      strip.text    = element_text(size = 12, face = "bold"),
      axis.text     = element_text(size = 10, color = "black"),
      axis.title    = element_text(size = 12, color = "black"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10),
      legend.position = "right"
    )
}

### Figure 2(a-c): Larval traits (SGR, prepupal weight, larval duration) ---------

summary_larvae <- summarize_traits(
  larval_pupal_data,
  c("SGR_per_day", "avg_pupal_weight", "larval_duration")
)

p_SGR <- plot_reaction(
  summary_larvae,
  "SGR_per_day",
  "Specific growth rate (per day)"
)

p_pupal_weight <- plot_reaction(
  summary_larvae,
  "avg_pupal_weight",
  "Mean prepupal weight (g)"
)

p_larval_duration <- plot_reaction(
  summary_larvae,
  "larval_duration",
  "Larval duration (days)"
)

# To check the plots:
p_SGR
p_pupal_weight
p_larval_duration

### Figure 2(d-f): Pupal traits (pupation rate, pupal duration, eclosion rate) ----

summary_pupal <- summarize_traits(
  larval_pupal_data,
  c("pupation_rate", "pupal_duration", "eclosion_rate")
)

p_pupation_rate <- plot_reaction(
  summary_pupal,
  "pupation_rate",
  "Pupation rate (%)"
)

p_pupal_duration <- plot_reaction(
  summary_pupal,
  "pupal_duration",
  "Pupal duration (days)"
)

p_eclosion_rate <- plot_reaction(
  summary_pupal,
  "eclosion_rate",
  "Eclosion rate (%)"
)

# To check the plots:
p_pupation_rate
p_pupal_duration
p_eclosion_rate

### Mixed linear models for larval and pupal traits ----------------------

model_SGR <- lmer(
  SGR_per_day ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_SGR, type = 3)

model_prepupal_wt <- lmer(
  avg_pupal_weight ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_prepupal_wt, type = 3)

model_larval_duration <- lmer(
  larval_duration ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_larval_duration, type = 3)

model_pupation_rate <- lmer(
  pupation_rate ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_pupation_rate, type = 3)

model_pupal_duration <- lmer(
  pupal_duration ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_pupal_duration, type = 3)

model_eclosion_rate <- lmer(
  eclosion_rate ~ Genetic_line * Treatment + (1 | Replicate),
  data = larval_pupal_data
)
anova(model_eclosion_rate, type = 3)

### Adult stage helpers ---------------------------------------------------

# Summarise a single adult trait by Genetic line, Treatment, Sex
summarize_trait_adult <- function(df, trait) {
  df %>%
    group_by(Genetic_line, Treatment, Sex) %>%
    summarize(
      Mean     = mean(.data[[trait]], na.rm = TRUE),
      SE       = sd(.data[[trait]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[trait]]))),
      CI_lower = Mean - qt(0.975, df = sum(!is.na(.data[[trait]])) - 1) * SE,
      CI_upper = Mean + qt(0.975, df = sum(!is.na(.data[[trait]])) - 1) * SE,
      .groups  = "drop"
    )
}

# Figure 2 (g-h): Plot with facets by Sex, same style as larval/pupal reactions
plot_reaction_adult <- function(df, y_lab) {
  ggplot(df,
         aes(x = Treatment, y = Mean,
             colour = Genetic_line, shape = Genetic_line,
             group = Genetic_line)) +
    geom_point(position = position_dodge(0.2), size = 2, stroke = 1) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_Upper),
                  position = position_dodge(0.2),
                  width = 0.8,
                  linewidth = 0.6) +
    geom_line(linetype = "dashed", linewidth = 0.9) +
    scale_color_manual(name = "Genetic line", values = custom_colors) +
    scale_shape_manual(name = "Genetic line", values = custom_shapes) +
    labs(x = "Diet", y = y_lab) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid    = element_blank(),
      strip.text    = element_text(size = 12, face = "bold"),
      axis.text     = element_text(size = 10, color = "black"),
      axis.title    = element_text(size = 12, color = "black"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10),
      legend.position = "right"
    ) +
    facet_wrap(~ Sex, ncol = 2)
}

#### Figure 2 (g): Adult body weight ----------------------------------------------------

summary_bw <- summarize_trait_adult(adult_repro_data, "Body.weight")

p_bw <- plot_reaction_adult(summary_bw,
                            y_lab = "Body weight (g)")

# To check:
p_bw

model_body <- lmer(
  Body.weight ~ Genetic_line * Treatment * Sex +
    (1 | replicate) + (1 | Operator),
  data = adult_repro_data
)
anova(model_body, type = 3)

#### Figure 2 (h):Adult lifespan -------------------------------------------------------

summary_lifespan <- summarize_trait_adult(adult_survival_data, "lifespan")

p_lifespan <- plot_reaction_adult(summary_lifespan,
                                  y_lab = "Lifespan (days)")

# To check:
p_lifespan

model_lifespan <- lmer(
  lifespan ~ Genetic_line * Treatment * Sex + (1 | Replicate),
  data = adult_survival_data
)
anova(model_lifespan, type = 3)


#### Figure 2 (i-l): Reproductive organs and glands --------------------------------------

# Helper: summarise one trait within a given sex
summarize_trait_sex_specific <- function(data, sex, trait) {
  data %>%
    filter(Sex == sex) %>%
    group_by(Genetic_line, Treatment) %>%
    summarize(
      Mean     = mean(.data[[trait]], na.rm = TRUE),
      SE       = sd(.data[[trait]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[trait]]))),
      CI_lower = Mean - qt(0.975, df = sum(!is.na(.data[[trait]])) - 1) * SE,
      CI_upper = Mean + qt(0.975, df = sum(!is.na(.data[[trait]])) - 1) * SE,
      .groups  = "drop"
    )
}

# Plotting helper for sex-specific traits (no facet)
plot_trait <- function(df, y_label, title = NULL,
                       remove_legend = FALSE, remove_x_label = FALSE) {
  ggplot(df,
         aes(x = Treatment, y = Mean,
             colour = Genetic_line, shape = Genetic_line,
             group  = Genetic_line)) +
    geom_point(position = position_dodge(0.2), size = 2, stroke = 1) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(0.2),
                  width = 0.8,
                  linewidth = 0.6) +
    geom_line(linetype = "dashed", linewidth = 0.9) +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(name = "Genetic line", values = custom_shapes) +
    labs(
      x     = if (remove_x_label) NULL else "Diet",
      y     = y_label,
      title = title
    ) +
    theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(hjust = 0.5),
      legend.position = if (remove_legend) "none" else "right"
    )
}

# Summaries for gonads and glands
summary_ow_f <- summarize_trait_sex_specific(adult_repro_data, "F", "ow")
summary_tw_m <- summarize_trait_sex_specific(adult_repro_data, "M", "tw")
summary_ag_f <- summarize_trait_sex_specific(adult_repro_data, "F", "ag")
summary_ag_m <- summarize_trait_sex_specific(adult_repro_data, "M", "ag")

# Four Figure 2 (i-l) panels (arrange into 2x2 later)
p_ow_f <- plot_trait(summary_ow_f, "Ovaries (mg)",   title = "Female",
                     remove_legend = TRUE,  remove_x_label = TRUE)
p_tw_m <- plot_trait(summary_tw_m, "Testes (mg)",    title = "Male",
                     remove_legend = TRUE,  remove_x_label = TRUE)
p_ag_f <- plot_trait(summary_ag_f, "Female AG (mg)",
                     remove_legend = TRUE,  remove_x_label = FALSE)
p_ag_m <- plot_trait(summary_ag_m, "Male AG (mg)",
                     remove_legend = FALSE, remove_x_label = FALSE)

# To check these:
p_ow_f
p_tw_m
p_ag_f
p_ag_m

# Mixed models for gonads and glands  -------------------------------------

data_female <- subset(adult_repro_data, Sex == "F")
data_male   <- subset(adult_repro_data, Sex == "M")

model_ow_f <- lmer(
  ow ~ Genetic_line * Treatment + (1 | replicate),
  data = data_female
)
anova(model_ow_f, type = 3)

model_tw_m <- lmer(
  tw ~ Genetic_line * Treatment +
    (1 | replicate) + (1 | Operator),
  data = data_male
)
anova(model_tw_m, type = 3)

model_ag_f <- lmer(
  ag ~ Genetic_line * Treatment + (1 | replicate),
  data = data_female
)
anova(model_ag_f, type = 3)

model_ag_m <- lmer(
  ag ~ Genetic_line * Treatment +
    (1 | replicate) + (1 | Operator),
  data = data_male
)
anova(model_ag_m, type = 3)


### Figure 1 - Growth rate ###

# Read the data
data_growth <- read.csv("D:/My research/Common garden project/data and analysis/Growth rate data.csv")

# Check required columns
if (!all(c("line", "diet", "day", "weight_g") %in% colnames(data_growth))) {
  stop("The dataset must contain: 'line', 'diet', 'day', 'weight_g'.")
}

# Exclude missing
if (any(is.na(data_growth))) {
  warning("Missing values detected; they will be removed before summary.")
}

# ----- Set ORDER of facet panels -----
# Put SW and WT_2 in the last row
data_growth$line <- factor(
  data_growth$line,
  levels = c("A","B","C","D","E","WT_1","WT_2","SW")
)

# Summaries
summary_data_growth <- data_growth %>% 
  group_by(line, diet, day) %>% 
  summarise(
    mean_biomass = mean(weight_g, na.rm = TRUE), 
    se_biomass   = sd(weight_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop_last"
  )

# Plot
ggplot(summary_data_growth, aes(
  x = day, 
  y = mean_biomass, 
  linetype = diet, 
  group = interaction(diet, line)
)) +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 1, color = "black") +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass,
                    ymax = mean_biomass + se_biomass),
                width = 0.5, size = 0.5, color = "black") +
  facet_wrap(~ line, scales = "free_y") +
  labs(
    x = "Time (Days)",
    y = "Biomass (g)"
  ) +
  scale_linetype_manual(values = c("CF" = "dashed", "FW" = "solid")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid   = element_blank(),
    strip.text   = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 10, color = "black"),
    axis.title   = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    legend.position     = c(0.95, 0.1),
    legend.justification = c(1.1, 0.6)
  )


####Table 1 GSI and AG to gonad ratio##

if(!"group" %in% names(newdata)){
  newdata <- newdata %>% mutate(group = interaction(Genetic_line, Treatment))
}

get_summary <- function(df, measure) {
  if(all(is.na(df[[measure]]))) {
    warning(paste("Column", measure, "is empty or has only NAs."))
    return(data.frame(group = unique(df$group), mean = NA, se = NA, letter = NA))
  }
  # Diagnostic QQ plot
  qqnorm(df[[measure]]); qqline(df[[measure]])
  
  # Run Bartlett's test (using df$group).
  bart <- bartlett.test(df[[measure]] ~ df$group)
  
  # Compute group-wise mean and SE.
  summ <- df %>% group_by(group) %>% 
    summarise(mean = mean(!!sym(measure), na.rm = TRUE),
              se = sd(!!sym(measure), na.rm = TRUE)/sqrt(n()),
              .groups = "drop")
  
  # If variances are homogeneous, use ANOVA and LSD.test; otherwise, use Kruskal–Wallis.
  if(bart$p.value > 0.05) {
    aov_model <- aov(as.formula(paste(measure, "~ group")), data = df)
    LSDres <- LSD.test(aov_model, "group", group = TRUE)
    letters <- LSDres$groups
  } else {
    kruskal_res <- kruskal(df[[measure]], df$group, group = TRUE, p.adj = "bonferroni")
    letters <- kruskal_res$groups
  }
  letters <- data.frame(group = rownames(letters), letter = letters$groups, stringsAsFactors = FALSE)
  left_join(summ, letters, by = "group")
}

# Separate data by Sex.
data_f <- newdata %>% filter(Sex == "F")
data_m <- newdata %>% filter(Sex == "M")

# For each sex, we compare genetic lines independently within each treatment.
# For females:
data_f_CF <- data_f %>% filter(Treatment == "CF") %>% mutate(group = Genetic_line)
data_f_FW <- data_f %>% filter(Treatment == "FW") %>% mutate(group = Genetic_line)

summ_f_GSI_CF <- get_summary(data_f_CF, "GSI") %>% mutate(Sex = "F", Treatment = "CF", Measure = "GSI")
summ_f_ag_CF  <- get_summary(data_f_CF, "ag_gonad_ratio") %>% mutate(Sex = "F", Treatment = "CF", Measure = "ag_to_gonad")
summ_f_GSI_FW <- get_summary(data_f_FW, "GSI") %>% mutate(Sex = "F", Treatment = "FW", Measure = "GSI")
summ_f_ag_FW  <- get_summary(data_f_FW, "ag_gonad_ratio") %>% mutate(Sex = "F", Treatment = "FW", Measure = "ag_to_gonad")

# For males:
data_m_CF <- data_m %>% filter(Treatment == "CF") %>% mutate(group = Genetic_line)
data_m_FW <- data_m %>% filter(Treatment == "FW") %>% mutate(group = Genetic_line)

summ_m_GSI_CF <- get_summary(data_m_CF, "GSI") %>% mutate(Sex = "M", Treatment = "CF", Measure = "GSI")
summ_m_ag_CF  <- get_summary(data_m_CF, "ag_gonad_ratio") %>% mutate(Sex = "M", Treatment = "CF", Measure = "ag_to_gonad")
summ_m_GSI_FW <- get_summary(data_m_FW, "GSI") %>% mutate(Sex = "M", Treatment = "FW", Measure = "GSI")
summ_m_ag_FW  <- get_summary(data_m_FW, "ag_gonad_ratio") %>% mutate(Sex = "M", Treatment = "FW", Measure = "ag_to_gonad")

# Combine all summaries into one final table.
final_table_2 <- bind_rows(summ_f_GSI_CF, summ_f_ag_CF, summ_f_GSI_FW, summ_f_ag_FW,
                           summ_m_GSI_CF, summ_m_ag_CF, summ_m_GSI_FW, summ_m_ag_FW) %>% 
  mutate(mean_se = paste0(round(mean, 2), " ± ", round(se, 2))) %>% 
  dplyr::select(Sex, Treatment, Measure, group, mean_se, letter)

print(final_table_2, n = Inf)