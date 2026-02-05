library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
###########  Dot Matrix Plot
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)
library(ggplot2)

assign_group <- function(x) {
  case_when(
    str_detect(x, regex("TestAge_10|ethnicity|charlson|IMDPercentile|LinkedSex", TRUE)) ~ "Demographics",
    str_detect(x, regex("InfecSource|csu|multi_infection|hospital_stays|icu|frailty", TRUE)) ~ "Healthcare & clinical history",
    str_detect(x, regex("time_to_test_yearly|chemo|total_neut_10|cancer_|\\bSex\\b", TRUE)) ~ "Cancer related",
    #str_detect(x, regex("charlson", TRUE)) ~ "Comorbidities",
    str_detect(x, regex("days_on_any_abx|days_on_tested_abx|time_gap_30", TRUE)) ~ "Antibiotic exposure",
    str_detect(x, regex(
      "prior_resistance", TRUE)) ~ "Previous resistance \nin corresponding \npathogens in any culture",
    str_detect(x, regex("num_uc_10\\b|num_uc_pt_bug\\b|total\\b|num_prev_bld\\b", TRUE)) ~
      "Microbiological cultures sent",
    x == "cld_days_yearly" ~ "Study timeline",
    TRUE ~ "Other"
  )
}

pretty_label <- function(x) {
  # Exact term → pretty label 
  explicit_map <- c(
    "IMDPercentile"                     = "Most recent IMD percentile",
    "days_on_tested_abx"            = "Per 1 additional day on the target antimicrobial in 1y",
    "num_uc_10"                     = "Per 10 additional UCs taken in 1y",
    "prev_bug_R_drug_bld"           = "The target drug R BC isolates in corresponding bug in 1y",
    "prev_bug_R_drug_non_bld"       = "The target drug R non-BC isolates in corresponding bug in 1y",
    "prev_bug_S_drug_non_bld"       = "The target drug S non-BC isolates in corresponding bug in 1y",
    "prev_bug_R_other_drug_non_bld" = "The non-target drug R non-BC isolates in corresponding bug in 1y",
    "prev_bug_S_other_drug_non_bld" = "The non-target drug S non-BC isolates in corresponding bug in 1y",
    "prev_bug_other_drug_R_bld" = "The non-target drug R BC isolates in corresponding bug in 1y",
    "prev_bug_other_drug_S_bld"  = "The non-target drug S BC isolates in corresponding bug in 1y",
    "prev_bug_S_drug_bld"           = "The target drug S BC isolates in corresponding bug in 1y",
    "cld_days_yearly"               = "Per 1y since 01/04/2015",
    "ethnicityNon-white"            = "Non-white vs white ethnicity",
    "ethnicityWhite"                = "White ethnicity",
    "ethnicityMissed"               = "Missing vs white ethnicity",
    "multi_infection"               = "Multiple infection (bld) in 1y",
    "time_gap_30"                    = "Per 1m since last antibiotic used in 1y",
    "ns_time_gap_30_1"                 = "Months since last antibiotic used in 1y",
    "ns_time_gap_30_2"                 = "Months since last antibiotic used in 1y",
    "num_prev_bld"                  = "Per 1 additional BC taken in 1y",
    "active_chemo"                  = "Active chemotherapy in 1y",
    "time_to_test_yearly"           = "Per 1y since the first cancer diagnostic code",
    "cancer_groupLymphoid/Haematopoietic" = "Lymphoid/Haematopoietic cancer",
    "days_since_last_active_chemo_30"     = "Months since last chemo in 1y",
    "ns(days_since_last_active_chemo_30)"     = "Months since last chemo in 1y",
    "days_on_any_abx_10" = "Per 10 additional days on non-target antimicrobials in 1y",
    "hospital_stays_10" = "Per 10 additiona days in hospital in 1y",
    "active_chemo_days_30" =  "Months on chemo in 1y",
    "LinkedSexM" = "Male",
    "LinkedSexF" = "Female vs Male",
    "total_neut_10" = "Per 10 additiona neutrophil tests in 1y",
    "num_uc_pt_bug" = "Number of positive UCs on this bug in 1y",
    "icu_stays" = "Per 1 additiona day in ICU in 1y"
  )
  
  y <- dplyr::recode(x, !!!explicit_map, .default = x)
  
  # General rules (apply only where not already replaced)
  # cancer_groupFoo -> Foo cancer
  y <- str_replace_all(y, fixed("cancer_groupLymphoid/Haematopoietic"), "Lymphoid/Haematopoietic vs colorectal cancer")
  mask_cancer <- str_detect(x, "^cancer_group") & (y == x)
  y[mask_cancer] <- paste0(str_replace(x[mask_cancer], "^cancer_group", ""), " vs colorectal cancer")
  y <- str_replace_all(y, fixed("Lymphoid/Haematopoietic cancer"), "Lymphoid/Haematopoietic vs colorectal cancer")
  # Specific semantic renames
  y <- str_replace_all(y, fixed("TestAge_10"), "Per 10y increase in age at the bacteraemia")
  y <- str_replace_all(y, fixed("InfecSourceHospital"),  "Hospital vs community infection source")
  y <- str_replace_all(y, fixed("InfecSourceNone"),      "Unknown vs community infection source")
  y <- str_replace_all(y, fixed("InfecSourceCommunity"), "Community infection source")
  y <- str_replace_all(y, fixed("charlson"), "Per 1 unit higher charlson score")
  y <- str_replace_all(y, fixed("frailty_score_10"), "Per 10 unit higher frailty scores")
  # Acronyms / tokens → readable text
  y <- str_replace_all(y, regex("\\bicu_stays\\b", ignore_case = TRUE), "Days in ICU in 1y")
  y <- str_replace_all(y, regex("\\bcsu\\b", ignore_case = TRUE), "Postive CSU in 1y")
  # Cosmetic clean-up
  y <- str_replace_all(y, "_", " ")
  y <- str_squish(y)
  y <- str_replace_all(y, fixed("prior resistanceR target"),      "R to target antimicrobial vs no AST/positive isolates in 1y")
  y <- str_replace_all(y, fixed("prior resistanceR other only"),      "Fully S to target, R to any non-target antimicrobials vs \nno AST/positive isolates in 1y")
  y <- str_replace_all(y, fixed("prior resistanceFully S"),      "S to all tested antimicrobials vs no AST/positive isolates in 1y")
  
  as.character(y)
}

factor_lookup <- factor_master %>%
  transmute(
    term,
    group = assign_group(term),
    label = pretty_label(term)
  )
desired_order <- c("TestAge_10", "ethnicity",
                   "num_prev_bld",
                   "num_uc_10"
)

factor_lookup <- factor_lookup %>%
  mutate(
    term = factor(
      term,
      levels = c(desired_order, setdiff(unique(term), desired_order))
    ),
    label = factor(label, levels = unique(factor_lookup$label))
  )
factor_lookup <- factor_lookup %>% arrange(term)

# Identify baseline terms (OR == NA in combined_df_1 after dropping NA term)
baseline_lookup <- factor_master %>%
  filter(is.na(OR)) %>%
  select(term) %>%
  left_join(factor_lookup, by = "term")
baseline_terms <- baseline_lookup$term

#  Bind all robust_* (sig-only; adjust rename() if your columns differ)  
sig_list <- list(
  "Co-amoxiclav & Enterobacterales" = tab_cs2,
  "Fluoroquinolone & Enterobacterales" = tab_cs3,
  "Trimethoprim/Sulfamethoxazole & Enterobacterales" = tab_cs4,
  "Third-generation cephalosporin & Enterobacterales" = tab_cs1,
  "Gentamicin & Enterobacterales" = tab_cs5,
  "Piperacillin-tazobactam & Enterobacterales" = tab_cs6,
  "Vancomycin & Enterococcus" = tab_cs7
) 

sig <- imap_dfr(sig_list, ~
                  .x %>%
                  rename(term = term, OR = OR, p = p.value) %>%
                  mutate(combination = .y) %>%
                  select(combination, term, OR, conf.low, conf.high, p) %>%
                  mutate(
                    term = case_when(
                      str_detect(coalesce(term, ""), fixed("ns_time_gap_30_1")) ~ "ns_time_gap_30_1",
                      str_detect(coalesce(term, ""), fixed("ns_time_gap_30_2")) ~
                        "ns_time_gap_30_2",
                      TRUE ~ term
                    )
                  )
)

combo_levels <- names(sig_list)

# Full grid: every factor × combination; keep NA where nothing happened ──
grid <- expand_grid(
  term        = factor_lookup$term,     # full factor list (includes baselines)
  combination = combo_levels
) %>%
  left_join(factor_lookup, by = "term") %>%                 # add group + label
  left_join(sig, by = c("term", "combination")) %>%         # bring OR/p (sig-only)
  mutate(
    label       = factor(label, levels = rev(levels(factor_lookup$label))),
    combination = factor(combination, levels = combo_levels),
    group       = factor(group, levels = c('Demographics', 'Cancer related', 'Comorbidities', 'Healthcare & clinical history', 'Antibiotic exposure', 'Previous resistance \nin corresponding \npathogens in any culture', 'Microbiological cultures sent', 'Study timeline')),
    LogOR = log(OR),
    
    # Discrete p bins (only these get dots)
    p_bin = case_when(
      is.na(p) ~ NA_character_,
      p < 0.001 ~ "<0.001",
      p >= 0.001 & p < 0.01 ~ "[0.001, 0.01)",
      p >= 0.01  & p <= 0.05  ~ "[0.01, 0.05]",
      p > 0.05  ~ ">0.05",
      TRUE ~ NA_character_           # also keep p > 0.05  
    ),
    p_bin = factor(
      p_bin,
      levels = c(">0.05", "[0.01, 0.05]", "[0.001, 0.01)", "<0.001"),
      ordered = TRUE
    )
  ) %>% filter(term != "(Intercept)")



grid <- grid %>%
  mutate(
    reln = ifelse(str_detect(term, "ns_"), "Non-linear", "Linear"),
    reln = factor(reln, levels = c("Linear", "Non-linear"))
  ) %>%
  filter(term != "LinkedSexM", term != "cancer_groupColorectal", term != "ethnicityWhite", term != "InfecSourceCommunity", term != "prior_resistance")



grid <- grid %>%
  mutate(
    .tag = case_when(
      str_detect(term, fixed("ns_time_gap_30_1")) ~ "time_gap_30",
      str_detect(term, fixed("ns_time_gap_30_2")) ~ "time_gap_30",
      TRUE ~ NA_character_
    ),
    .row = row_number()              # to restore original order later
  ) %>%
  group_by(combination, .tag) %>%
  filter(is.na(.tag) | row_number() == 1) %>%  # keep first per tag/combination
  ungroup() %>%
  arrange(.row) %>%
  select(-.tag, -.row)



grid <- grid %>%
  group_by(label) %>%
  filter(!all(is.na(LogOR) & is.na(p_bin))) %>%
  ungroup()


############# No nonlinear effects 
library(dplyr)
grid_lin    <- grid %>% filter(reln == "Linear")
p_sig  <- grid_lin %>% filter(!is.na(p_bin), p_bin != ">0.05", !is.na(LogOR))
grid_nonlin <- grid %>% filter(reln == "Non-linear", !is.na(LogOR), !is.na(p_bin))
p_nsig <- grid_lin %>% filter(!is.na(p_bin), p_bin == ">0.05", !is.na(LogOR))
p_all <- grid %>% filter(!is.na(p_bin), !is.na(LogOR))
lo <- min(p_all$LogOR, na.rm = TRUE)
hi <- max(p_all$LogOR, na.rm = TRUE)
stopifnot(is.finite(lo), is.finite(hi))

# small band around 0
eps <- max(1e-6, 0.02 * max(abs(c(lo, hi))))

# colours: neg (blue) -> pos (red)
cols  <- c("#0571b0", "#cfe8ff", "#fde0e6", "#ca0020")
stops <- c(lo, -eps, eps, hi)

# clip stops to [lo, hi] so values stay in [0,1] even if 0 is outside range
stops  <- pmin(pmax(stops, lo), hi)

# >>> key fix: call the right rescale
values <- scales::rescale(stops, to = c(0, 1), from = c(lo, hi))
values <- pmin(pmax(values, 0), 1)  # just in case

# if you have only positives or only negatives, you can simplify:
if (lo >= 0) { cols <- c("#fde0e6", "#ca0020"); values <- c(0, 1) }
if (hi <= 0) { cols <- c("#0571b0", "#cfe8ff"); values <- c(0, 1) }

# then in ggplot:
# scale_fill_gradientn(colors = cols, values = values, limits = c(lo, hi), oob = scales::squish)


full_figure <- ggplot() +
  geom_blank(
    data = grid %>% dplyr::distinct(combination, label, group),
    aes(x = combination, y = label)
  ) +
  # Significant: filled by LogOR, white border
  geom_point(
    data = p_sig,
    aes(x = combination, y = label, shape = p_bin, size = p_bin, fill = LogOR),
    colour = "white", stroke = 0.7, alpha = 0.95
  ) +
  # Non-significant: different shape, white fill, grey border
  geom_point(
    data = p_nsig,
    aes(x = combination, y = label, shape = p_bin, size = p_bin, fill = LogOR),
    colour = "white", stroke = 0.7, alpha = 0.95
  ) +
  # Shapes: 21 = filled circle; 24 = filled triangle (pick any 21–25 for filled symbols)
  scale_shape_manual(
    values = c(
      ">0.05"         = 21,
      "[0.01, 0.05]"  = 21,
      "[0.001, 0.01)" = 21,
      "<0.001"        = 21
    ),
    drop = FALSE, name = "p-value"
  ) +
  scale_size_manual(
    values = c(
      ">0.05"         = 2,
      "[0.01, 0.05]"  = 4,
      "[0.001, 0.01)" = 6,
      "<0.001"        = 8
    ),
    drop  = FALSE, name = "p-value",
  ) +
  guides(shape = guide_legend(override.aes = list(fill = "white", colour = "grey60", stroke = 0.7)),
         size = guide_legend(override.aes = list(fill = "white", colour = "grey60", stroke = 0.7 ))
  )+
  # Fill gradient only applies to the significant layer (since ns layer sets fill outside aes)
  scale_fill_gradientn(
    colours = cols,
    values  = values,          # near-discontinuous break around 0
    limits  = c(lo, hi),
    oob     = squish,
    name    = "Log of odds ratio"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Pathogen-antimicrobial combination", y = NULL) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )
print(full_figure)
ggsave("/home/anniehu/AMR Cancer/full figure.png", full_figure, width = 15, height = 11, dpi = 300)

################ plot only the previous resistance variables
# --- the five variables to keep (exact string matches) ---
vars_to_keep <- c( "prior_resistanceR_target", 
                   "prior_resistanceR_other_only", "prior_resistanceFully_S"
                   
)
# Prepare the plotting data
plot_data <- expand_grid(
  term        = factor_lookup$term,     # full factor list (includes baselines)
  combination = combo_levels
) %>%
  left_join(factor_lookup, by = "term") %>%                 # add group + label
  left_join(sig, by = c("term", "combination")) %>%         # bring OR/p (sig-only)
  mutate(
    label       = factor(label, levels = rev(unique(factor_lookup$label))),
    combination = factor(combination, levels = combo_levels),
    group       = factor(group, levels = c('Demographics', 'Cancer related', 'Healthcare & clinical history', 'Antibiotic exposure', 'Previous resistance \nin corresponding \npathogens in any culture', 'Microbiological cultures sent', 'Study timeline')),
  ) %>% filter(term != "(Intercept)")%>% filter(term != "ns_time_gap_30_1") %>% filter(term != "ns_time_gap_30_2") %>%
  filter(term %in% vars_to_keep) %>%
  # ensure 'vars' facets are in the requested order
  mutate(term = factor(as.character(term), levels = vars_to_keep),
         # ensure the y-axis (label) orders the 7 combos consistently within each facet
         combination = factor(as.character(combination), levels = rev(combo_levels))
  )

#Make comparisons updated plot plot only the previous resistance variables ---
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

# Mapping from term to panel_label 
panel_map <- tibble(
  term = c(
    "prior_resistanceR_target", 
    "prior_resistanceR_other_only", "prior_resistanceFully_S"
  ),
  panel_label = c("Resistance to the target antimicrobial vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y",
                  "Fully susceptible to target, resistant to any non-target antimicrobials vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y",
                  "Susceptibility to all tested antimicrobials vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y"
                  
  )
)

# Prepare plotting dataframe 
plot_plotting <- plot_data %>%
  inner_join(panel_map, by = "term") %>%
  mutate(
    panel_label = factor(panel_label,
                         levels = c(
                           "Resistance to the target antimicrobial vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y",
                           "Fully susceptible to target, resistant to any non-target antimicrobials vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y",
                           "Susceptibility to all tested antimicrobials vs no AST/positive isolates, \nin corresponding pathogens in any cultures in 1y"
                         )),
    combination2 = factor(as.character(combination), levels = rev(combo_levels))
  ) %>%
  filter(!is.na(OR) & !is.na(conf.low) & !is.na(conf.high))

#  Create centers for each combination WITHIN each panel 
centers <- plot_plotting %>%
  group_by(panel_label) %>%
  distinct(combination2) %>%
  arrange(combination2) %>%
  mutate(center_idx = row_number()) %>%
  ungroup()

# Merge centers back into plotting data
offset <- 0.15  # Smaller offset = tighter spacing within each combination
plot_for_plot <- plot_plotting %>%
  inner_join(centers, by = c("panel_label", "combination2")) %>%
  mutate(
    y = center_idx + offset
  )

# y-axis labels (one per combination center) 
label_df <- centers %>%
  arrange(panel_label, center_idx)

# Color mapping
cols <- c("#2C7FB8")

# Plot
p_final <- ggplot(plot_for_plot, aes(x = OR, y = y, colour = cols)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.08, linewidth = 0.8, alpha = 0.95) +
  geom_point(size = 2.0) +
  scale_color_manual(values = cols, labels = NULL, name = NULL) +
  facet_wrap(~ panel_label, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
    labels = c("0.1","0.2","0.5","1","2","5","10", "20", "50"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = function(limits) {
      # Get unique center indices for this panel
      unique(label_df$center_idx[label_df$center_idx >= limits[1] & label_df$center_idx <= limits[2]])
    },
    labels = function(breaks) {
      # Return combination labels for these breaks
      sapply(breaks, function(b) {
        combo <- label_df$combination2[label_df$center_idx == b][1]
        as.character(combo)
      })
    },
    expand = expansion(add = c(0.4, 0.4))
  ) +
  labs(x = "Odds ratio (95% CI)", y = NULL) +
  theme_minimal(base_size = 11) +
  guides(colour = "none") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 9, hjust = 0),
    axis.title.x = element_text(size = 9),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(8, 12, 8, 12)
  ) 


print(p_final)