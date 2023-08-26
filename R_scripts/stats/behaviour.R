library(tidyverse) # Handling large datasets
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(emmeans) # Comparisons among groups
library(fishualize) # Nice fish colours
library(bayestestR) # Gives access to the describe_posterior() function
# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")

raw_data <- readxl::read_excel("Overwinter Data 2018 Master Sheet.xlsx", sheet = "Sheet1") %>% 
  filter(!is.na(`TimeInCue(%)`)) %>% 
  rename(Timepoint = `Time Point`, Tank_ID = `Tank#`, Cue_Proportion = `TimeInCue(%)`, Trial_ID = `Trial ID`) %>% 
  mutate(Treatment = case_when(Treatment == "pCO2" ~ "pCO2",
                               Treatment == "Temp + PCO2" ~ "Combo", 
                               Treatment == "control" ~ "Control",
                               Treatment == "temp" ~ "Temp")) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Combo", "pCO2", "Temp", "Control")), Timepoint = factor(Timepoint, levels = c("Zero", "24hr"))) %>% 
  mutate(treatment_combined = paste0(Treatment, "_", Timepoint)) %>% 
  mutate(Date_Numeric = as.numeric(Date)) %>% 
  mutate(Date_Factor = as.factor(case_when(Date == 1518134400 ~ "Feb09",
                                           Date == 1518393600 ~ "Feb12",
                                           Date == 1518480000 ~ "Feb13",
                                           Date == 1518652800 ~ "Feb15"))) %>% 
  select(-Date_Numeric) %>% 
  mutate(Tank_ID = as.factor(letters[Tank_ID])) %>% 
  filter(is.na(Notes)) %>% 
  mutate(Cue_Proportion_Fraction = Cue_Proportion/100) %>% 
  mutate(treatment_combined = factor(treatment_combined, levels = c("Control_Zero", "Control_24hr", "Temp_Zero", "Temp_24hr", "pCO2_Zero", "pCO2_24hr", "Combo_Zero", "Combo_24hr")))

raw_data$treatmemt_combined

# Take a look at the raw data
ggplot(data = raw_data, aes(x = Timepoint, y = Cue_Proportion)) +
  geom_point() +
  facet_wrap(~Treatment) +
  theme_bw()


cue_model <- brm(data = raw_data, family = skew_normal,
                      bf(Cue_Proportion ~ (Treatment | Timepoint) + (1 | CueSide / Tank_ID) + (1 | Date_Factor),
                         sigma ~ treatment_combined), 
                      prior = c(prior(normal(50, 40), class = Intercept)),
                      iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                      control = list(adapt_delta = 0.99, max_treedepth = 12))

summary(cue_model)

get_variables(cue_model)

hypothesis(cue_model, "treatment_combinedCombo_24hr < treatment_combinedCombo_Zero")
hypothesis(cue_model, hypothesis = "treatment_combinedTemp_24hr > treatment_combinedTemp_Zero", class = "b")


describe_posterior(cue_model, centrality = "mean", dispersion = TRUE, test = "bayesfactor")

# Look at the posterior distribution. It was left skewed when using a Gaussian family distribution, but was fixed by using skew_normal instead
pp_check(cue_model, ndraws = 100) +
  theme_bw(base_size = 20)

plot(conditional_effects(cue_model), points = T)



# Use gather emmeans draws to plot the model results, showing different posterior distributions for each treatment at each timepoint
gather_emmeans_draws(emmeans(cue_model, ~ Treatment | Timepoint, re_formula = NULL)) |>
  ggplot( aes(x = .value, y = Timepoint, fill = Treatment, colour = Treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab("Proportion of Time in Cue (%)") +
  ylab("Timepoint (hours)") +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  #xlim(-10, 100) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.91, 0.14), legend.background = element_blank())

emmeansdraws <- gather_emmeans_draws(emmeans(cue_model, ~ Treatment | Timepoint, re_formula = NULL)) |>
  #separate(treatment_combined, into = c("treatment", "time"), sep = "_") |> 
  mutate(Treatment = factor(Treatment, levels = c("Combo", "pCO2", "Temp", "Control"))) |> 
  mutate(Timepoint = factor(Timepoint, levels = c("Zero", "24hr"))) 


cue_plot <- ggplot(data = emmeansdraws, aes(x = .value, y = Timepoint, fill = Treatment, colour = Treatment)) +
  stat_halfeye(orientation = "horizontal", alpha = 0.5) +
  scale_fill_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  scale_colour_fish_d(option = "Lampris_guttatus", direction = -1, labels = c(bquote(italic(p)*CO[2]*"+"*Temp), bquote(italic(p)*CO[2]), "Temperature", "Control")) +
  xlab("Proportion of Time in Alarm Cue (%)") +
  ylab("Timepoint (hours)") +
  scale_y_discrete(labels = c("0", "24")) +
  guides(fill = guide_legend(title="Treatment"), colour = guide_legend(title="Treatment")) +
  #xlim(-10, 100) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.87, 0.18), legend.background = element_blank(), plot.margin = margin(t = 5.6, r = 14, b = 5.6, l = 5.6, "points"))

ggsave(filename = "cue_plot.pdf", plot = cue_plot, dpi = 2000)


# Pull pairwise results from the emmeans function with tidybayes, and make halfeye plots
plot(pairs(emmeans(test, ~ Treatment | Timepoint, re_formula = NULL)))
pairwise_table <- pairs(emmeans(cue_model, ~ Treatment | Timepoint, re_formula = NULL))
pairwise_table <- summary(pairwise_table)
write_tsv(x = pairwise_table, file = "2018_overwinter_cue_emmeans.txt")


# pairwise_groups <- emmeans(cue_model, pairwise ~ treatment_combined)
# 
# bayesfactor_parameters(pairwise_groups, prior = cue_model)
# 
# prior_model <- unupdate(cue_model)
# prior_emmgrid <- emmeans(prior_model, pairwise ~ treatment_combined)
# modelbased::estimate_contrasts(cue_model, test = "bf", bf_prior = prior_emmgrid)
