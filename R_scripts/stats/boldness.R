library(tidyverse) # Handling large datasets
library(brms) # Run the Bayesian models
library(tidybayes) # Plot Bayesian models
library(emmeans) # Comparisons among groups
library(fishualize) # Nice fish colours
library(bayestestR) # Gives access to the describe_posterior() function
# Faster Hamiltonian Monte Carlo sampling
options(brms.backend = "cmdstanr")


# Read in the raw data. Replace + signs with _ to keep brms happy, create a new combined treatment column, and turn tank IDs into letters so that brms does not think they are numeric values
raw_data <- readxl::read_excel("Individual Mastersheet Feb 10 2020_MTadditions.xlsm", sheet = "Sheet1") %>% 
  mutate_if(is.character, str_replace_all, pattern = "\\+", replacement = "_") %>% 
  mutate(treatment_combined = paste0(Fish_Treatment, "_", Trial_Treatment)) %>% 
  mutate(Tank = as.factor(letters[Tank]), TreatTankID = as.factor(letters[TreatTankID]))

# Take a quick look at the different data distributions
hist(raw_data$Time_Near_Total)
hist(raw_data$Time_Thigmotaxis_Total)
hist(raw_data$Total_Activity_Total)
hist(raw_data$Activity_Near_Total)

# Take a look at what priors would look for the model
get_prior(data = raw_data, family = skew_normal,
          bf(Total_Activity_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g),
             sigma ~ treatment_combined))

# Run a model on total activity
totalactivity_model <- brm(data = raw_data, family = skew_normal,
                 bf(Total_Activity_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g),
                    sigma ~ treatment_combined), 
                 prior = c(prior(normal(20, 50), class = Intercept),
                           #prior(normal(-50, 60), class = b, coef = "Weight_g"),
                           prior(normal(20, 50), class = b),
                           prior(cauchy(0, 5), class = sd)),
                 iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                 control = list(adapt_delta = 0.99, max_treedepth = 12))

# Look at model results
summary(totalactivity_model)

# Look at model variables
get_variables(totalactivity_model)

# A posterior predictive check to see if the model-generated data matched the actual data
pp_check(totalactivity_model, ndraws = 100) +
  theme_bw(base_size = 20)

# Plot model results
plot(totalactivity_model, ask = F)
plot(conditional_effects(totalactivity_model), points = T, theme = theme_bw())

# Take a look at emmeans
pairs(emmeans(totalactivity_model, ~ treatment_combined | Cue, re_formula = NULL))
plot(emmeans(totalactivity_model, ~ treatment_combined | Cue, re_formula = NULL))

# Pull model posterior draws with gather emmeans draws
totalactivity_draws <- gather_emmeans_draws(emmeans(totalactivity_model, ~ treatment_combined | Cue, re_formula = NULL)) %>% 
  mutate(Cue = factor(Cue, levels = rev(c("Blank", "Food", "Alarm", "Food_Alarm"))))

# Plot posterior 95% credible intervals from the model draws
totalactivity_plot <- ggplot(totalactivity_draws, aes(x = .value, y = treatment_combined, colour = Cue, fill = Cue)) + 
  stat_pointinterval(position = position_dodge(width = 1)) +
  scale_fill_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  xlab("Total Activity (Total # of Body Length Changes in Position)") +
  ylab("Rearing Group, Trial Group") +
  scale_y_discrete(labels = c("Temp_Temp" = "Temp, Temp", "Temp_Control" = "Temp, Control", "pCO2_Temp_pCO2_Temp" = bquote(italic(p)*CO[2]*"+Temp,"~italic(p)*CO[2]*"+"*Temp), "pCO2_Temp_Control" = bquote(italic(p)*CO[2]*"+"*'Temp, Control'), "pCO2_pCO2" = bquote(italic(p)*CO[2]*','~italic(p)*CO[2]), "pCO2_Control" = bquote(italic(p)*CO[2]*', Control'), "Control_Temp" = "Control, Temp", "Control_pCO2_Temp" = bquote("Control, "*italic(p)*CO[2]*'+Temp'), "Control_pCO2" = bquote("Control, "*italic(p)*CO[2]), "Control_Control" = "Control, Control"), breaks = c(levels(totalactivity_draws$treatment_combined)), limits = c("Control_Control", "skip", "Control_Temp", "skip", "Control_pCO2", "skip", "Control_pCO2_Temp", "skip", "Temp_Control", "skip", "Temp_Temp", "skip", "pCO2_Control", "skip", "pCO2_pCO2", "skip", "pCO2_Temp_Control", "skip", "pCO2_Temp_pCO2_Temp")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.16, 0.11), legend.background = element_blank(), legend.key = element_blank())

#totalactivity_plot

# Write the plot to a pdf
ggsave(filename = "totalactivity.pdf", plot = totalactivity_plot, dpi = 2000, width = 8.03, height = 7)
# Also write out the contrasts
write_tsv(file = "totalactivity_emmeans.txt", x = summary(pairs(emmeans(totalactivity_model, ~ treatment_combined | Cue, re_formula = NULL))))




totalactivity_model_t1 <- brm(data = raw_data, family = skew_normal,
                           bf(Total_Activity_1st ~ Cue * treatment_combined + (1 | Tank)), 
                           prior = c(prior(normal(0, 30), class = Intercept),
                                     prior(normal(0, 10), class = b)),
                           iter = 5000, warmup = 1000, chains = 4, cores = 4,  
                           control = list(adapt_delta = 0.99, max_treedepth = 12))

plot(conditional_effects(totalactivity_model_t1), points = F)

totalactivity_model_t2 <- brm(data = raw_data, family = skew_normal,
                              bf(Total_Activity_2nd ~ Cue * treatment_combined + (1 | Tank)), 
                              prior = c(prior(normal(0, 30), class = Intercept),
                                        prior(normal(0, 10), class = b)),
                              iter = 5000, warmup = 1000, chains = 4, cores = 4,  
                              control = list(adapt_delta = 0.99, max_treedepth = 12))

plot(conditional_effects(totalactivity_model_t2), points = F)

totalactivity_model_t3 <- brm(data = raw_data, family = skew_normal,
                              bf(Total_Activity_3rd ~ Cue * treatment_combined + (1 | Tank)), 
                              prior = c(prior(normal(0, 30), class = Intercept),
                                        prior(normal(0, 10), class = b)),
                              iter = 5000, warmup = 1000, chains = 4, cores = 4,  
                              control = list(adapt_delta = 0.99, max_treedepth = 12))

plot(conditional_effects(totalactivity_model_t3), points = F)



# Look at a model of time near the stormtrooper
timeneartotal_model <- brm(data = raw_data, family = skew_normal,
                           bf(Time_Near_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g)), 
                           prior = c(prior(normal(0, 100), class = Intercept),
                                     prior(normal(0, 100), class = b),
                                     prior(cauchy(0, 10), class = sd)),
                           iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                           control = list(adapt_delta = 0.99, max_treedepth = 12))

summary(timeneartotal_model)

get_variables(totalactivity_model)

pp_check(timeneartotal_model, ndraws = 100) +
  theme_bw(base_size = 20)

plot(conditional_effects(timeneartotal_model), points = F, theme = theme_bw())

plot(emmeans(timeneartotal_model, ~ treatment_combined | Cue, re_formula = NULL))

# Pull model posterior draws with gather emmeans draws
timeneartotal_draws <- gather_emmeans_draws(emmeans(timeneartotal_model, ~ treatment_combined | Cue, re_formula = NULL)) %>% 
  mutate(Cue = factor(Cue, levels = rev(c("Blank", "Food", "Alarm", "Food_Alarm"))))

# Plot posterior 95% credible intervals from the model draws
timeneartotal_plot <- ggplot(timeneartotal_draws, aes(x = .value, y = treatment_combined, colour = Cue, fill = Cue)) + 
  stat_pointinterval(position = position_dodge(width = 1)) +
  scale_fill_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  xlab("Time Near Novel Object (%)") +
  ylab("Rearing Group, Trial Group") +
  scale_y_discrete(labels = c("Temp_Temp" = "Temp, Temp", "Temp_Control" = "Temp, Control", "pCO2_Temp_pCO2_Temp" = bquote(italic(p)*CO[2]*"+Temp,"~italic(p)*CO[2]*"+"*Temp), "pCO2_Temp_Control" = bquote(italic(p)*CO[2]*"+"*'Temp, Control'), "pCO2_pCO2" = bquote(italic(p)*CO[2]*','~italic(p)*CO[2]), "pCO2_Control" = bquote(italic(p)*CO[2]*', Control'), "Control_Temp" = "Control, Temp", "Control_pCO2_Temp" = bquote("Control, "*italic(p)*CO[2]*'+Temp'), "Control_pCO2" = bquote("Control, "*italic(p)*CO[2]), "Control_Control" = "Control, Control"), breaks = c(levels(activityneartotal_draws$treatment_combined)), limits = c("Control_Control", "skip", "Control_Temp", "skip", "Control_pCO2", "skip", "Control_pCO2_Temp", "skip", "Temp_Control", "skip", "Temp_Temp", "skip", "pCO2_Control", "skip", "pCO2_pCO2", "skip", "pCO2_Temp_Control", "skip", "pCO2_Temp_pCO2_Temp")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.14, 0.10), legend.background = element_blank(), legend.key = element_blank())

timeneartotal_plot

# Write the plot to a pdf
ggsave(filename = "timeneartotal.pdf", plot = timeneartotal_plot, dpi = 2000, width = 8.03, height = 7)
# Also write the emmeans contrasts
write_tsv(file = "timeneartotal_emmeans.txt", x = summary(pairs(emmeans(timeneartotal_model, ~ treatment_combined | Cue, re_formula = NULL))))


timethigmotaxistotal_model <- brm(data = raw_data, family = skew_normal,
                           bf(Time_Thigmotaxis_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g),
                              sigma ~ treatment_combined), 
                           prior = c(prior(normal(0, 100), class = Intercept),
                                     prior(normal(0, 100), class = b),
                                     prior(cauchy(0, 10), class = sd)),
                           iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                           control = list(adapt_delta = 0.99, max_treedepth = 12))

summary(timethigmotaxistotal_model)

get_variables(timethigmotaxistotal_model)

pp_check(timethigmotaxistotal_model, ndraws = 100) +
  #xlim(-10, 200) +
  theme_bw(base_size = 20) 


plot(emmeans(timethigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL))

# Pull model posterior draws with gather emmeans draws
timethigmotaxistotal_draws <- gather_emmeans_draws(emmeans(timethigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL)) %>% 
  mutate(Cue = factor(Cue, levels = rev(c("Blank", "Food", "Alarm", "Food_Alarm"))))

# Plot posterior 95% credible intervals from the model draws
timethigmotaxistotal_plot <- ggplot(timethigmotaxistotal_draws, aes(x = .value, y = treatment_combined, colour = Cue, fill = Cue)) + 
  stat_pointinterval(position = position_dodge(width = 1)) +
  scale_fill_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  xlab("Time in Thigmotaxis Zone (%)") +
  ylab("Rearing Group, Trial Group") +
  scale_y_discrete(labels = c("Temp_Temp" = "Temp, Temp", "Temp_Control" = "Temp, Control", "pCO2_Temp_pCO2_Temp" = bquote(italic(p)*CO[2]*"+Temp,"~italic(p)*CO[2]*"+"*Temp), "pCO2_Temp_Control" = bquote(italic(p)*CO[2]*"+"*'Temp, Control'), "pCO2_pCO2" = bquote(italic(p)*CO[2]*','~italic(p)*CO[2]), "pCO2_Control" = bquote(italic(p)*CO[2]*', Control'), "Control_Temp" = "Control, Temp", "Control_pCO2_Temp" = bquote("Control, "*italic(p)*CO[2]*'+Temp'), "Control_pCO2" = bquote("Control, "*italic(p)*CO[2]), "Control_Control" = "Control, Control"), breaks = c(levels(activityneartotal_draws$treatment_combined)), limits = c("Control_Control", "skip", "Control_Temp", "skip", "Control_pCO2", "skip", "Control_pCO2_Temp", "skip", "Temp_Control", "skip", "Temp_Temp", "skip", "pCO2_Control", "skip", "pCO2_pCO2", "skip", "pCO2_Temp_Control", "skip", "pCO2_Temp_pCO2_Temp")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.16, 0.11), legend.background = element_blank(), legend.key = element_blank())

timethigmotaxistotal_plot

# Write the plot to a pdf
ggsave(filename = "timethigmotaxistotal.pdf", plot = timethigmotaxistotal_plot, dpi = 2000, width = 8.03, height = 7)
# Also write the emmeans contrasts
write_tsv(file = "timethigmotaxistotal_emmeans.txt", x = summary(pairs(emmeans(timethigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL))))

activityneartotal_model <- brm(data = raw_data, family = skew_normal,
                                  bf(Activity_Near_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g)), 
                                  prior = c(prior(normal(0, 10), class = Intercept),
                                            prior(normal(0, 10), class = b),
                                            prior(cauchy(0, 5), class = sd)),
                                  iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                                  control = list(adapt_delta = 0.99, max_treedepth = 12))

summary(activityneartotal_model)

get_variables(activityneartotal_model)

plot(activityneartotal_model)

pp_check(activityneartotal_model, ndraws = 100) +
  theme_bw(base_size = 20) 


plot(emmeans(activityneartotal_model, ~ treatment_combined | Cue, re_formula = NULL))

# Pull model posterior draws with gather emmeans draws
activityneartotal_draws <- gather_emmeans_draws(emmeans(activityneartotal_model, ~ treatment_combined | Cue, re_formula = NULL)) %>% 
  mutate(Cue = factor(Cue, levels = rev(c("Blank", "Food", "Alarm", "Food_Alarm"))))

# Plot posterior 95% credible intervals from the model draws
activityneartotal_plot <- ggplot(activityneartotal_draws, aes(x = .value, y = treatment_combined, colour = Cue, fill = Cue)) + 
  stat_pointinterval(position = position_dodge(width = 1)) +
  scale_fill_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  xlab("Total Activity Near Novel Object") +
  ylab("Rearing Group, Trial Group") +
  scale_y_discrete(labels = c("Temp_Temp" = "Temp, Temp", "Temp_Control" = "Temp, Control", "pCO2_Temp_pCO2_Temp" = bquote(italic(p)*CO[2]*"+Temp,"~italic(p)*CO[2]*"+"*Temp), "pCO2_Temp_Control" = bquote(italic(p)*CO[2]*"+"*'Temp, Control'), "pCO2_pCO2" = bquote(italic(p)*CO[2]*','~italic(p)*CO[2]), "pCO2_Control" = bquote(italic(p)*CO[2]*', Control'), "Control_Temp" = "Control, Temp", "Control_pCO2_Temp" = bquote("Control, "*italic(p)*CO[2]*'+Temp'), "Control_pCO2" = bquote("Control, "*italic(p)*CO[2]), "Control_Control" = "Control, Control"), breaks = c(levels(activityneartotal_draws$treatment_combined)), limits = c("Control_Control", "skip", "Control_Temp", "skip", "Control_pCO2", "skip", "Control_pCO2_Temp", "skip", "Temp_Control", "skip", "Temp_Temp", "skip", "pCO2_Control", "skip", "pCO2_pCO2", "skip", "pCO2_Temp_Control", "skip", "pCO2_Temp_pCO2_Temp")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.14, 0.1), legend.background = element_blank(), legend.key = element_blank())

activityneartotal_plot

# Write the plot to a pdf
ggsave(filename = "activityneartotal.pdf", plot = activityneartotal_plot, dpi = 2000, width = 8.03, height = 7)
# Also write the emmeans contrasts
write_tsv(file = "activityneartotal_emmeans.txt", x = summary(pairs(emmeans(activityneartotal_model, ~ treatment_combined | Cue, re_formula = NULL))))



activitythigmotaxistotal_model <- brm(data = raw_data, family = skew_normal,
                               bf(Activity_Thigmotaxis_Total ~ (treatment_combined | Cue) + (1 | Tank / TreatTankID) + (Length_mm * Weight_g),
                                  sigma ~ treatment_combined), 
                               prior = c(prior(normal(0, 20), class = Intercept),
                                         prior(normal(0, 20), class = b),
                                         prior(cauchy(0, 5), class = sd)),
                               iter = 20000, warmup = 5000, chains = 4, cores = 4,  
                               control = list(adapt_delta = 0.99, max_treedepth = 12))

summary(activitythigmotaxistotal_model)

get_variables(activitythigmotaxistotal_model)

pp_check(activitythigmotaxistotal_model, ndraws = 100) +
  theme_bw(base_size = 20) 


plot(emmeans(activitythigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL))

# Pull model posterior draws with gather emmeans draws
activitythigmotaxistotal_draws <- gather_emmeans_draws(emmeans(activitythigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL)) %>% 
  mutate(Cue = factor(Cue, levels = rev(c("Blank", "Food", "Alarm", "Food_Alarm"))))

# Plot posterior 95% credible intervals from the model draws
activitythigmotaxistotal_plot <- ggplot(activitythigmotaxistotal_draws, aes(x = .value, y = treatment_combined, colour = Cue, fill = Cue)) + 
  stat_pointinterval(position = position_dodge(width = 1)) +
  scale_fill_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  scale_colour_fish_d(option = "Epinephelus_lanceolatus", labels = c("Alarm" = "Alarm", "Blank" = "Blank", "Food" = "Food", "Food_Alarm" = "Food & Alarm"), limits = c("Blank", "Food", "Alarm", "Food_Alarm"), direction = -1) +
  xlab("Total Activity in Thigmotaxis Zone") +
  ylab("Rearing Group, Trial Group") +
  scale_y_discrete(labels = c("Temp_Temp" = "Temp, Temp", "Temp_Control" = "Temp, Control", "pCO2_Temp_pCO2_Temp" = bquote(italic(p)*CO[2]*"+Temp,"~italic(p)*CO[2]*"+"*Temp), "pCO2_Temp_Control" = bquote(italic(p)*CO[2]*"+"*'Temp, Control'), "pCO2_pCO2" = bquote(italic(p)*CO[2]*','~italic(p)*CO[2]), "pCO2_Control" = bquote(italic(p)*CO[2]*', Control'), "Control_Temp" = "Control, Temp", "Control_pCO2_Temp" = bquote("Control, "*italic(p)*CO[2]*'+Temp'), "Control_pCO2" = bquote("Control, "*italic(p)*CO[2]), "Control_Control" = "Control, Control"), breaks = c(levels(activityneartotal_draws$treatment_combined)), limits = c("Control_Control", "skip", "Control_Temp", "skip", "Control_pCO2", "skip", "Control_pCO2_Temp", "skip", "Temp_Control", "skip", "Temp_Temp", "skip", "pCO2_Control", "skip", "pCO2_pCO2", "skip", "pCO2_Temp_Control", "skip", "pCO2_Temp_pCO2_Temp")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.14, 0.52), legend.background = element_blank(), legend.key = element_blank())

activitythigmotaxistotal_plot

# Write the plot to a pdf
ggsave(filename = "activitythigmotaxistotal.pdf", plot = activitythigmotaxistotal_plot, dpi = 2000, width = 8.03, height = 7)
# Also write the emmeans contrasts
write_tsv(file = "activitythigmotaxistotal_emmeans.txt", x = summary(pairs(emmeans(activitythigmotaxistotal_model, ~ treatment_combined | Cue, re_formula = NULL))))
