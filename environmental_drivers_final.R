#####################################
## Environmental drivers modelling ##
#####################################
library(tidyverse)
library(ggeffects)
library(patchwork)


## Model drivers

dat <- 
  read_csv("2026-03-16_extractEnv_data.csv") %>% 
  mutate(stranded = as.integer(stranded),
         species = as_factor(species),
         state = as_factor(state),
         match_id = as_factor(event_id)) %>% 
  filter(!is.na(sst_anom),
         !is.na(wind_stress_onshore),
         !is.na(current_onshore),
         !is.na(species),
         !is.na(state),
         !is.na(match_id))


## Distribution plots

p1 <-
  dat %>% mutate(ss = case_when(stranded %in% 0 ~ "Control", TRUE ~ "Stranded")) %>% 
  ggplot(aes(x = sst_anom, y = after_stat(scaled), fill = ss, group = ss)) +
  geom_hline(yintercept = 0, lwd = 0.2) +
  geom_density(alpha = 0.3, col = NA) +
  geom_point(aes(y = ifelse(ss == "Control", -0.015, -0.04), color = ss), 
             alpha = 0.5, pch = "|", cex = 2, show.legend = F) +
  scale_y_continuous(expand = expansion(mult = c(0.015, 0.02))) +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.5) +
  theme_classic() +
  labs(x = "Sea Surface Temperature Anomaly (˚C)", y = "Scaled density") +
  scale_fill_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  scale_color_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  theme(legend.position = "inside", legend.position.inside = c(0.15, 0.9))

p2 <-
  dat %>% mutate(ss = case_when(stranded %in% 0 ~ "Control", TRUE ~ "Stranded")) %>% 
  ggplot(aes(x = wind_stress_onshore, y = after_stat(scaled), fill = factor(stranded), group = stranded)) +
  geom_hline(yintercept = 0, lwd = 0.2) +
  geom_density(alpha = 0.3, col = NA) +
  geom_point(aes(y = ifelse(stranded == 0, -0.015, -0.04), 
                 color = factor(stranded)), alpha = 0.5, pch = "|", cex = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.015, 0.02))) +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.5) +
  theme_classic() +
  labs(x = expression("Onshore Wind Stress (N "*m^-2*")"), y = NULL) +
  scale_fill_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  scale_color_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  scale_x_continuous(breaks = seq(-16, 16, by = 4)) +
  theme(legend.position = "none")

p3 <-
  dat %>% mutate(ss = case_when(stranded %in% 0 ~ "Control", TRUE ~ "Stranded")) %>% 
  ggplot(aes(x = current_onshore, y = after_stat(scaled), fill = factor(stranded), group = stranded)) +
  geom_hline(yintercept = 0, lwd = 0.2) +
  geom_density(alpha = 0.3, col = NA) +
  geom_point(aes(y = ifelse(stranded == 0, -0.015, -0.04), 
                 color = factor(stranded)), alpha = 0.5, pch = "|", cex = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.015, 0.02))) +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.5) +
  theme_classic() +
  labs(x = expression("Onshore Current (m "*s^-1*")"), y = NULL) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  scale_fill_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  scale_color_brewer(name = NULL, palette = "Set1", labels = c("Control", "Stranded"), direction = -1) +
  theme(legend.position = "none")

dist <- p1 + p2 + p3

ggsave("Figures/extractEnv_distribution_plots.png", dist, width = 14, height = 5)



## GLM modelling

glm_model <-
  glm(stranded ~ sst_anom + wind_stress_onshore + current_onshore,
      family = binomial(link = "logit"),
      data = dat)

p1 <-
  predict_response(glm_model, terms = c("sst_anom [all]")) %>%
  plot(show_data = T, show_title = F, color = "coral", dot_alpha = 0.2) +
  labs(x = "Sea Surface Temperature Anomaly (˚C)",
       y = "Probability of stranding") +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  geom_vline(xintercept = 0, lwd = 0.3, lty = 2, color = "black") +
  theme_classic()

p2 <-
  predict_response(glm_model, terms = c("wind_stress_onshore [all]")) %>%
  plot(show_data = T, show_title = F, color = "steelblue3") +
  labs(x = expression("Onshore Wind stress (N "*m^-2*")"),
       y = NULL) +
  scale_x_continuous(breaks = seq(-15, 15, by = 4)) +
  geom_vline(xintercept = 0, lwd = 0.3, lty = 2, color = "black") +
  theme_classic()

p3 <-
  predict_response(glm_model, terms = c("current_onshore [all]")) %>%
  plot(show_data = T, show_title = F, color = "gold3") +
  labs(x = expression("Onshore Current (m "*s^-1*")"),
       y = NULL) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  geom_vline(xintercept = 0, lwd = 0.3, lty = 2, color = "black") +
  theme_classic()

glm_out <- p1 + p2 + p3

ggsave("Figures/extractEnv_glm_output.png", glm_out, width = 14, height = 5)


#Statistical Reporting
nrow(dat)
table(dat$stranded)
summary(glm_model)

