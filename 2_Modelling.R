### PA Priorization
### Tjaden-McClement et al., 2024
### Modelling

library(tidyverse)
library(sjPlot)
library(AICcmodavg)

theme_set(theme_classic())

# load flattened yearly wdpa dataset made in 1_Matching
load("G:/PA_Prioritization/Key_R_files/pa_by_year_flat.RData")

##### LOTW #####

# Load in matched LOTW data (exact on country and biome)
load("Key_R_files/Matching/LOTW/lotw_match_data_cb_cal.RData")
head(lotw_match_data_cb_cal)

panel_fun <- function(yearly_input, datID_input){
  datID_input %>% 
    mutate(year = max(yearly_input$year), 
           protected = ifelse(datID %in% yearly_input$datID, 1, 0))
}

lotw_match_datIDs <- select(lotw_match_data_cb_cal, "datID", "lotw")

lotw_panel <- lapply(pa_by_year_flat, FUN = panel_fun,
                     datID_input = lotw_match_datIDs)
lotw_panel <- do.call(rbind, lotw_panel)

lotw_panel <- arrange(lotw_panel, datID, year)

# add "trt_est" column that's 1 if after lotw established (>=2002) and 0 if before
lotw_panel$trt_est <- ifelse(lotw_panel$year >= 2002, 1, 0)
head(lotw_panel)

lotw_panel <- lotw_panel %>% 
  mutate(lotw = as.factor(lotw),
         trt_est = as.factor(trt_est),
         protected = as.factor(protected))
str(lotw_panel)

save(lotw_panel, file = "Key_R_files/Modeling/lotw_panel.RData")

lotw_agg <- lotw_panel %>% 
  mutate(protected = if_else(protected == 0, 0, 1)) %>% 
  group_by(year, lotw) %>% 
  summarise(prop_protected = mean(protected)) %>% 
  mutate(year_centered = year - 2002,
         trt_est = if_else(year_centered > 0, 1, 0))
head(lotw_agg)

ggplot(lotw_agg, aes(x = year,
                     y = prop_protected, 
                     colour = lotw, 
                     group = lotw)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 2002, colour = "red",
             linetype = "dashed")

### difference in differences model with immediate and trend changes
lotw_agg_did_model <- lm(prop_protected ~ trt_est*lotw*year_centered,
                         data = lotw_agg)
summary(lotw_agg_did_model)

lotw_agg <- tibble(lotw_agg, 
                   predicted_protection = as.vector(predict(lotw_agg_did_model)))

ggplot(lotw_agg, aes(x = year,
                     y = prop_protected, 
                     colour = lotw, 
                     group = lotw)) +
  geom_vline(xintercept = 2002, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#9769C1"),
                     name = "Treatment", 
                     labels = c("Control", "Last of the Wild")) +
  geom_line(aes(x = year,
                y = predicted_protection),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Null model that doesn't split the data into before and after LOTW established

lotw_agg_null <- lm(prop_protected ~ lotw + year,
                    data = lotw_agg)
lotw_null_predict <- tibble(lotw_agg, 
                            predicted_protection_null = as.vector(predict(lotw_agg_null)))

ggplot(lotw_agg, aes(x = year,
                     y = prop_protected, 
                     colour = lotw, 
                     group = lotw)) +
  geom_vline(xintercept = 2002, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#9769C1"),
                     name = "Treatment", 
                     labels = c("Control", "Last of the Wild")) +
  geom_line(data = lotw_null_predict,
            aes(x = year,
                y = predicted_protection_null),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Model with trend effects only (no immediate effect of treatment)

lotw_agg_trend <- lm(prop_protected ~ lotw + year_centered + trt_est:year_centered +
                       lotw:year_centered + trt_est:lotw:year_centered,
                     data = lotw_agg)
summary(lotw_agg_trend)

lotw_trend_predict <- tibble(lotw_agg, 
                             predicted_protection_trend = as.vector(predict(lotw_agg_trend)))

ggplot(lotw_agg, aes(x = year,
                     y = prop_protected, 
                     colour = lotw, 
                     group = lotw)) +
  geom_vline(xintercept = 2002, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#9769C1"),
                     name = "Treatment", 
                     labels = c("Control", "Last of the Wild")) +
  geom_line(data = lotw_trend_predict,
            aes(x = year,
                y = predicted_protection_trend),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Comparing models with AIC
lotw_cand_models <- list(lotw_agg_did_model, lotw_agg_trend, lotw_agg_null)
model_names <- c("full", "trend_only", "null")

aictab(lotw_cand_models, model_names)
# full model is best by far


##### Hotspots #####

# Load in matched Hotspots data (exact on country and biome)
load("Key_R_files/Matching/Hotspots/hotspot_match_data_cb_cal.RData")
head(hotspot_match_data_cb_cal)

hotspots_match_datIDs <- select(hotspot_match_data_cb_cal, "datID", "hotspot")

hotspots_panel <- lapply(pa_by_year_flat, FUN = panel_fun,
                         datID_input = hotspots_match_datIDs)
hotspots_panel <- do.call(rbind, hotspots_panel)

hotspots_panel <- arrange(hotspots_panel, datID, year)

# add "trt_est" column that's 1 if after hotspots established (>=2000) and 0 if before
hotspots_panel$trt_est <- ifelse(hotspots_panel$year >= 2000, 1, 0)
head(hotspots_panel)

hotspots_panel <- hotspots_panel %>% 
  mutate(hotspot = as.factor(hotspot),
         trt_est = as.factor(trt_est),
         protected = as.factor(protected))
str(hotspots_panel)

#save(hotspots_panel, file = "Key_R_files/Modeling/hotspots_panel.RData")
#load("Key_R_files/Modeling/hotspots_panel.RData")

hotspots_agg <- hotspots_panel %>% 
  mutate(protected = if_else(protected == 0, 0, 1)) %>% 
  group_by(year, hotspot) %>% 
  summarise(prop_protected = mean(protected)) %>% 
  mutate(year_centered = year - 2000,
         trt_est = if_else(year_centered > 0, 1, 0))
head(hotspots_agg)

# overall trend
ggplot(hotspots_agg, aes(x = year,
                         y = prop_protected, 
                         colour = hotspot, 
                         group = hotspot)) +
  geom_vline(xintercept = 2000, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#1F77B4FF"),
                     name = "Treatment", 
                     labels = c("Control", "Hotspots")) +
  labs(x = "Year", y = "Proportion area protected")

# difference in differences model with immediate and trend changes
hotspot_agg_did_model <- lm(prop_protected ~ trt_est*hotspot*year_centered,
                            data = hotspots_agg)
summary(hotspot_agg_did_model)

hotspots_agg <- tibble(hotspots_agg, 
                       predicted_protection = as.vector(predict(hotspot_agg_did_model)))

ggplot(hotspots_agg, aes(x = year,
                         y = prop_protected, 
                         colour = hotspot, 
                         group = hotspot)) +
  geom_vline(xintercept = 2000, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#1F77B4FF"),
                     name = "Treatment", 
                     labels = c("Control", "Hotspots")) +
  geom_line(aes(x = year,
                y = predicted_protection),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Null model that doesn't split the data into before and after Hotspots established
hotspots_agg_null <- lm(prop_protected ~ hotspot + year,
                        data = hotspot_agg)
hotspots_null_predict <- tibble(hotspot_agg, 
                                predicted_protection_null = as.vector(predict(hotspots_agg_null)))

ggplot(hotspot_agg, aes(x = year,
                        y = prop_protected, 
                        colour = hotspot, 
                        group = hotspot)) +
  geom_vline(xintercept = 2000, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#1F77B4FF"),
                     name = "Treatment", 
                     labels = c("Control", "Hotspots")) +
  geom_line(data = hotspots_null_predict,
            aes(x = year,
                y = predicted_protection_null),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Model with trend effects only (no immediate effect of treatment)
hotspots_agg_trend <- lm(prop_protected ~ hotspot + year_centered + trt_est:year_centered +
                           hotspot:year_centered + trt_est:hotspot:year_centered,
                         data = hotspot_agg)
summary(hotspots_agg_trend)

hotspots_trend_predict <- tibble(hotspot_agg, 
                                 predicted_protection_trend = as.vector(predict(hotspots_agg_trend)))

ggplot(hotspot_agg, aes(x = year,
                        y = prop_protected, 
                        colour = hotspot, 
                        group = hotspot)) +
  geom_vline(xintercept = 2000, colour = "#818181",
             linetype = "dashed", size = 0.7) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = c("black", "#1F77B4FF"),
                     name = "Treatment", 
                     labels = c("Control", "Hotspots")) +
  geom_line(data = hotspots_trend_predict,
            aes(x = year,
                y = predicted_protection_trend),
            colour = "red", size = 0.7) +
  labs(x = "Year", y = "Proportion area protected")

### Compare models with AIC
hotspots_models <- list(hotspot_agg_did_model, hotspots_agg_trend, hotspots_agg_null)

hotspots_aic_table <- as_tibble(aictab(hotspots_models, model_names))

