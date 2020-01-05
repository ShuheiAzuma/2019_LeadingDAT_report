rm(list = ls())

# load library ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(KFAS)
library(knitr)

# load data ---------------------------------------------------------------

path_csv <- file.path(".", "002__input", "data_1.csv")

x <- read_csv(path_csv) %>% 
  glimpse()

colnames(x) <- c("year_month",
                 "total",
                 "china",
                 "hongkong",
                 "taipei",
                 "korea",
                 "filipinas",
                 "thai",
                 "singapole",
                 "malaysia",
                 "indonesia",
                 "india",
                 "vietnam",
                 "UK",
                 "germany",
                 "france",
                 "rosia",
                 "italy",
                 "USA",
                 "canada",
                 "australia")

x$year_month <- ymd(x$year_month)

# x <- x %>% mutate_if(is.numeric, log)

germany_ts <- x$germany

# AIC function ----

AIC_ammend_nonGaussian <- function(x, degree){

  return(AIC_ <- 2 * x$optim.out$value + 2 * degree)

  }

# For NA ------------------------------------------------------------------

set.seed(123)

model_1_s <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA) +
                   SSMseasonal(12, q),
                 H = NA)
  
  if(is.na(q)) {
    fit_deg <- 3
  } else {
    fit_deg <- 2
  }
  
  fit <- fitSSM(mod, numeric(fit_deg), method = "SANN")
  kfs <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

result_model_1_s <- model_1_s(germany_ts, NA)

res_obs <- rstandard(result_model_1_s$kfs, "pearson")
res_state <- rstandard(result_model_1_s$kfs, "state")[, 1]

germany_res <- x %>% 
  select(year_month, germany) %>% 
  bind_cols(res_obs = res_obs,
            res_state = res_state)

res_obs_out <- germany_res %>%
  filter(abs(res_obs) > 2) %>% 
  select(-res_state)

res_state_out <- germany_res %>% 
  filter(abs(res_state) > 2) %>% 
  select(-res_obs)

# write_csv(res_obs_out, "./003__output/res_obs_out.csv")
# write_csv(res_obs_out, "./003__output/res_state_out.csv")

remove_index <- abs(germany_res$res_obs) > 1.96

germany_ts_rm <- germany_ts
germany_ts_rm[remove_index] <- NA


# model 2 Poisson seasonal * smoothed trend * ----

model_1_s_pois <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA) +
                   SSMseasonal(12, NA),
                 distribution = "poisson",
                 u = 1)
  
  fit <- fitSSM(mod, c(-15, -10), nsim = 1000, method = "BFGS")
  kfs <- KFS(fit$model,  nsim = 10000)
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# Culculate --------------------------------------------------------------

# germany_ts <- x$germany %>% magrittr::multiply_by(10^-3)

result_pois <- model_1_s_pois(germany_ts_rm, NA)

pre_Interval_Pois <- predict(result_pois$fit$model, interval = "prediction", level = 0.95, nsim = 1000) %>% 
  as_tibble()

AIC_pois_n30_n20 <- AIC_ammend_nonGaussian(result_pois$fit, degree = 1 + 12)

# plot --------------------------------------------------------------------

pred <- tibble(pred_level = result_pois$kfs$alphahat[, "level"] %>% exp() %>% as.numeric()) %>% 
  mutate(seasonal = result_pois$kfs$muhat - pred_level,
         pred_mu = result_pois$kfs$muhat) %>% 
  bind_cols(year_month = x$year_month,
            germany = x$germany,
            interval = pre_Interval_Pois)



theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12,
                                           angle = 90)
            )
)


plot_pred <- ggplot(pred,
                    aes(x = year_month, y = seasonal)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.6) +
  
  geom_line(aes(y = upr),
            lty = "dashed",
            alpha = 0.2) +
  
  geom_line(aes(y = lwr),
            lty = "dashed",
            alpha = 0.2) +
  
  # Germany
  
  geom_point(aes(y = germany),
             alpha = 0.6) +

  geom_line(aes(y = germany),
            alpha = 0.6) +

  # only Calendar effect
  geom_line(alpha = 0.6,
            colour = "tomato") +
  
  # Prediction of level 
  geom_line(aes(y = pred_level),
            colour = "blue",
            alpha = 0.8) +
  
  # geom_line(aes(y = pred)) +
  
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  labs(x = "Year/Month",
       y = "Number of foreign visitors") +
  NULL


plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_pois_model_sep_effect_level.png",
       width = 6,
       height = 4)



# AIC_curve ---------------------------------------------------------------

df_AIC_curve <- tibble(name = factor(c("Initial_c(-30, -20)",
                                       "Initial_c(-15, -10)",
                                       "Initial_c(  0,   0)",
                                       "Initial_c( 15,  10)"),
                                     levels = c("Initial_c(-30, -20)",
                                                "Initial_c(-15, -10)",
                                                "Initial_c(  0,   0)",
                                                "Initial_c( 15,  10)")),
                       value = c(AIC_pois_n30_n20,
                                 AIC_pois_n15_n10,
                                 AIC_pois_0_0,
                                 AIC_pois_p15_p10)
                       )


theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 20,
                                            face = "bold"),
                  axis.text = element_text(size = 20,
                                             angle = 90)
            )
)


g <- ggplot(df_AIC_curve,
       aes(x = name, y = value)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = value - 200,
                label = value %>% round()),
            color = "white",
            fontface = "bold",
            size = 12) +
  scale_y_continuous(breaks = seq(2500, 4500, 500)) +
  coord_cartesian(ylim = c(2500, 4500)) +
  labs(x = "Initial Value",
       y = "AIC") +
  NULL

ggsave(g, filename = "./003__output/AIC_comparison_poisson.png",
       width = 8,
       height = 6)
