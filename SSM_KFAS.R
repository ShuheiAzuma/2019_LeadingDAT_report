
rm(list = ls())

# load library ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(KFAS)


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

germany_ts <- x$germany %>% ts

# model, fit, KFS ----

mod <- SSModel(germany_ts ~ SSMtrend(1, Q = NA), H = NA)
fit <- fitSSM(mod, numeric(2), method = "SANN")
kfs <- KFS(fit$model)

# model 1 * local level * ----

model_1 <- function(x){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA), H = NA)
  fit <- fitSSM(mod, numeric(2), method = "SANN")
  kfs <- KFS(fit$model)

  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# model 2 * smoothed trend * ----

model_2 <- function(x){
  
  mod <- SSModel(x ~ SSMtrend(2, Q = list(0, NA)), H = NA)
  fit <- fitSSM(mod, numeric(2), method = "SANN")
  kfs <- KFS(fit$model)
  
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# model 3 * local linear trend * ----

model_3 <- function(x){
  
  mod <- SSModel(x ~ SSMtrend(2, Q = list(NA, NA)), H = NA)
  fit <- fitSSM(mod, numeric(3), method = "SANN")
  kfs <- KFS(fit$model)
  
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# model 1 seasonal * local level * ----

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
  kfs <- KFS(fit$model)
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# model 2  seasonal * smoothed trend * ----

model_2_s <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(2, Q = list(0, NA)) +
                   SSMseasonal(12, q),
                 H = NA)
  
  if(is.na(q)) {
    fit_deg <- 3
  } else {
    fit_deg <- 2
  }
  
  fit <- fitSSM(mod, numeric(fit_deg), method = "SANN")
  kfs <- KFS(fit$model)
  
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

# model 3 seasonal * local linear trend * ----

model_3_s <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(2, Q = list(NA, NA)) +
                   SSMseasonal(12, q),
                 H = NA)
  
  if(is.na(q)) {
    fit_deg <- 4
  } else {
    fit_deg <- 3
  }
  
  fit <- fitSSM(mod, numeric(fit_deg), method = "SANN")
  kfs <- KFS(fit$model)
  
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}



# calculate model , fit, kfs ----

result_model_1 <- model_1(germany_ts)
result_model_2 <- model_2(germany_ts)
result_model_3 <- model_3(germany_ts)

result_model_1_s_fixed <- model_1_s(germany_ts, 0)
result_model_1_s_fluct <- model_1_s(germany_ts, NA)
result_model_2_s_fixed <- model_2_s(germany_ts, 0)
result_model_2_s_fluct <- model_2_s(germany_ts, NA)
result_model_3_s_fixed <- model_3_s(germany_ts, 0)
result_model_3_s_fluct <- model_3_s(germany_ts, NA)

# log likelyhood ----

AIC_ammend <- function(x, degree){
  loglik <- x$logLik - sum(x$Finf > 0) * log(2*pi) / 2
  return(AIC_ <- - 2 * loglik + 2 * degree)
}

df_AIC <- tibble(
  
  model_name = c(
    "AIC_1",
    "AIC_2",
    "AIC_3",
    "AIC_1_s_fixed",
    "AIC_1_s_fluct",
    "AIC_2_s_fixed",
    "AIC_2_s_fluct",
    "AIC_3_s_fixed",
    "AIC_3_s_fluct"
  ),
  
  AIC_culc = c(AIC_ammend(result_model_1$kfs, degree = 2 + 1),
               AIC_ammend(result_model_2$kfs, degree = 2 + 2),
               AIC_ammend(result_model_3$kfs, degree = 3 + 2),
               AIC_ammend(result_model_1_s_fixed$kfs, degree = 2 + 1 + 11),
               AIC_ammend(result_model_1_s_fluct$kfs, degree = 3 + 1 + 11),
               AIC_ammend(result_model_2_s_fixed$kfs, degree = 2 + 2 + 11),
               AIC_ammend(result_model_2_s_fluct$kfs, degree = 3 + 2 + 11),
               AIC_ammend(result_model_3_s_fixed$kfs, degree = 3 + 2 + 11), 
               AIC_ammend(result_model_3_s_fluct$kfs, degree = 4 + 2 + 11)
  )
  
)

plot_aic <- ggplot(df_AIC,
       aes(x = model_name, y = AIC_culc)) +
  geom_col(aes(x = model_name, y = AIC_culc),
           alpha = 0.9,
           fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(aes(y = AIC_culc - 100,
                label = AIC_culc %>% round()),
            size = 5) +
  geom_hline(aes(yintercept = min(AIC_culc)),
             colour = "tomato",
             lty = "dashed") +
  coord_cartesian(ylim = c(3000, 3750)) +
  NULL

ggsave(plot_aic, filename = "./003__output/plot_aic.png",
       width = 6,
       height = 3)

# predict ----
alphahatconf <- predict(result_model_1_s_fluct$fit$model, interval = "prediction", level = 0.95) %>% as_tibble()

x <- x %>% bind_cols(alphahatconf)

# plot ----

theme_set(theme_minimal() +
          theme(axis.title = element_text(size = 12),
                axis.text = element_text(size = 12,
                                         angle = 90)
                )
          )

plot_pred <- ggplot(x,
       aes(x = year_month, y = germany)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.4) +
  geom_line(aes(y = fit),
            colour = "lightblue",
            alpha = 0.8) +
  
  geom_line(alpha = 0.6) +
  geom_point(size = 2,
             alpha = 0.6) +
  
  expand_limits(y = 0) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  labs(x = "Year/Month",
       y = "Number of foreign visitors from Germany") +
  NULL


ggsave(plot_pred, filename = "./003__output/plot_pred.png",
       width = 6,
       height = 4)


# zansa bunseki  ----------------------------------------------------------

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


plot_germany <- ggplot(germany_res %>% na.omit()) +
  geom_line(aes(x = year_month,
                y = germany/1000),
            alpha = 0.6) +
  expand_limits(y = 0) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  # scale_y_continuous(breaks = seq(-10, 10, 2)) + 
  # coord_cartesian(ylim = c(-8, 8)) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 8)) +
  labs(x = "year/month",
       y = "Number of foreign visitors\n from Geramany / 1000") +
  NULL

  ggsave(plot_germany, filename = "./003__output/plot_germany.png",
         width = 8, height = 2)
  
plot_res <- function(x = germany_res, y_string) {
  
  ggplot(germany_res %>% na.omit()) +
    geom_line(aes_string(x = "year_month",
                  y = y_string)) +
    expand_limits(y = 0) +
    scale_x_date(date_labels = "%y/%m",
                 date_breaks = "1 year") +
    geom_hline(yintercept = -1.96, lty = "dashed") +
    geom_hline(yintercept = 1.96, lty = "dashed") +
    scale_y_continuous(breaks = seq(-10, 10, 2)) + 
    coord_cartesian(ylim = c(-8, 8)) +
    theme(axis.text.x = element_text(angle = 90)) +
    NULL

}

theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  axis.text.y = element_text(size = 8)
            )
)


plot_res_obs <- plot_res(y_string = "res_obs")
plot_res_state <- plot_res(y_string = "res_state")

ggsave(plot_res_obs, filename = "./003__output/plot_res_obs.png",
       width = 8, height = 2)

ggsave(plot_res_state, filename = "./003__output/plot_res_state.png",
       width = 8, height = 2)


res_obs_out <- germany_res %>%
  filter(abs(res_obs) > 2) %>% 
  select(-res_state)

res_state_out <- germany_res %>% 
  filter(abs(res_state) > 2) %>% 
  select(-res_obs)


# ar model ----------------------------------------------------------------------


model_1_s_ar <- function(x){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA) +
                   SSMseasonal(12, NA) +
                   SSMarima(ar = 0,  Q = NA),
                 H = NA)

  updatefn <- function(pars, model){
    SSModel(x ~ SSMtrend(1, Q = exp(pars[1])) +
              SSMseasonal(12, exp(pars[2])) +
              SSMarima(ar = artransform(pars[3]), Q = exp(pars[4])),
            H = exp(pars[5]))
  }
  
  fit <- fitSSM(mod, numeric(5), updatefn, method = "SANN")
  kfs <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

germany_ts <- germany_ts %>% na.omit()

result_model_1_s_ar <- model_1_s_ar(germany_ts)

AIC_ar <- AIC_ammend(result_model_1_s_ar$kfs, degree = 5 + 1 + 11)


alphahatconf <- predict(result_model_1_s_ar$fit$model, interval = "prediction", level = 0.95) %>% as_tibble()

x_add <- x %>%
  na.omit() %>%
  bind_cols(alphahatconf)

# plot ----

theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12,
                                           angle = 90)
            )
)

plot_pred <- ggplot(x_add,
                    aes(x = year_month, y = germany)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.4) +
  geom_line(aes(y = fit),
            colour = "lightblue",
            size = 0.5,
            alpha = 0.8) +
  
  geom_line(alpha = 0.6) +
  geom_point(size = 2,
             alpha = 0.6) +
  
  expand_limits(y = 0) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  labs(x = "Year/Month",
       y = "Number of foreign visitors from Germany") +
  NULL


ggsave(plot_pred, filename = "./003__output/plot_pred_ar_model.png",
       width = 6,
       height = 4)


# NA ka ----

remove_index <- abs(germany_res$res_obs) > 1.96

germany_ts_rm <- germany_ts
germany_ts_rm[remove_index] <- NA

germany_ts_rm <- germany_ts_rm[c(1:192)]

result_model_1_s_fluct_na <- model_1_s(germany_ts_rm, NA)

AIC_na <- AIC_ammend(result_model_1_s_fluct_na$kfs, degree = 3 + 1 + 11)


alphahatconf <- predict(result_model_1_s_fluct_na$fit$model,
                        interval = "prediction",
                        level = 0.95) %>% as_tibble()

level <- tibble(pred_level = result_model_1_s_fluct_na$kfs$alphahat[, "level"] %>% as.numeric())

x_add <- x %>%
  na.omit() %>% 
  select(year_month, germany) %>% 
  # mutate_at(vars(germany), ~ log(.)) %>%
  bind_cols(alphahatconf,
            level)

plot_pred <- ggplot(x_add,
                    aes(x = year_month, y = germany)) +

  geom_line(aes(y = pred_level),
            colour = "blue",
            size = 0.1,
            alpha = 0.8) +
  
  geom_line(alpha = 0.6) +
  
  geom_line(aes(y = upr),
            lty = "dashed",
            alpha = 0.6) +
  
  geom_line(aes(y = lwr),
            lty = "dashed",
            alpha = 0.6) +
  
  geom_point(size = 2,
             alpha = 0.6) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.8) +
  
  # expand_limits(y = 0) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
               # limit=c(as.Date("2010-01-01"),as.Date("2015-01-01"))) +

  labs(x = "Year/Month",
       y = "Number of foreign visitors from Germany") +
  
  NULL

plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_na_model_with_pred.png",
       width = 6,
       height = 4)



# kansho hennsu dounyu ----


germany_ts <- germany_ts %>% na.omit()
shiftlevel <- (1:NROW(germany_ts) >= which(x$year_month == as.Date("2011-01-01")))
shiftSlope <- 1:NROW(germany_ts) - which(x$year_month == as.Date("2011-01-01"))
shiftSlope <- ifelse(shiftSlope < 0, 0, shiftSlope)

model_1_s_interferance <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA) +
                   SSMseasonal(12, q) +
                   shiftSlope,
                 H = NA)
  
  fit <- fitSSM(mod, numeric(3), method = "SANN")
  kfs <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
}

result_model_1_s_interferance <- model_1_s_interferance(germany_ts, NA)


AIC_interf <- AIC_ammend(result_model_1_s_interferance$kfs, degree = 3 + 1 + 1 + 11)


alphahatconf <- predict(result_model_1_s_interferance$fit$model,
                        interval = "prediction",
                        level = 0.95) %>% as_tibble()

level <- tibble(pred_level = result_model_1_s_interferance$kfs$alphahat[, "level"] %>% as.numeric())

x_add <- x %>%
  na.omit() %>% 
  select(year_month, germany) %>% 
  bind_cols(alphahatconf,
            level)

plot_pred <- ggplot(x_add,
                    aes(x = year_month, y = germany)) +
  
  geom_line(aes(y = pred_level),
            colour = "blue",
            size = 0.1,
            alpha = 0.8) +
  
  geom_line(alpha = 0.6) +
  
  geom_line(aes(y = upr),
            lty = "dashed",
            alpha = 0.6) +
  
  geom_line(aes(y = lwr),
            lty = "dashed",
            alpha = 0.6) +
  
  geom_point(size = 2,
             alpha = 0.6) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.8) +
  
  # expand_limits(y = 0) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  labs(x = "Year/Month",
       y = "Number of foreign visitors from Germany") +

  NULL

plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_na_model_with_pred.png",
       width = 6,
       height = 4)


# AIC 

df_AIC_2 <- tibble(model_name = c("3-1_NA_model", "3-2_AR_model", "3-3_Interferance_slope_model"),
                   AIC_culc = c(AIC_na, AIC_ar, AIC_interf))

# write_csv(df_AIC_2, path = "./003__output/df_AIC_2.csv")

plot_aic <- ggplot(df_AIC_2,
                   aes(x = model_name, y = AIC_culc)) +
  geom_col(aes(x = model_name, y = AIC_culc),
           alpha = 0.9,
           fill = "lightblue") +
  
  geom_text(aes(y = AIC_culc - 100,
                label = AIC_culc %>% round()),
            size = 5) +
  
  geom_hline(aes(yintercept = min(AIC_culc)),
             colour = "tomato",
             lty = "dashed") +
  
  coord_cartesian(ylim = c(2500, 3500)) +
  
  theme(axis.text.x = element_text(angle = 90)) +
  
  NULL

plot_aic

ggsave(plot_aic, filename = "./003__output/plot_aic_2.png",
       width = 4,
       height = 6)

# Pred by NA model --------------------------------------------------------

germany_pred <- c(germany_ts_rm, rep(NA, 12))

result_model_1_s_fluct_na_pred <- model_1_s(germany_pred, NA)

AIC_na_pred <- AIC_ammend(result_model_1_s_fluct_na_pred$kfs, degree = 3 + 1 + 11)


alphahatconf <- predict(result_model_1_s_fluct_na_pred$fit$model,
                        interval = "prediction",
                        level = 0.95) %>% as_tibble()

level <- tibble(pred_level = result_model_1_s_fluct_na_pred$kfs$alphahat[, "level"] %>% as.numeric())

x_add <- x %>%
  select(year_month, germany) %>% 
  bind_cols(alphahatconf,
            level)

plot_pred <- ggplot(x_add,
                    aes(x = year_month, y = germany)) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "lightblue",
              alpha = 0.8) +
  
  geom_line(aes(y = pred_level),
            colour = "blue",
            size = 0.1,
            alpha = 0.8) +
  
  # geom_line(alpha = 0.6) +
  
  geom_line(aes(y = upr),
            lty = "dashed",
            alpha = 0.2) +
  
  geom_line(aes(y = lwr),
            lty = "dashed",
            alpha = 0.2) +
  
  geom_point(size = 2,
             alpha = 0.6) +
  
  expand_limits(y = 0) +
  scale_y_continuous(breaks = seq(0, 30000, 5000)) +
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "2 month",
  limit=c(as.Date("2018-01-01"), as.Date("2019-12-01"))) +
  
  labs(x = "Year/Month",
       y = "Number of foreign visitors from Germany") +
  
  theme(axis.text.y = element_text(size = 10)) +
  
  NULL

plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_na_model_with_pred_from2018-.png",
       width = 6,
       height = 4)



# table of prediction -----------------------------------------------------

x_add %>% select(year_month,
                 fit,
                 lwr,
                 upr) %>% 
  filter(year_month >= as.Date("2019-01-01")) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  kable()


