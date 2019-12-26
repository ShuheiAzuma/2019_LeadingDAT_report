
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

# Create a variable -------------------------------------------------------

leap <- leap_year(x$year_month)

dates <- seq(x$year_month[1], x$year_month[NROW(x)], 1)
weeks <- table(str_sub(dates, 1, 7), wday(dates, label = TRUE, locale = "US"))

sun <- weeks[, "Sun"]

calendar <- tibble(

  mon = weeks[, "Mon"] - sun,
  tue = weeks[, "Tue"] - sun,
  wed = weeks[, "Wed"] - sun,
  thu = weeks[, "Thu"] - sun,
  fri = weeks[, "Fri"] - sun,
  sat = weeks[, "Sat"] - sun
  
) %>% as.matrix()

# for Na model ----

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

result_model_1_s_fluct <- model_1_s(germany_ts, NA)

res_obs <- rstandard(result_model_1_s_fluct$kfs, "pearson")
res_state <- rstandard(result_model_1_s_fluct$kfs, "state")[, 1]


germany_res <- x %>% 
  select(year_month, germany) %>% 
  bind_cols(res_obs = res_obs,
            res_state = res_state)

remove_index <- abs(germany_res$res_obs) > 1.96

germany_ts_rm <- germany_ts
germany_ts_rm[remove_index] <- NA

# germany_ts_rm <- germany_ts_rm[c(1:192)]


# difine model ------------------------------------------------------------

model_1_s_cale <- function(x, q){
  
  mod <- SSModel(x ~ SSMtrend(1, Q = NA) +
                   SSMseasonal(12, q) +
                   leap +
                   calendar,
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

result_cale <- model_1_s_cale(germany_ts_rm, NA)

AIC_ammend <- function(x, degree){
  loglik <- x$logLik - sum(x$Finf > 0) * log(2*pi) / 2
  return(AIC_ <- - 2 * loglik + 2 * degree)
}


AIC_calendar <- AIC_ammend(result_cale$kfs, degree = 1 + 3 + 11 + 1 + 6)


pred <- tibble(pred_level = result_cale$kfs$alphahat[, "level"] %>% as.numeric()) %>% 
  mutate(calendar = result_cale$kfs$muhat - pred_level,
         pred = result_cale$kfs$muhat) %>% 
  bind_cols(year_month = x$year_month,
            germany = x$germany)



theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12,
                                           angle = 90)
            )
)

plot_pred <- ggplot(pred,
                    aes(x = year_month, y = calendar)) +
  # geom_ribbon(aes(ymin = lwr, ymax = upr),
  #             fill = "lightblue",
  #             alpha = 0.4) +
  
  geom_line(aes(y = germany),
            colour = "tomato",
            alpha = 0.6, 
            size = 2) +
  
  geom_line(alpha = 0.6) +
  
  geom_line(aes(y = pred_level),
            colour = "lightblue",
            size = 3,
            alpha = 0.8) +
  
  geom_line(aes(y = pred)) +

  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  labs(x = "Year/Month",
       y = "Number of foreign visitors") +
  NULL

plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_cale_model_sep_effect_level.png",
       width = 6,
       height = 4)



alphahatconf <- predict(result_cale$fit$model, interval = "prediction", level = 0.95) %>% as_tibble()

x_add <- x %>%
  bind_cols(alphahatconf)


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
            size = 3,
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


plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_calendar_model.png",
       width = 6,
       height = 4)
