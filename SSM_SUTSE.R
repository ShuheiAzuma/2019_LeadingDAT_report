
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

germany_ts <- x$germany%>%
  scale() %>%
  as_vector()

total <- x$total %>%
  scale() %>%
  as_vector()

sutse <- cbind(total = total, germany = germany_ts) %>% na.omit()

# difine AIC ammend function ----------------------------------------------

AIC_ammend <- function(x, degree){
  loglik <- x$logLik - sum(x$Finf > 0) * log(2*pi) / 2
  return(AIC_ <- - 2 * loglik + 2 * degree)
}

# difine SUTSE model ------------------------------------------------------

model_1_s_sutse <- function(x1, x2){
  
  mod <- SSModel(cbind(x1, x2) ~ SSMtrend(1, Q = matrix(NA, 2, 2)),
                 H = matrix(NA, 2, 2))
  
  fit <- fitSSM(mod, numeric(6), method = "SANN")
  kfs <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
  
  return(list(mod = mod,
              fit = fit,
              kfs = kfs))
  
}

result_sutse <- model_1_s_sutse(x1 = total, x2 = germany_ts)


# culculate Relation coefficient ------------------------------------------

# H kannsoku ti kakurannkou bunsann 

h_sutse <- result_sutse$fit$model$H[,, 1] %>% 
  as.matrix(nrow = 2)

r_h <- h_sutse[1, 2] / (sqrt(h_sutse[1, 1]) * sqrt(h_sutse[2, 2]))

# Q jyoutai kakurankou bunnsann

q_sutse <- result_sutse$fit$model$Q[,, 1] %>% 
  as.matrix(nrow = 2)

r_q <- q_sutse[1, 2] / (sqrt(q_sutse[1, 1]) * sqrt(q_sutse[2, 2]))


# AIC ---------------------------------------------------------------------

AIC_sutse <- AIC_ammend(result_sutse$kfs, degree = 6 + 1)

alphahatconf <- predict(result_sutse$fit$model, interval = "confidence", level = 0.95) 

alphahatconf_1 <- alphahatconf[[1]] %>% as_tibble()
alphahatconf_2 <- alphahatconf[[2]] %>% as_tibble()

colnames_alphahat <- result_sutse$kfs$alphahat %>% colnames()

pred <- tibble(pred_level_1 = result_sutse$kfs$alphahat[, colnames_alphahat[1]] %>% as.numeric(),
               pred_level_2 = result_sutse$kfs$alphahat[, colnames_alphahat[2]] %>% as.numeric()) %>% 

  bind_cols(year_month = x$year_month,
            total = x$total %>% scale(),
            germany = x$germany %>% scale(),
            conf_1 = alphahatconf_1,
            conf_2 = alphahatconf_2)

theme_set(theme_minimal() +
            theme(axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12,
                                           angle = 90)
            )
)

plot_pred <- ggplot(pred,
                    aes(x = year_month)) +
  
  geom_line(aes(y = pred_level_1),
            alpha = 0.6,
            color = "tomato",
            size = 2) +
  
  geom_line(aes(y = pred_level_2),
            alpha = 0.6,
            color = "purple4",
            size = 2) +
  
  # Germany
  geom_point(aes(y = total),
             alpha = 0.6,
             color = "tomato") +
  
  geom_point(aes(y = germany),
             alpha = 0.6,
             color = "purple4") +
  
  scale_x_date(date_labels = "%y/%m",
               date_breaks = "1 year") +
  
  scale_y_continuous(breaks = seq(-3, 3, 1)) +
  labs(x = "Year/Month",
       y = "Number of foreign visitors") +
  
  coord_cartesian(ylim = c(-3, 3)) +
  
  NULL

plot_pred

ggsave(plot_pred, filename = "./003__output/plot_pred_SUSTE_model_Scale.png",
       width = 6,
       height = 6)


