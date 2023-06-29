### Libraries ____________________________________________________________________________
library(fBasics)
library(dplyr)
library(here)
library(fpp2)
library(xts)
library(tseries)
library(rJava)
library(RJDemetra)
library(FinTS)
library(rugarch)
library(Metrics)
library(DataCombine)
library(keras)
library(tensorflow)
library(forecast) 

# Seed
set.seed(as.integer(411974))
tf$random$set_seed(as.integer(411974))
np <- reticulate::import("numpy")
np$random$seed(as.integer(411974))

# Environment settings
options(scipen = 100)
Sys.setenv(LANG = "en")
Sys.setenv("OMP_NUM_THREADS" = "1")
Sys.setenv("TF_CPP_MIN_LOG_LEVEL" = "2")
Sys.setenv("CUDA_VISIBLE_DEVICES" = "")
keras::backend()$clear_session()

### Functions ____________________________________________________________________________
TestADF <- function(variable, adf_order){
  library(urca)
  library(lmtest)
  results_adf <- data.frame(order = -1, adf = 0,
                            p_adf = "", bgodfrey = 0, p_bg = 0)
  variable <- variable[!is.na(variable)]
  
  for(order in 0:adf_order){
    df.test_ <- ur.df(variable, type = c("drift"), lags = order)
    df_ <- df.test_@teststat[1]
    df_crit <- df.test_@cval[1, ]
    df_crit <- (df_ < df_crit) * 1
    p_adf <- ifelse (sum(df_crit) == 0,
                     ">10pct",
                     paste("<", names(df_crit)[min(which(df_crit == 1))],
                           sep = ""))
    resids_ <- df.test_@testreg$residuals
    bgtest_ <- bgtest(resids_ ~ 1, order = 1)
    bgodfrey <- bgtest_$statistic
    names(bgodfrey)<-NULL
    p_bg <- bgtest_$p.value
    
    results_adf <- rbind(results_adf,
                         data.frame(order = order,
                                    adf = df_,
                                    p_adf = p_adf,
                                    bgodfrey = bgodfrey,
                                    p_bg = p_bg)
    )
  }
  
  results_adf<-results_adf[results_adf$order>=0,]
  
  plot(variable,
       type = "l",
       col = "blue",
       lwd = 2,
       main = "Plot of the examined variable")
  
  return(results_adf)
}

EstimateARIMA <- function(p, q, n_opt = n_in){
  
  arima <- arima(in_sample$CPI_diff_log[(n_in - n_opt + 1):n_in],
                 order = c(p, 0, q),
                 include.mean = FALSE) 
  
  if(n_opt == n_in){
    model <- paste0("arima", p, 1, q) 
  } else {
    model <- paste0("arima_opt", p, 1, q)
  }
  
  forecasts[, model] <<- c(rep(NA, (n_in - n_opt)), c(data$CPI_diff_log[(n_in - n_opt + 1):n_in] - residuals(arima), rep(NA, n_out)))
  
  # rolling window
  for (j in 1:n_out){
    
    # time window = n_opt
    input <- data[(n_in - n_opt + j):(n_in + j - 1), ]
    print(paste("data is from", n_in - n_opt + j, "to", n_in + j - 1))
    
    forecast <- arima(input$CPI_diff_log,
                      order = c(p, 0, q),
                      include.mean = FALSE,
                      method = "ML") %>% 
      forecast(h = 1)
    forecasts[n_in + j, model] <<- forecast$mean[1]
    print(paste("saving forecast number", n_in + j))
    
  }
  
  i <- sum(!is.na(forecast_errors$method)) + 1
  forecast_errors$method[i] <<- model
  forecast_errors$mae[i] <<- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model], n = n_out))
  forecast_errors$rmse[i] <<- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model], n = n_out))
  
}

EstimateGARCH <- function(p, q, P, Q, distr, mean = FALSE, n_opt = n_in){
  
  if(mean == FALSE){
    spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                             garchOrder = c(p, q)),
                       mean.model = list(armaOrder = c(P, Q),
                                         include.mean = FALSE),
                       distribution.model = distr) 
    if(n_opt == n_in){
      model <- paste0("garch", p, q, "_arma", P, Q, "_", distr)
    } else {
      model <- paste0("garch_opt", p, q, "_arma", P, Q, "_", distr)
    }
    
  } else if (mean == TRUE){
    spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                             garchOrder = c(p, q)),
                       mean.model = list(armaOrder = c(P, Q),
                                         archm = TRUE, 
                                         archpow = 1,
                                         include.mean = FALSE),
                       distribution.model = distr) 
    
    if(n_opt == n_in){
      model <- paste0("garch_m", p, q, "_arma", P, Q, "_", distr)
    } else {
      model <- paste0("garch_m_opt", p, q, "_arma", P, Q, "_", distr)
    }
  }
  
  garch.fit <- ugarchfit(spec = spec,
                         data = na.omit(data$CPI_diff_log[(n_in - n_opt + 1):(n_in + n_out)]), 
                         out.sample = n_out)
  
  forecasts[, model] <<- c(rep(NA, (n_in - n_opt)), c(garch.fit@fit$fitted.values, rep(NA, n_out)))
  
  # rolling window
  for (j in 1:n_out){
    
    # time window = n_out
    input <- data[(n_in - n_opt + j):(n_in + j - 1), ]
    print(paste("data is from", n_in - n_opt + j, "to", n_in + j - 1))
    
    garch.fit <- ugarchfit(spec = spec,
                           data = input$CPI_diff_log) 
    forecast <- ugarchforecast(garch.fit, n.ahead = 1)
    
    forecasts[n_in + j, model] <<- forecast@forecast$seriesFor
    print(paste("saving forecast number", n_in + j))
    
  }
  
  i <- sum(!is.na(forecast_errors$method)) + 1
  forecast_errors$method[i] <<- model
  forecast_errors$mae[i] <<- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model], n = n_out))
  forecast_errors$rmse[i] <<- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model], n = n_out))
  
}

EstimateNN <- function(type = "rnn", window_size, batch_size, units, activation, units_two = NULL, activation_two = NULL, epochs = 20){
  
  # Preparing the matrix input
  exch_matrix <- matrix(0, nrow = length(scaled_train_data) - window_size, ncol = window_size + 1)
  for(i in 1:nrow(exch_matrix)){
    exch_matrix[i, ] <- scaled_train_data[i:(i + window_size)]
  }
  
  # Splitting the matrix into lags and target values
  x_train <- exch_matrix[, -ncol(exch_matrix)]
  y_train <- exch_matrix[, ncol(exch_matrix)]
  
  # Reshaping the data
  x_train <- array_reshape(x_train, dim = c((length(scaled_train_data) - window_size), window_size, 1))
  
  # Data for forecast calculation
  all_data <- append(scaled_train_data, scaled_test_data)
  exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
  for(i in 1:nrow(exch_matrix_all)){
    exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
  }
  
  x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
  y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
  x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
  
  # Model architecture
  model <- keras_model_sequential()
  
  if(is.null(units_two) & is.null(activation_two)){
    
    if(type == "rnn"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_simple_rnn(units = units, activation = activation, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
      
    } else if (type == "lstm"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_lstm(units = units, activation = activation, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
      
    } else if (type == "gru"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_gru(units = units, activation = activation, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
    }
    
    model_name <- paste0(type, "1_", window_size, activation, units)
    
  } else{
    
    if(type == "rnn"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_simple_rnn(units = units, activation = activation, kernel_initializer = "glorot_uniform", return_sequences = TRUE) %>% 
        layer_simple_rnn(units = units_two, activation = activation_two, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
      
    } else if (type == "lstm"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_lstm(units = units, activation = activation, kernel_initializer = "glorot_uniform", return_sequences = TRUE) %>% 
        layer_lstm(units = units_two, activation = activation_two, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
      
    } else if (type == "gru"){
      model %>% 
        layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_initializer = "glorot_uniform",  kernel_regularizer = regularizer_l2(0.01)) %>% 
        layer_gru(units = units, activation = activation, kernel_initializer = "glorot_uniform", return_sequences = TRUE) %>% 
        layer_gru(units = units_two, activation = activation_two, kernel_initializer = "glorot_uniform") %>% 
        layer_dense(units = 1, kernel_initializer = "glorot_uniform")
    }
    
    model_name <- paste0(type, "2_", window_size, activation, units)
  }
  
  # Model training
  model %>% compile(
    loss = "mse",
    optimizer= "adam",
    metric = "mae" 
  )
  
  history <- model %>% 
    fit(x_train, y_train, epochs = epochs, batch_size = batch_size)
  
  model_name <- paste0(model_name, "_l2")
  
  # Predictions for all observations
  forecasts_all_scaled <- model %>% predict(x_train_all)
  forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
  forecasts[, model_name] <<- append(rep(NA, window_size), forecasts_all)
  
  i <- sum(!is.na(forecast_errors$method)) + 1
  forecast_errors$method[i] <<- model_name
  forecast_errors$mae[i] <<- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model_name], n = n_out))
  forecast_errors$rmse[i] <<- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts[, model_name], n = n_out))
  
  keras::backend()$clear_session()
}

### Data loading ____________________________________________________________________________
data <- read.csv(file = here("CPIAUCNS.csv")) %>%
  mutate(DATE = as.Date(DATE, format = "%Y-%m-%d")) %>%
  filter(DATE >= "1913-01-01") %>%
  rename("Date" = "DATE",
         "CPI" = "CPIAUCNS")

### Data visualization ____________________________________________________________________________
{
data %>%
  ggplot(aes(x = Date,
             y = CPI)) + 
  geom_line() +
  ggtitle("US CPI monthly")

# Differenced logarythim of unseasoned CPI - RJDemetra
tramo <- tramoseats(ts(data$CPI, start = c(1913, 1), frequency = 12), spec = "RSAfull")
CPI_unseasoned <- tramo[["final"]][["series"]][, 2]
data$CPI_diff_log <- 1200 * diff.xts(log(CPI_unseasoned))

data %>%
  ggplot(aes(x = Date,
             y = CPI_diff_log)) + 
  geom_line() +
  ggtitle("US CPI diff log monthly seasoned")

data %>%
  ggplot(aes(x = CPI_diff_log)) + 
  geom_histogram(aes(y = ..density..), bins = 50) +
  geom_density(col = "red") +
  ggtitle("Histogram of differenciated logarithm of US CPI monthly seasoned") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  xlab("US CPI diff log monthly seasoned")

data <- na.omit(data)
}
### Data preparation and analysis ____________________________________________________________________________
{
# Basic statistics
basicStats(data$CPI_diff_log) # excess kurtosis > 0 (leptokurtic distribution)

# Normality test
jarque.bera.test(na.omit(data$CPI_diff_log)) # p-value < 5%, H0 about normality rejected

# Stationarity
TestADF(variable = data$CPI_diff_log, adf_order = 3) # stationary

pp.test <- ur.pp(data$CPI_diff_log, type = c("Z-tau"), model = c("constant"))
summary(pp.test) # Z-tau less than 5% critical value - H0 about non-stationarity rejected - stationary

kpss.test <- ur.kpss(data$CPI_diff_log, type = c("mu")) # type = c("mu", "tau")
summary(kpss.test) # statistic less than 5% critical value  - H0 about stationarity not rejected - stationary

# Data splitting
in_sample = data[which(data$Date < "2007-01-01"), ] # ~ 85%
out_of_sample = data[-which(data$Date < "2007-01-01"), ] # ~ 15%
n_in <- nrow(in_sample) # 1127
n_out <- nrow(out_of_sample) # 196
n_val <- nrow(data[which(data$Date < "1991-01-01"), ]) # 935 
# train set: n_val = 935, val set: n_in - n_val = 192, test set: n_out = 196

# ACF and PACF
par(mfrow = c(2, 1))
acf(in_sample$CPI_diff_log, lag.max = 36,
    ylim = c(-0.15, 0.5),
    lwd = 7, col = "dark green") 
pacf(in_sample$CPI_diff_log, lag.max = 36,       
     lwd = 7, col = "dark green")
par(mfrow = c(1, 1))

# Results tables preparation
method <- rep(NA, 30)
forecasts <- data.frame(Date = data$Date, CPI_diff_log = data$CPI_diff_log)
forecast_errors <- data.frame(method, "mae" = NA, "rmse" = NA)
}
### Simple models ____________________________________________________________________________
{
# Historical average
forecasts$HA <- mean(na.omit(in_sample$CPI_diff_log))

forecast_errors$method[1] <- "HA"
forecast_errors$mae[1] <- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$HA, n = n_out))
forecast_errors$rmse[1] <- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$HA, n = n_out))

# Moving average 1 year
forecasts$MA1y <- rollmean(data$CPI_diff_log, 12, fill = NA, align = "right") %>% shift(-1, reminder = FALSE)

i <- sum(!is.na(forecast_errors$method)) + 1
forecast_errors$method[i] <- "MA1y" 
forecast_errors$mae[i] <- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA1y, n = n_out))
forecast_errors$rmse[i] <- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA1y, n = n_out))

# Moving average 3 years
forecasts$MA3y <- rollmean(data$CPI_diff_log, 36, fill = NA, align = "right") %>% shift(-1, reminder = FALSE)

i <- sum(!is.na(forecast_errors$method)) + 1
forecast_errors$method[i] <- "MA3y" 
forecast_errors$mae[i] <- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA3y, n = n_out))
forecast_errors$rmse[i] <- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA3y, n = n_out))

# Moving average 6 years
forecasts$MA6y <- rollmean(data$CPI_diff_log, 72, fill = NA, align = "right") %>% shift(-1, reminder = FALSE)

i <- sum(!is.na(forecast_errors$method)) + 1
forecast_errors$method[i] <- "MA6y" 
forecast_errors$mae[i] <- mae(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA6y, n = n_out))
forecast_errors$rmse[i] <- rmse(tail(forecasts$CPI_diff_log, n = n_out), tail(forecasts$MA6y, n = n_out))
}
### ARIMA models ____________________________________________________________________________
{
# Auto ARIMA
auto_arima <- auto.arima(in_sample$CPI_diff_log, seasonal = FALSE, d = 0)
auto_arima # 102

Box.test(resid(auto_arima), type = "Ljung-Box", lag = 24) # p-value < 5%, H0 rejected - residuals are not white noise
ArchTest(resid(auto_arima), lags = 24) # p-value < 5%, H0 rejected - ARCH effects present
  
# Finding the best model
if(0){
  method <- rep(NA, 25)
  forecasts_arima <- data.frame(Date = in_sample$Date, CPI_diff_log = in_sample$CPI_diff_log)
  forecast_errors_arima <- data.frame(method, "mae" = NA, "rmse" = NA, "aic" = NA, "bic" = NA)
  for (p in c(0, 1, 2, 3, 4)){
    for (q in c(0, 1, 2, 3, 4)){
      
      arima <- arima(in_sample$CPI_diff_log[1:n_val],
                     order = c(p, 0, q),
                     include.mean = FALSE) 
      
      model <- paste0("arima", p, 1, q)
      
      forecasts_arima[, model] <- c(data$CPI_diff_log[1:n_val] - residuals(arima), rep(NA,  (n_in - n_val)))
      
      k <- 1
      for (i in 1:(n_in - n_val)){
        
        # for now time window = n_val
        input <- data[i:(n_val + i), ]
        print(paste("data is from", i, "to", n_val + i))
        
        forecast <- arima(input$CPI_diff_log,
                          order = c(p, 0, q),
                          include.mean = FALSE) %>% 
          forecast(h = k)
        
        forecasts_arima[n_val + i, model] <- forecast$mean[k]
        print(paste("saving forecast number", n_val + i))
        
      }
      
      i <- sum(!is.na(forecast_errors_arima$method)) + 1
      forecast_errors_arima$method[i] <- model
      forecast_errors_arima$mae[i] <- mae(forecasts_arima$CPI_diff_log[(n_val + 1):n_in], forecasts_arima[(n_val + 1):n_in, model])
      forecast_errors_arima$rmse[i] <- rmse(forecasts_arima$CPI_diff_log[(n_val + 1):n_in], forecasts_arima[(n_val + 1):n_in, model])
      forecast_errors_arima$aic[i] <- arima$aic
      forecast_errors_arima$bic[i] <- AIC(arima, k = log(n_val))
    }
  }
  saveRDS(forecast_errors_arima, here("forecast_errors_arima.RData"))
} # 110, 313
forecast_errors_arima <- readRDS(here("forecast_errors_arima.RData"))

# Best ARIMA estimation 
EstimateARIMA(p = 1, q = 0) # rmse, mae
EstimateARIMA(p = 3, q = 3) # aic, bic

# Optimal window size
if(0){
  min_size <- 60
  max_size <- 900
  window <- rep(NA, max_size - min_size + 1)
  rmse <- rep(NA, max_size - min_size + 1)
  mae <- rep(NA, max_size - min_size + 1)
  forecasts_arima_window <- data.frame(Date = data$Date, CPI_diff_log = data$CPI_diff_log, arima = NA)
  forecast_errors_arima_window <- data.frame(window, rmse, mae)
  
  iterator <- 1
  
  for (n in max_size:min_size) {
    print(paste("window size:", n))
    
    for (i in 1:(n_in - n_val)) {
      
      input <- data[(n_val - n + i):(n_val + i), ]
      
      # ARIMA
      try(arima <- arima(input$CPI_diff_log, order = c(3, 0, 3), method = "ML") %>% 
            forecast(h = 1))
      try(forecasts_arima_window$arima[n_val + i]<- arima$mean[k])
    }
    
    forecast_errors_arima_window$window[iterator] <- n
    forecast_errors_arima_window$rmse[iterator] <- rmse(forecasts_arima_window$CPI_diff_log[(n_val + 1):n_in], forecasts_arima_window$arima[(n_val + 1):n_in]) 
    forecast_errors_arima_window$mae[iterator] <- mae(forecasts_arima_window$CPI_diff_log[(n_val + 1):n_in], forecasts_arima_window$arima[(n_val + 1):n_in])
    
    iterator <- iterator + 1
  } 
  saveRDS(forecast_errors_arima_window, here("forecast_errors_arima_window.RData"))
}
forecast_errors_arima_window <- readRDS(here("forecast_errors_arima_window.RData"))
plot(forecast_errors_arima_window$window, forecast_errors_arima_window$mae) # 337
plot(forecast_errors_arima_window$window, forecast_errors_arima_window$rmse) # 337

# Best ARIMA estimation with optimal window
EstimateARIMA(p = 1, q = 0, n_opt = 250)
EstimateARIMA(p = 3, q = 3, n_opt = 250)
}
### GARCH models ____________________________________________________________________________
{
# Finding the best model
if(0){
  start <- Sys.time()
  method <- rep(NA, 864)
  forecasts_garch <- data.frame(Date = in_sample$Date, CPI_diff_log = in_sample$CPI_diff_log)
  forecast_errors_garch <- data.frame(method, "mae" = NA, "rmse" = NA, "aic" = NA, "bic" = NA)
  for (distr in c("norm", "ged", "std")){
    for (p in c(0, 1, 2)){
      for (q in c(0, 1, 2)){
        for (P in c(0, 1, 2, 3)){
          for (Q in c(0, 1, 2, 3)){
            for (mean in c(TRUE, FALSE)){
              
              if(p == 0 & q == 0 & P == 0 & Q == 0){
                break
              }
              
              if(mean == FALSE){
                spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                                         garchOrder = c(p, q)),
                                   mean.model = list(armaOrder = c(P, Q),
                                                     include.mean = FALSE),
                                   distribution.model = distr) 
                
                model <- paste0("garch", p, q, "_arma", P, Q, "_", distr)
                
              } else if (mean == TRUE){
                spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                                         garchOrder = c(p, q)),
                                   mean.model = list(armaOrder = c(P, Q),
                                                     archm = TRUE, 
                                                     archpow = 1,
                                                     include.mean = FALSE),
                                   distribution.model = distr) 
                
                model <- paste0("garch_m", p, q, "_arma", P, Q, "_", distr)
              }
              
              print(paste("estimation of", model))
              
              check <- try(garch.fit <- ugarchfit(spec = spec,
                                                  data = na.omit(in_sample$CPI_diff_log),
                                                  out.sample = (n_in - n_val),
                                                  solver = 'hybrid'), silent = TRUE)
              
              if(class(check) != "try-error"){
                
                check <- try(forecasts_out <- ugarchforecast(garch.fit, n.ahead = 1, n.roll = n_in - n_val - 1), silent = TRUE)
                if(class(check) != "try-error"){
                  forecasts_garch[, model] <- append(garch.fit@fit$fitted.values, forecasts_out@forecast$seriesFor)
                } else{break}
                
                if(!any(is.na(forecasts_garch[, model]))){
                  i <- sum(!is.na(forecast_errors_garch$method)) + 1
                  forecast_errors_garch$method[i] <- model
                  forecast_errors_garch$mae[i] <- mae(forecasts_garch$CPI_diff_log[(n_val + 1):n_in], forecasts_garch[(n_val + 1):n_in, model])
                  forecast_errors_garch$rmse[i] <- rmse(forecasts_garch$CPI_diff_log[(n_val + 1):n_in], forecasts_garch[(n_val + 1):n_in, model])
                  forecast_errors_garch$aic[i] <- infocriteria(garch.fit)[1]
                  forecast_errors_garch$bic[i] <- infocriteria(garch.fit)[2]
                }
              }
            }
          }
        }
      }
    }
  }
  end <- Sys.time()
  print(end - start)
  saveRDS(forecast_errors_garch, here("forecast_errors_garch.RData"))
}
forecast_errors_garch <- readRDS(here("forecast_errors_garch.RData")) 

# Best GARCH estimation
EstimateGARCH(0, 2, 0, 1, "std", TRUE) # mae
EstimateGARCH(0, 1, 0, 1, "norm", TRUE) # rmse
try(EstimateGARCH(1, 0, 3, 1, "norm", FALSE)) # mae, rmse - standard garch
EstimateGARCH(1, 0, 0, 0, "ged", FALSE) # aic, bic (no mean model)

# Optimal window size
if(0){
  min_size <- 12
  max_size <- 900
  window <- rep(NA, max_size - min_size + 1)
  rmse <- rep(NA, max_size - min_size + 1)
  mae <- rep(NA, max_size - min_size + 1)
  forecasts_garch_window <- data.frame(Date = data$Date, CPI_diff_log = data$CPI_diff_log, garch = NA)
  forecast_errors_garch_window <- data.frame(window, rmse, mae)
  
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(0, 1)),
                     mean.model = list(armaOrder = c(0, 1),
                                       archm = TRUE, 
                                       archpow = 1, # 1 - standard deviation, 2 - variance
                                       include.mean = FALSE), 
                     distribution.model = "std")
  
  iterator <- 1
  for (n in max_size:min_size) {
    
    print(paste("window size:", n))
    input <- data[(n_val - n + 1):n_in, ]
    
    check <- try(garch.fit <- ugarchfit(spec = spec,
                                        data = na.omit(input$CPI_diff_log),
                                        out.sample = (n_in - n_val),
                                        solver = 'hybrid'), silent = TRUE)
    
    if(class(check) != "try-error"){
      
      check <- try(forecasts_out <- ugarchforecast(garch.fit, n.ahead = 1, n.roll = n_in - n_val - 1), silent = TRUE)
      if(class(check) != "try-error"){
        forecasts_garch_window$garch <- c(rep(NA, n_val - n), garch.fit@fit$fitted.values, forecasts_out@forecast$seriesFor, rep(NA, n_out))
      } else{break}
      
      forecast_errors_garch_window$window[iterator] <- n 
      forecast_errors_garch_window$rmse[iterator] <- rmse(forecasts_garch_window$CPI_diff_log[(n_val + 1):n_in], forecasts_garch_window$garch[(n_val + 1):n_in]) 
      forecast_errors_garch_window$mae[iterator] <- mae(forecasts_garch_window$CPI_diff_log[(n_val + 1):n_in], forecasts_garch_window$garch[(n_val + 1):n_in])
    }
    
    iterator <- iterator + 1
  } 
  saveRDS(forecast_errors_garch_window, here("forecast_errors_garch_window.RData"))
}
forecast_errors_garch_window <- readRDS(here("forecast_errors_garch_window.RData"))
plot(forecast_errors_garch_window$window, forecast_errors_garch_window$mae) # 777
plot(forecast_errors_garch_window$window, forecast_errors_garch_window$rmse) # 882

# Best GARCH estimation with optimal window
EstimateGARCH(0, 2, 0, 1, "std", TRUE, n_opt = 880)
EstimateGARCH(0, 1, 0, 1, "norm", TRUE, n_opt = 880)
try(EstimateGARCH(1, 0, 3, 1, "norm", FALSE, n_opt = 880))
EstimateGARCH(1, 0, 0, 0, "ged", FALSE, n_opt = 880) 
}
### Neural networks ____________________________________________________________________________
{
# Data preparation
train_data <- ts(in_sample$CPI_diff_log)
test_data <- ts(out_of_sample$CPI_diff_log)

max_train <- max(train_data)
min_train <- min(train_data)
scaled_train_data <- (train_data - min_train)/(max_train - min_train)
scaled_test_data <- (test_data - min_train)/(max_train - min_train)

if(0){
  rnn_one_errors <- data.frame("method" = rep(NA, 1152), "mae" = NA, "rmse" = NA)
  for (window_size in c(12, 24, 36, 60, 120, 240)){
    for (units in c(3, 5, 10)){
      for (batch_size in c(16, 32, 64, 128)){
        for (activation in c("relu", "sigmoid", "tanh")){
          for (l2 in c(T, F)){
            for(bn in c(T, F)){
              
              # Data preparation
              val_data <- scaled_train_data[1:n_val]
              exch_matrix <- matrix(0, nrow = length(val_data) - window_size, ncol = window_size + 1)
              for(i in 1:nrow(exch_matrix)){
                exch_matrix[i, ] <- val_data[i:(i + window_size)]
              }
              
              x_train <- exch_matrix[, -ncol(exch_matrix)]
              y_train <- exch_matrix[, ncol(exch_matrix)]
              x_train <- array_reshape(x_train, dim = c((length(val_data) - window_size), window_size, 1))
              
              all_data <- scaled_train_data
              exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
              for(i in 1:nrow(exch_matrix_all)){
                exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
              }
              
              x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
              y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
              x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
              
              # Neural network
              model <- keras_model_sequential()
              
              if(l2 == T & bn == T){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_regularizer = regularizer_l2()) %>% 
                  layer_simple_rnn(units = units, activation = activation) %>% 
                  layer_batch_normalization() %>%
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_l2_bn")
              } else if(l2 == T & bn == F){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_regularizer = regularizer_l2(0.01)) %>% 
                  layer_simple_rnn(units = units, activation = activation) %>% 
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_l2")
              } else if(l2 == F & bn == T){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_simple_rnn(units = units, activation = activation) %>% 
                  layer_batch_normalization() %>%
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_bn")
              } else{
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_simple_rnn(units = units, activation = activation) %>% 
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation)
              }
              
              summary(model)
              model %>% compile(
                loss = "mse",
                optimizer= "adam",
                metric = "mae" 
              )
              
              history <- model %>% 
                fit(x_train, y_train, epochs = 20, batch_size = batch_size)
              
              forecasts_all_scaled <- model %>% predict(x_train_all)
              forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
              
              i <- sum(!is.na(rnn_one_errors$method)) + 1
              rnn_one_errors$method[i] <- model_name
              rnn_one_errors$mae[i] <- mae(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              rnn_one_errors$rmse[i] <- rmse(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              
            }
          }
        }
      }
    }
  }
  saveRDS(rnn_one_errors, here("rnn_one_errors.RData"))
  
  rnn_two_errors <- data.frame("method" = rep(NA, 1215), "mae" = NA, "rmse" = NA)
  for (window_size in c(24, 36, 60, 120, 240)){
    for (units_one in c(5, 10, 15))
      for (units_two in c(2, 5, 10)){
        for (batch_size in c(32, 64, 128)){
          for (activation_one in c("relu", "sigmoid", "tanh")){
            for (activation_two in c("relu", "sigmoid", "tanh")){
              
              if (units_one > units_two){
                
                # Data preparation
                val_data <- scaled_train_data[1:n_val]
                exch_matrix <- matrix(0, nrow = length(val_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix)){
                  exch_matrix[i, ] <- val_data[i:(i + window_size)]
                }
                
                x_train <- exch_matrix[, -ncol(exch_matrix)]
                y_train <- exch_matrix[, ncol(exch_matrix)]
                x_train <- array_reshape(x_train, dim = c((length(val_data) - window_size), window_size, 1))
                
                all_data <- scaled_train_data
                exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix_all)){
                  exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
                }
                
                x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
                y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
                x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
                
                # Neural network
                model <- keras_model_sequential()
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_simple_rnn(units = units_one, activation = activation_one, return_sequences = TRUE) %>% 
                  layer_simple_rnn(units = units_two, activation = activation_two) %>% 
                  layer_dense(units = 1)
                summary(model)
                
                model %>% compile(
                  loss = "mse",
                  optimizer= "adam",
                  metric = "mae" 
                )
                
                history <- model %>% 
                  fit(x_train, y_train, epochs = 20, batch_size = batch_size)
                
                forecasts_all_scaled <- model %>% predict(x_train_all)
                forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
                
                i <- sum(!is.na(rnn_two_errors$method)) + 1
                rnn_two_errors$method[i] <- paste(window_size, units_one, activation_one, units_two, activation_two, batch_size)
                rnn_two_errors$mae[i] <- mae(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
                rnn_two_errors$rmse[i] <- rmse(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              }
            }
          }
        }
      }
    }
  saveRDS(rnn_two_errors, here("rnn_two_errors.RData"))
  
  lstm_two_errors <- data.frame("method" = rep(NA, 1215), "mae" = NA, "rmse" = NA)
  for (window_size in c(24, 36, 60, 120, 240)){
    for (units_one in c(5, 10, 15))
      for (units_two in c(2, 5, 10)){
        for (batch_size in c(32, 64, 128)){
          for (activation_one in c("relu", "sigmoid", "tanh")){
            for (activation_two in c("relu", "sigmoid", "tanh")){
              
              if (units_one > units_two){
                
                # Data preparation
                val_data <- scaled_train_data[1:n_val]
                exch_matrix <- matrix(0, nrow = length(val_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix)){
                  exch_matrix[i, ] <- val_data[i:(i + window_size)]
                }
                
                x_train <- exch_matrix[, -ncol(exch_matrix)]
                y_train <- exch_matrix[, ncol(exch_matrix)]
                x_train <- array_reshape(x_train, dim = c((length(val_data) - window_size), window_size, 1))
                
                all_data <- scaled_train_data
                exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix_all)){
                  exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
                }
                
                x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
                y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
                x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
                
                # Neural network
                model <- keras_model_sequential()
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_lstm(units = units_one, activation = activation, return_sequences = TRUE) %>% 
                  layer_lstm(units = units_two, activation = activation) %>% 
                  layer_dense(units = 1)
                summary(model)
                
                model %>% compile(
                  loss = "mse",
                  optimizer= "adam",
                  metric = "mae" 
                )
                
                history <- model %>% 
                  fit(x_train, y_train, epochs = 20, batch_size = batch_size)
                
                forecasts_all_scaled <- model %>% predict(x_train_all)
                forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
                
                i <- sum(!is.na(lstm_two_errors$method)) + 1
                lstm_two_errors$method[i] <-  paste(window_size, units_one, activation_one, units_two, activation_two, batch_size)
                lstm_two_errors$mae[i] <- mae(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
                lstm_two_errors$rmse[i] <- rmse(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              }
            }
          }
        }
      }
  }
  saveRDS(lstm_two_errors, here("lstm_two_errors.RData"))
  
  gru_one_errors <- data.frame("method" = rep(NA, 1152), "mae" = NA, "rmse" = NA)
  for (window_size in c(12, 24, 36, 60, 120, 240)){
    for (units in c(3, 5, 10)){
      for (batch_size in c(16, 32, 64, 128)){
        for (activation in c("relu", "sigmoid", "tanh")){
          for (l2 in c(T, F)){
            for(bn in c(T, F)){
              
              # Data preparation
              val_data <- scaled_train_data[1:n_val]
              exch_matrix <- matrix(0, nrow = length(val_data) - window_size, ncol = window_size + 1)
              for(i in 1:nrow(exch_matrix)){
                exch_matrix[i, ] <- val_data[i:(i + window_size)]
              }
              
              x_train <- exch_matrix[, -ncol(exch_matrix)]
              y_train <- exch_matrix[, ncol(exch_matrix)]
              x_train <- array_reshape(x_train, dim = c((length(val_data) - window_size), window_size, 1))
              
              all_data <- scaled_train_data
              exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
              for(i in 1:nrow(exch_matrix_all)){
                exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
              }
              
              x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
              y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
              x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
              
              # Neural network
              model <- keras_model_sequential()
              
              if(l2 == T & bn == T){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_regularizer = regularizer_l2()) %>% 
                  layer_gru(units = units, activation = activation) %>% 
                  layer_batch_normalization() %>%
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_l2_bn")
              } else if(l2 == T & bn == F){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size, kernel_regularizer = regularizer_l2(0.01)) %>% 
                  layer_gru(units = units, activation = activation) %>% 
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_l2")
              } else if(l2 == F & bn == T){
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_gru(units = units, activation = activation) %>% 
                  layer_batch_normalization() %>%
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation, "_bn")
              } else{
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_gru(units = units, activation = activation) %>% 
                  layer_dense(units = 1)
                model_name <- paste(window_size, units, batch_size, activation)
              }
              
              summary(model)
              model %>% compile(
                loss = "mse",
                optimizer= "adam",
                metric = "mae" 
              )
              
              history <- model %>% 
                fit(x_train, y_train, epochs = 20, batch_size = batch_size)
              
              forecasts_all_scaled <- model %>% predict(x_train_all)
              forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
              
              i <- sum(!is.na(gru_one_errors$method)) + 1
              gru_one_errors$method[i] <- model_name
              gru_one_errors$mae[i] <- mae(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              gru_one_errors$rmse[i] <- rmse(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              
            }
          }
        }
      }
    }
  }
  saveRDS(gru_one_errors, here("gru_one_errors.RData"))
  
  gru_two_errors <- data.frame("method" = rep(NA, 1215), "mae" = NA, "rmse" = NA)
  for (window_size in c(24, 36, 60, 120, 240)){
    for (units_one in c(5, 10, 15))
      for (units_two in c(2, 5, 10)){
        for (batch_size in c(32, 64, 128)){
          for (activation_one in c("relu", "sigmoid", "tanh")){
            for (activation_two in c("relu", "sigmoid", "tanh")){
              
              if (units_one > units_two){
                
                # Data preparation
                val_data <- scaled_train_data[1:n_val]
                exch_matrix <- matrix(0, nrow = length(val_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix)){
                  exch_matrix[i, ] <- val_data[i:(i + window_size)]
                }
                
                x_train <- exch_matrix[, -ncol(exch_matrix)]
                y_train <- exch_matrix[, ncol(exch_matrix)]
                x_train <- array_reshape(x_train, dim = c((length(val_data) - window_size), window_size, 1))
                
                all_data <- scaled_train_data
                exch_matrix_all <- matrix(0, nrow = length(all_data) - window_size, ncol = window_size + 1)
                for(i in 1:nrow(exch_matrix_all)){
                  exch_matrix_all[i, ] <- all_data[i:(i + window_size)]
                }
                
                x_train_all <- exch_matrix_all[, -ncol(exch_matrix_all)]
                y_train_all <- exch_matrix_all[, ncol(exch_matrix_all)]
                x_train_all <- array_reshape(x_train_all, dim = c((length(all_data) - window_size), window_size, 1))
                
                # Neural network
                model <- keras_model_sequential()
                model %>% 
                  layer_dense(input_shape = dim(x_train)[-1], units = window_size) %>% 
                  layer_gru(units = units_one, activation = activation_one, return_sequences = TRUE) %>% 
                  layer_gru(units = units_two, activation = activation_two) %>% 
                  layer_dense(units = 1)
                summary(model)
                
                model %>% compile(
                  loss = "mse",
                  optimizer= "adam",
                  metric = "mae" 
                )
                
                history <- model %>% 
                  fit(x_train, y_train, epochs = 20, batch_size = batch_size)
                
                forecasts_all_scaled <- model %>% predict(x_train_all)
                forecasts_all <- forecasts_all_scaled * (max_train - min_train) + min_train
                
                i <- sum(!is.na(gru_two_errors$method)) + 1
                gru_two_errors$method[i] <- paste(window_size, units_one, activation_one, units_two, activation_two, batch_size)
                gru_two_errors$mae[i] <- mae(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
                gru_two_errors$rmse[i] <- rmse(forecasts$CPI_diff_log[(n_val + 1):n_in], tail(forecasts_all, n = n_in - n_val))
              }
            }
          }
        }
      }
    }
  saveRDS(gru_two_errors, here("gru_two_errors.RData"))
}
rnn_one_errors <- readRDS(here("rnn_one_errors.RData")) # 120 3 32 tanh
rnn_two_errors <- readRDS(here("rnn_two_errors.RData")) # 24 10 sigmoid 5 tanh 32, 36 5 sigmoid 2 relu 32
lstm_one_errors <- readRDS(here("lstm_one_errors.RData")) # 60 5 32 sigmoid
lstm_two_errors <- readRDS(here("lstm_two_errors.RData")) # 120 10 tanh 2 sigmoid 64?, 60 5 relu 2 sigmoid 64
gru_one_errors <- readRDS(here("gru_one_errors.RData")) # 36 10 32 sigmoid

# Estimation of best models
EstimateNN(type = "rnn", 
           window_size = 120, batch_size = 32, 
           units = 3, activation = "tanh", 
           units_two = NULL, activation_two = NULL, 
           epochs = 25) 

EstimateNN(type = "rnn", 
           window_size = 24, batch_size = 32, 
           units = 10, activation = "sigmoid", 
           units_two = 5, activation_two = "tanh", 
           epochs = 25)

EstimateNN(type = "rnn", 
           window_size = 36, batch_size = 32, 
           units = 5, activation = "sigmoid", 
           units_two = 2, activation_two = "relu", 
           epochs = 25)

EstimateNN(type = "lstm", 
           window_size = 60, batch_size = 32, 
           units = 5, activation = "sigmoid", 
           units_two = NULL, activation_two = NULL, 
           epochs = 25)

EstimateNN(type = "lstm", 
           window_size = 120, batch_size = 64, 
           units = 10, activation = "tanh", 
           units_two = 2, activation_two = "sigmoid", 
           epochs = 25)

EstimateNN(type = "lstm", 
           window_size = 60, batch_size = 64, 
           units = 5, activation = "relu", 
           units_two = 2, activation_two = "sigmoid", 
           epochs = 25)

EstimateNN(type = "gru", 
           window_size = 36, batch_size = 32, 
           units = 10, activation = "sigmoid", 
           units_two = NULL, activation_two = NULL, 
           epochs = 25)
}
### Forecasts visualization ____________________________________________________________________________
forecasts_plot <- ggplot(forecasts, aes(x = Date)) + 
  geom_line(aes(y = CPI_diff_log), color = "black") + 
  geom_line(aes(y = MA6y), color = "blue") + 
  geom_line(aes(y = garch_m01_arma01_norm), color = "red") + 
  geom_line(aes(y = lstm1_60sigmoid5), color = "green") +
  geom_line(aes(y = arima_opt313), color = "yellow") 
forecasts_plot

### DM test ____________________________________________________________________________
errors <- data.frame(Date = data$Date, CPI_diff_log = data$CPI_diff_log)
errors$ARIMA <- forecasts$CPI_diff_log - forecasts$arima_opt313
errors$NN <- forecasts$CPI_diff_log - forecasts$lstm1_60sigmoid5
errors$GARCH <- forecasts$CPI_diff_log - forecasts$garch_m01_arma01_norm

# RNN vs ARIMA
dm.test(tail(errors$NN, n = n_out), tail(errors$ARIMA, n = n_out)) 

# RNN vs GARCH
dm.test(tail(errors$NN, n = n_out), tail(errors$GARCH, n = n_out)) 

# ARIMA vs GARCH
dm.test(tail(errors$ARIMA, n = n_out), tail(errors$GARCH, n = n_out)) 
