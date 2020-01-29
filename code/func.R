
# Plot Time Series ----------------------

plotTS <- function(data,
                   x,
                   x_label,
                   y_limit)
{
  # removed seasonality from each temperature record via sinusoidal regression
  d <- data %>%
    select_('date', x) %>%
    mutate(JDay = rep(x = 1:365, nrow(data) / 365),
           X = 1,
           time = time(date))

  fit <- cosinor.lm(formula = as.formula(paste0(x, ' ~ time(JDay) + X + amp.acro(X)')),
                    data = as.data.frame(d),
                    period = 365,
                    na.action = na.omit)

  d$JDay <- as.factor(d$JDay)
  d$time <- as.factor(d$time)
  fit_value <- d %>%
    dplyr::select(JDay, time) %>%
    left_join(y = data.frame(time = as.factor(names(fit$fit$fitted.values)),
                             f_v = fit$fit$fitted.values),
              by = 'time') %>%
    group_by(JDay) %>%
    summarise(fit_value = mean(x = f_v, na.rm = TRUE))

  # seasonality is removed from each temperature record (x_sr = seasonality removed)
  d$x_sr = d[[x]] -
    rep(x = fit_value$fit_value, nrow(d)/365) +
    rep(x = fit$coefficients['(Intercept)'], nrow(d))

  # plot QR time serise
  P <- ggplot(data = d) +
    theme_bw() +
    geom_line(mapping = aes(x = date, y = x_sr),
              size = 0.6,
              col = 'gray55') +
    geom_quantile(mapping = aes(x = date, y = x_sr),
                  quantiles = c(0.05),
                  colour = "black",
                  linetype = "longdash",
                  size = 0.8) +
    geom_quantile(mapping = aes(x = date, y = x_sr),
                  quantiles = c(0.95),
                  colour = "black",
                  linetype = "dotdash",
                  size = 0.8) +
    geom_quantile(mapping = aes(x = date, y = x_sr),
                  quantiles = c(0.50),
                  colour = "black",
                  linetype = "solid",
                  size = 0.8) +
    xlab(label = "") +
    ylab(label = x_label) +
    ylim(y_limit) +
    theme(text = element_text(size = 12))

  # save plot QR time serise
  st <- unique(data$station)
  ggsave(filename = paste(st, ".png", sep = ""),
         plot = P,
         width = 6,
         height = 3,
         units = "in",
         dpi = 1200,
         path = "./result/Time Serise Plot/")
}

# Plot Quantile Slopes ----------------------------------------------------

QuantileSlopesPlot <- function(data,
                               x,
                               x_label,
                               y_limit)
{
  from = 0.05
  nfrom = 1

  to = 0.95
  nto = 91

  by = 0.01
  nby = 1

  # removed seasonality from each temperature record via sinusoidal regression
  d <- data %>%
    select_('date', x) %>%
    mutate(JDay = rep(x = 1:365, nrow(data) / 365),
           X = 1,
           time = time(date))

  fit <- cosinor.lm(formula = as.formula(paste0(x, ' ~ time(JDay) + X + amp.acro(X)')),
                    data = as.data.frame(d),
                    period = 365,
                    na.action = na.omit)

  d$JDay <- as.factor(d$JDay)
  d$time <- as.factor(d$time)
  fit_value <- d %>%
    select(JDay, time) %>%
    left_join(y = data.frame(time = as.factor(names(fit$fit$fitted.values)),
                             f_v = fit$fit$fitted.values),
              by = 'time') %>%
    group_by(JDay) %>%
    summarise(fit_value = mean(x = f_v, na.rm = TRUE))

  # seasonality is removed from each temperature record (x_sr = seasonality removed)
  d$x_sr = d[[x]] -
    rep(x = fit_value$fit_value, nrow(d)/365) +
    rep(x = fit$coefficients['(Intercept)'], nrow(d))

  # calculate QR coef.
  QR <- summary(object = rq(formula = x_sr ~ date,
                            tau = seq(from = from, to = to, by = by),
                            data = d,
                            ci = FALSE),
                se = 'boot',
                covariance = TRUE,
                R = 200,
                bsmethod = "mcmb")

  # calculate OLS coef.
  OLS <- summary(lm(formula = x_sr ~ date, data = d))

  # plot type 01:
  # extract some quantile slop and standard error (0.05, 0.08, 0.11, ..., 0.95).
  coef.se.vector <- lapply(X = seq(from = nfrom, to = nto, by = nby),
                           FUN = function(x) data.frame(coef = QR[[x]][["coefficients"]][2,1] * 10 * 365,
                                                        se = QR[[x]][["coefficients"]][2,2] * 10 * 365)) %>%
    unlist()

  # convert to data frame and add tau quantile.
  coef.se.dataframe <- data.frame(coef = coef.se.vector[names(coef.se.vector) == "coef"],
                                  se = coef.se.vector[names(coef.se.vector) == "se"])
  coef.se.dataframe$tau <- seq(from = from, to = to, by = by)

  # plot quantile slopes and corresponding standard errors and ordinary least squares slope.
  P1 <- ggplot(data = coef.se.dataframe, mapping = aes(x = tau, y = coef)) +
    geom_hline(mapping = aes(yintercept = OLS[["coefficients"]][2,1] * 10 * 365),
               linetype = "dashed",
               color = "red",
               size = 0.5) +
    geom_pointrange(mapping = aes(ymin = coef - se, ymax = coef + se)) +
    xlab(label = "quantile") +
    ylab(label = "°C/decade") +
    xlim(0,1) +
    ylim(y_limit)

  # save plot
  st <- unique(data$station)
  ggsave(filename = paste(st, "_type01.png", sep = ""),
         plot = P1,
         width = 8,
         height = 4,
         units = "in",
         dpi = 1200,
         path = "./result/Quantile Slopes Plot/")

  # plot type 02
  # extract some quantile slop and standard error (0.05, 0.06, 0.07, ..., 0.95).
  con.vector <- lapply(X = seq(from = nfrom, to = nto, by = nby),
                       FUN = function(x) data.frame(coef = QR[[x]][["coefficients"]][2,1] * 10 * 365,
                                                     Low = t(apply(QR[[x]][["B"]], 2, quantile, c(0.025,0.975)))[2,1] * 10 * 365,
                                                     High = t(apply(QR[[x]][["B"]], 2, quantile, c(0.025,0.975)))[2,2] * 10 * 365)) %>%
    unlist()

  # convert to data frame and add tau quantile.
  con.dataframe <- data.frame(coef = con.vector[names(con.vector) == "coef"],
                              Low = con.vector[names(con.vector) == "Low"],
                              High = con.vector[names(con.vector) == "High"])
  con.dataframe$tau <- seq(from = from, to = to, by = by)

  # plot quantile slopes and corresponding standard errors and ordinary least squares slope.
  data.plot <- data.frame(x = c(con.dataframe$tau,rev(con.dataframe$tau)),
                          y = c(con.dataframe$Low,rev(con.dataframe$High)))
  P2 <- ggplot(data = data.plot) +
    geom_polygon(mapping = aes(x = x,
                               y = y), fill = "gray55") +
    geom_hline(mapping = aes(yintercept = OLS[["coefficients"]][2,1] * 10 * 365),
               linetype = "dashed",
               color = "red",
               size = 0.5) +
    geom_point(data = con.dataframe, mapping = aes(x = tau, y = coef), color = "black") +
    xlab(label = "quantile") +
    ylab(label = "°C/decade") +
    xlim(0,1) +
    ylim(y_limit)

  # save plot
  st <- unique(data$station)
  ggsave(filename = paste(st, "_type02.png", sep = ""),
         plot = P2,
         width = 8,
         height = 4,
         units = "in",
         dpi = 1200,
         path = "./result/Quantile Slopes Plot/")
}

# Calculate Quantile Slopes ---------------------------------------------------------

calculateQR <- function(data,
                        x)
{
  # removed seasonality from each temperature record via sinusoidal regression
  d <- data %>%
    select_('date', x) %>%
    mutate(JDay = rep(x = 1:365, nrow(data) / 365),
           X = 1,
           time = time(date))

  fit <- cosinor.lm(formula = as.formula(paste0(x, ' ~ time(JDay) + X + amp.acro(X)')),
                    data = as.data.frame(d),
                    period = 365,
                    na.action = na.omit)

  d$JDay <- as.factor(d$JDay)
  d$time <- as.factor(d$time)
  fit_value <- d %>%
    select(JDay, time) %>%
    left_join(y = data.frame(time = as.factor(names(fit$fit$fitted.values)),
                             f_v = fit$fit$fitted.values),
              by = 'time') %>%
    group_by(JDay) %>%
    summarise(fit_value = mean(x = f_v, na.rm = TRUE))

  # seasonality is removed from each temperature record (x_sr = seasonality removed)
  d$x_sr = d[[x]] -
    rep(x = fit_value$fit_value, nrow(d)/365) +
    rep(x = fit$coefficients['(Intercept)'], nrow(d))

  # calculate QR coef.
  QR <- summary(object = rq(formula = x_sr ~ date,
                            tau = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
                            data = d,
                            ci = FALSE),
                se = 'boot',
                covariance = TRUE,
                R = 200,
                bsmethod = "mcmb")

  # calculate OLS coef.
  OLS <- summary(lm(formula = x_sr ~ date, data = d))

  # Quantile Slopes Table
  p.vector <- lapply(X = 1:7,
                     FUN = function(x) data.frame(coef = QR[[x]][["coefficients"]][2,1] * 10 * 365,
                                                  se = QR[[x]][["coefficients"]][2,2] * 10 * 365,
                                                  pv = QR[[x]][["coefficients"]][2,4])) %>%
    unlist()

  # convert to data frame and add row names.
  p.dataframe <- data.frame(coef = p.vector[names(p.vector) == "coef"],
                            se = p.vector[names(p.vector) == "se"],
                            pv = p.vector[names(p.vector) == "pv"])
  rownames(p.dataframe) <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
  p.dataframe <- as.data.frame(t(p.dataframe))

  # extract slops and standard error and statistically signiﬁcant for ordinary least squares
  p.dataframe$`ols`[1] <- OLS[["coefficients"]][2,1] * 10 * 365
  p.dataframe$`ols`[2] <- OLS[["coefficients"]][2,2] * 10 * 365
  p.dataframe$`ols`[3] <- OLS[["coefficients"]][2,4]

  # store results
  return(p.dataframe)
}

# Bootstrap Distribution Of Quantile Slopes -------------------------------

Bootstrap_Distribution_QR <- function(data,
                                      x,
                                      x_label,
                                      x_limit)
{
  # removed seasonality from each temperature record via sinusoidal regression
  d <- data %>%
    select_('date', x) %>%
    mutate(JDay = rep(x = 1:365, nrow(data) / 365),
           X = 1,
           time = time(date))

  st <- unique(data$station)

  fit <- cosinor.lm(formula = as.formula(paste0(x, ' ~ time(JDay) + X + amp.acro(X)')),
                    data = as.data.frame(d),
                    period = 365,
                    na.action = na.omit)

  d$JDay <- as.factor(d$JDay)
  d$time <- as.factor(d$time)
  fit_value <- d %>%
    select(JDay, time) %>%
    left_join(y = data.frame(time = as.factor(names(fit$fit$fitted.values)),
                             f_v = fit$fit$fitted.values),
              by = 'time') %>%
    group_by(JDay) %>%
    summarise(fit_value = mean(x = f_v, na.rm = TRUE))

  # seasonality is removed from each temperature record (x_sr = seasonality removed)
  d$x_sr = d[[x]] -
    rep(x = fit_value$fit_value, nrow(d)/365) +
    rep(x = fit$coefficients['(Intercept)'], nrow(d))

  # calculate QR coef.
  QR <- summary(object = rq(formula = x_sr ~ date,
                            tau = c(0.05, 0.50, 0.95),
                            data = d,
                            ci = FALSE),
                se = 'boot',
                covariance = TRUE,
                R = 200,
                bsmethod = "mcmb")

  # Bootstrap Distribution Quantile Slopes Plot
  for (t in 1:3)
  {
    data <- as.data.frame(QR[[t]][["B"]][,2] * 10 * 365)
    names(data) <- "x"

    P <- ggplot() +
      geom_histogram(data = data, mapping = aes(x = x), binwidth = 0.04) +
      geom_vline(mapping = aes(xintercept = QR[[t]][["coefficients"]][2,1] * 10 * 365),
                 linetype = "dashed",
                 color = "red",
                 size = 1) +
      xlim(x_limit) +
      ylim(c(0, 200)) +
      xlab(label = x_label) +
      ylab(label = "Frequency")

    ggsave(filename = paste(st, "_", c(0.05, 0.50, 0.95)[t], ".png", sep = ""),
           plot = P,
           width = 6,
           height = 3,
           units = "in",
           dpi = 1200,
           path = "./result/Bootstrap Distribution Quantile Slopes Plot/")
  }
}
