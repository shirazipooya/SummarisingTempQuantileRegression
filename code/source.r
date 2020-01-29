
# Requirement -------------------------------------------------------------

  # remove objects from a specified environment
  rm(list = ls())

  # shuts down all open graphics devices
  graphics.off()

  # load require functions
  source(file = "code/func.R")

  # load requaired packages
  wants = c("pacman", "dplyr", "tibble", "lubridate", "ggplot2",
            "ClimClass", "rgdal", "ggmap", "ggrepel", "quantreg",
            "cosinor", "pbapply", "parallel", "ggdendro", "dendextend")

  has = wants %in% rownames(x = installed.packages())

  if (any(!has)) install.packages(wants[!has])

  pacman::p_load(char = wants)

  # remove unwanted objects
  rm(has, wants)

# Select Required Stations ------------------------------------------------

  # load yearly data
  dataYearly <- readRDS(file = "data/dataYearly.RDS")

  # select stations (to: 2018, number of yers: 30)
  toYear = 2018
  nYear = 60
  fromYear = toYear - nYear + 1

  st <- dataYearly %>%
    select(Site, Date, t) %>%
    na.omit() %>%
    group_by(Site) %>%
    summarise(minYear = year(x = min(Date)),
              maxYear = year(x = max(Date))) %>%
    mutate(numYear = maxYear - minYear + 1) %>%
    filter(maxYear == toYear & numYear >= nYear) %>%
    select(Site) %>%
    unlist()

  names(st) <- NULL

  # remove unwanted objects
  rm(dataYearly)

# Load And Select Data ----------------------------------------------------

  dataDaily <- readRDS(file = "data/dataDaily.RDS") %>%
    filter(Site %in% st) %>%
    filter((year(x = Date) >= fromYear) & (year(x = Date) <= toYear)) %>%
    rename(station = Site, date = Date)

  # detecting missing values of t
  dataDailyMissing <- dataDaily %>%
    group_by(station) %>%
    summarise(missing = (round(x = 100 * sum(is.na(x = t)) / length(x = t),
                               digits = 2)))

  dataMonthly <- readRDS(file = "data/dataMonthly.RDS") %>%
    filter(Site %in% st) %>%
    filter((year(x = Date) >= fromYear) & (year(x = Date) <= toYear)) %>%
    rename(station = Site, date = Date)

  # remove unwanted objects
  rm(fromYear, toYear, nYear, st)

# Classification Of Climate According To Koeppen - Geiger -----------------

  dataKG <- dataMonthly %>%
    mutate(year = year(x = date),
           month = month(x = date)) %>%
    select(station, year, month, tmin, tmax, t, rrr) %>%
    na.omit()

  # dataKG convert to list
  dataKG <- split(x = dataKG %>% select(-station),
                  f = dataKG %>% select(station))

  # creates climate mean monthly values
  newNames <- c('year', 'month', 'Tn', 'Tx', 'Tm', 'P')

  climateNormal <- lapply(X = lapply(X = dataKG,
                                     FUN = setNames, newNames),
                          FUN = ClimClass::climate, max.perc.missing = 50)

  # calculate Koeppen - Geiger’s climate classification
  KGCC <- do.call(what = rbind.data.frame,
                  args = lapply(X = climateNormal,
                                FUN = ClimClass::koeppen_geiger,
                                      A_B_C_special_sub.classes = FALSE,
                                      clim.resume_verbose = TRUE,
                                      class.nr = FALSE)) %>%
    rownames_to_column("station") %>%
    select(station, class) %>%
    left_join(y = read.csv(file = 'data/stations.csv') %>%
                mutate(station = as.character(station)),
              by = "station") %>%
    select(station, name, latitude, longitude, elevation, class) %>%
    left_join(y = dataDailyMissing,
              by = "station")

  # export Koeppen - Geiger’s climate classification
  write.csv(x = KGCC,
            file = "result/Koeppen_Geigers_Climate_Classification.csv",
            row.names = FALSE)

  # remove unwanted objects
  rm(dataKG, newNames, climateNormal, dataDailyMissing)

# Plot Case Study ---------------------------------------------------------

  # load shapefile and converted to a dataframe:
  iran <- rgdal::readOGR(dsn = "data/shapefiles/Iran", layer = "Iran") %>%
    ggplot2::fortify()

  sea <- rgdal::readOGR(dsn = "data/shapefiles/Sea", layer = "Sea") %>%
    ggplot2::fortify()

  # plot
  iran_map <- ggplot() +
    theme_bw() +
    geom_polygon(data = iran,
                 mapping = aes(x = long, y = lat, group = group),
                 color = 'black', fill = 'white', size = .2) +
    geom_polygon(data = sea,
                 mapping = aes(x = long, y = lat, group = group),
                 color = 'steelblue1', fill = 'steelblue1', size = .2) +
    geom_point(data = KGCC,
               mapping = aes(x = longitude, y = latitude, pch = class),
               size = 3) +
    geom_text_repel(data = KGCC,
                    mapping = aes(x = longitude, y = latitude, label = station)) +
    coord_map() +
    xlab(label = "") +
    ylab(label = "") +
    theme(text = element_text(size = 18))

  # save plot
  ggsave(filename = "Case_Study_Map.png",
         plot = iran_map,
         width = 10,
         height = 10,
         units = "in",
         dpi = 1200,
         path = "./result/")

  # remove unwanted objects
  rm(iran, sea, iran_map)

# Quantile slopes ---------------------------------------------------------

  # data preparation
  data <- dataDaily %>%
    filter(!(month(date) == 2 & day(date) == 29)) %>% # remove day */02/29 from data
    split(f = as.factor(.$station)) # convert data to list base station column

  # initiate cluster
  cl <- makeCluster(detectCores())

  # apply operations using clusters
  clusterEvalQ(cl = cl,
               expr = {
                 library(dplyr)
                 library(lubridate)
                 library(quantreg)
                 library(cosinor)
               })

  clusterExport(cl = cl,
                c('calculateQR', 'data'))

  # calculate QR
  QR_TS <- parLapply(cl = cl,
                     X = data,
                     fun = function(x) calculateQR(data = x,
                                                   x = 't'))

  stopCluster(cl)

  # convert result to data.frame
  slope <- names(x = QR_TS[[1]])
  nStation <- length(QR_TS)
  result <- data.frame(station = names(QR_TS))

  for (i in slope) {
    result[paste0("Q", i)] <- unlist(x = lapply(X = names(QR_TS), FUN = function(x) round(x = QR_TS[[x]][[i]][1], digits = 3)))
    result[paste0("SE", i)] <- unlist(x = lapply(X = names(QR_TS), FUN = function(x) round(x = QR_TS[[x]][[i]][2], digits = 3)))
    result[paste0("PV", i)] <- unlist(x = lapply(X = names(QR_TS), FUN = function(x) round(x = QR_TS[[x]][[i]][3], digits = 3)))
  }

  result <- left_join(x = KGCC,
                      y = result,
                      by = 'station')

  rownames(result) <- result$station

  write.csv(x = result,
            file = "result/Quantile_Slopes_Table.csv",
            row.names = FALSE,
            col.names = TRUE)

  # remove unwanted objects
  rm(i, nStation, slope, cl, QR_TS)

# Time Series Plot + Quantile Slopes Plot ---------------------------------

# initiate cluster
invisible(utils::memory.limit(16384))
cl <- makeCluster(detectCores())

# apply operations using clusters
clusterEvalQ(cl = cl,
             expr = {
               library(dplyr)
               library(lubridate)
               library(quantreg)
               library(cosinor)
               library(ggplot2)
             })

clusterExport(cl = cl,
              c('QuantileSlopesPlot', 'data'))

# plot
st <- c("40718", "40736", "40745", "40856")

for (i in st) {
  lapply(X = data[i],
         FUN = function(x) plotTS(data = x,
                                  x = 't',
                                  x_label = 'Temperature (°C)',
                                  y_limit = c(-5, 35)))

  parLapply(cl = cl,
            X = data[i],
            fun = function(x) QuantileSlopesPlot(data = x,
                                                 x = 't',
                                                 x_label = 'Temperature (°C)',
                                                 y_limit = c(0.15,0.55)))

lapply(X = data[i],
       FUN = function(x) Bootstrap_Distribution_QR(data = x,
                                                   x = 't',
                                                   x_limit = c(0.0, 0.8),
                                                   x_label = 'Slope (°C/decade)'))
}

stopCluster(cl)

# Dendrogram --------------------------------------------------------------

for (i in c("Q0.05", "Q0.5", "Q0.95", "Qols")) {
  # open jpeg file
  png(filename = paste0("./result/Dendrogram/", i, "_Dendrogram.png"),
      width = 8, height = 6, units = "in", res = 1200)

  # compute distances and hierarchical clustering
  result[c(i)] %>%
    dist(method = "euclidean") %>%
    hclust(method = "average") %>%
    as.dendrogram %>%
    set("leaves_pch", 19) %>%
    set("branches_k_color", value = c("red", "green", "blue", "orange"), k = 4) %>%
    set("branches_lwd", 2) %>%
    plot(horiz = FALSE)

  # close the file
  dev.off()
}
