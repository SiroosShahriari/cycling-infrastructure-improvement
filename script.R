library("GWmodel")
library("tidyverse")
library("sf")
library("mapview")
library("lubridate")
library("reshape")
library("fastDummies")
library("geosphere")
library("forecast")
library("dplyr")
library("ggplot2")
library("dodgr")
library("spatialreg")
library("janitor")
library("splm")
library("plm")
library("geodist")
library("plotly")


# Lyon -----------------

# read data
l.data <- read.csv(here::here("data/Lyon/input-data.csv"), header = TRUE)
l.description <- read.csv(here::here("data/Lyon/input-description.csv"))

# data cleaning
## group by genre
# monthly aggregate data
l.data.m <- l.data[, 3:127]
l.data.m$time.m <- dmy(l.data.m$time.m)
l.monthly.agg <- aggregate(l.data.m[-1],
  by = list(l.data.m$time.m),
  FUN = mean
)


# find duplicates in data
dups <- which(duplicated(l.data))


# select data from March 2018 till March 2020 (pre-Covid era)
l.monthly.agg <- l.monthly.agg[1:25, ]

# remove columns with 0 in the first
ones <- which(l.monthly.agg == 0, arr.ind = T)
l.monthly.agg <- l.monthly.agg[-c(unique(ones[, 2]))]

# remove columns with na
na <- (which(colSums(is.na(l.monthly.agg)) > 0))
l.monthly.agg <- l.monthly.agg[-c(unique(na))]


# standardize panel data
df <- l.monthly.agg %>% mutate_all(~ (scale(.) %>% as.vector()))
df[, 1] <- l.monthly.agg[, 1]
l.monthly.agg <- df


# remove duplicates
l.monthly.agg <- l.monthly.agg[!duplicated(as.list(l.monthly.agg))]

# unpivot data
l.monthly.agg.unpiv <- melt(l.monthly.agg, id = c("Group.1"))

# add improvement effect
l.monthly.agg.unpiv$improvement <- 0

# Quai Claude Bernard    31/12/2018  Added physically segregated cycle lane
l.monthly.agg.unpiv$improvement <- l.monthly.agg.unpiv$improvement +
  ifelse(l.monthly.agg.unpiv$Group.1 >= "2019-01-1" &
    l.monthly.agg.unpiv$variable == "X102027497", 1, 0)

l.monthly.agg.unpiv$improvement <- l.monthly.agg.unpiv$improvement +
  ifelse(l.monthly.agg.unpiv$Group.1 >= "2019-01-1" &
    l.monthly.agg.unpiv$variable == "X101027497", 1, 0)


# add time dummy variable
l.monthly.agg.unpiv <- dummy_cols(l.monthly.agg.unpiv, select_columns = "Group.1")
l.monthly.agg.unpiv <- dummy_cols(l.monthly.agg.unpiv, select_columns = "variable")


#| echo: false
lyon_streetnet <- dodgr_streetnet("Lyon, France")


l.bike_usage <- l.monthly.agg.unpiv %>%
  janitor::clean_names() %>%
  dplyr::rename_with(~ gsub("group_1", "date", .x)) %>%
  dplyr::rename_with(~ gsub("^variable_x", "counter_", .x)) %>%
  dplyr::rename(channel_id = "variable") %>%
  mutate(channel_id = gsub("^X", "", channel_id))
head(l.bike_usage)




l.counters_with_physical_improvement <- c(102027497, 101027497)
l.counters_with_paint_improvement <- c()
l.counters_with_shared_lane_improvement <- c()



l.bike_usage <- l.bike_usage %>%
  mutate(
    improvement_physical = ifelse(channel_id %in% l.counters_with_physical_improvement & improvement == 1, 1, 0),
    improvement_paint = ifelse(channel_id %in% l.counters_with_paint_improvement & improvement == 1, 1, 0),
    improvement_shared_lane = ifelse(channel_id %in% l.counters_with_shared_lane_improvement & improvement == 1, 1, 0)
  )

l.bike_usage <- l.bike_usage %>% filter(!channel_id %in% c(101026725, 102032824, 102055357, 101055357, 102049670, 101049670, 104042249, 353309131, 102030274, 101030274, 101030269, 353264965, 101063321, 101026724, 102047694, 353264065, 102030102, 101030102))



l.channels <- read.csv(here::here("data/Lyon/channels.csv")) %>% janitor::clean_names()
l.counters_sf <- read.csv(here::here("data/Lyon/sites.csv")) %>%
  select(site_id, lon, lat) %>%
  merge(select(l.channels, channel_id, site_id), by = "site_id") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  filter(channel_id %in% l.bike_usage$channel_id)
l.counters_sf


l.lyon_streetnet <- dodgr_streetnet(pts = st_coordinates(l.counters_sf))


# Match points to graph
l.graph <- weight_streetnet(l.lyon_streetnet, wt_profile = "bicycle")
verts <- dodgr_vertices(l.graph)
pts <- match_points_to_graph(l.graph, l.counters_sf)


# Shortest distance between counters
l.d <- dodgr_dists(l.graph,
  from = l.graph[pts, ]$from_id,
  to = l.graph[pts, ]$from_id
)

# linear with x km cut-off
# @param cutoff default to 1000m.
compute_linear_weight <- function(l.d, cutoff) {
  stopifnot(cutoff > 0)
  for (row in seq_len(nrow(l.d))) {
    d_row_cutoff <- ifelse(l.d[row, ] <= cutoff, l.d[row, ], NA)
    d_row_cutoff <- (cutoff - d_row_cutoff) / cutoff
    d_row_cutoff <- ifelse(is.na(d_row_cutoff), 0, d_row_cutoff)
    l.d[row, ] <- d_row_cutoff
  }
  l.d
}

l.counters_spatial_weights <- compute_linear_weight(l.d, cutoff = 1000)

# convert distance to kilometers
l.dd <- l.d / 1000


# different weighing systems
# counters_spatial_weights <- GWmodel::gw.weight(dd,bw=24.99574,kernel="exponential",adaptive=FALSE)
l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = max(l.dd), kernel = "gaussian", adaptive = FALSE)
# l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = max(l.dd), kernel = "bisquare", adaptive = FALSE)
# l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = max(l.dd), kernel = "tricube", adaptive = FALSE)
# l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = max(l.dd), kernel = "boxcar", adaptive = FALSE)

# Fit Spatial Panel Linear Models
l.results <- matrix(0, 16, 2)
l.charlist <- c("2019-02-01", "2019-03-01", "2019-04-01", "2019-05-01", "2019-06-01", "2019-07-01", "2019-08-01", "2019-09-01", "2019-10-01", "2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01")

for (i in l.charlist)
{
  i <- c("2024-02-01")
  l.bike_usage_panel_df <- l.bike_usage %>%
    filter(date < i) %>%
    select(channel_id, date, value, starts_with("improvement")) %>%
    plm::pdata.frame(index = c("channel_id", "date"))

  l.sdid_ind <-
    splm::spml(
      value ~ improvement_physical + factor(date),
      data = l.bike_usage_panel_df,
      listw = spdep::mat2listw(l.counters_spatial_weights),
      model = "within",
      spatial.error = "b",
      Hess = FALSE,
      lag = TRUE,
      effect = "individual"
    )

  summary(l.sdid_ind)
  l.results[which(l.charlist == i), 1:2] <- l.sdid_ind$coefficients[3:4]
}

summary(l.sdid_ind)




#############################################################################################################################################
### PARIS
# read data
p.data <- read.csv(here::here("data/Paris/input-data.csv"), header = TRUE)
p.description <- read.csv(here::here("data/Paris/input-description.csv"))


# data cleaning
## group by genre
# monthly aggregate data
p.data.m <- p.data[, 3:127]
p.data.m$time.m <- dmy(p.data.m$time.m)
p.monthly.agg <- aggregate(p.data.m[-1],
  by = list(p.data.m$time.m),
  FUN = mean
)

# find duplicates in data
dups <- which(duplicated(p.data))


# monthly.agg <- weekday.peakh.monthly.agg

# select data from March 2018 till March 2020 (pre-Covid era)
# monthly.agg <- monthly.agg[7:47,]

p.monthly.agg <- p.monthly.agg[6:34, ]

# remove columns with 0 in the first
ones <- which(p.monthly.agg == 0, arr.ind = T)
p.monthly.agg <- p.monthly.agg[-c(unique(ones[, 2]))]

# remove columns with na
na <- (which(colSums(is.na(p.monthly.agg)) > 0))
p.monthly.agg <- p.monthly.agg[-c(unique(na))]


# standardize panel data
df4 <- p.monthly.agg %>% mutate_all(~ (scale(.) %>% as.vector()))
df4[, 1] <- p.monthly.agg[, 1]
p.monthly.agg <- df4


# remove duplicates
p.monthly.agg <- p.monthly.agg[!duplicated(as.list(p.monthly.agg))]
# unpivot data
p.monthly.agg.unpiv <- melt(p.monthly.agg, id = c("Group.1"))


# add improvement effect
p.monthly.agg.unpiv$improvement <- 0

# boulevard Voltaire   31/08/2019    Added physically segregated cycle lane
p.monthly.agg.unpiv$improvement <- p.monthly.agg.unpiv$improvement +
  ifelse(p.monthly.agg.unpiv$Group.1 >= "2019-09-1" &
    p.monthly.agg.unpiv$variable == "X100044506.SC", 1, 0)


# add time dummy variable
p.monthly.agg.unpiv <- dummy_cols(p.monthly.agg.unpiv, select_columns = "Group.1")
p.monthly.agg.unpiv <- dummy_cols(p.monthly.agg.unpiv, select_columns = "variable")



# add inv-distance as reggressor

p.data.lon.lat <- p.monthly.agg.unpiv
p.description.ll <- p.description
p.description.ll$id <- sub("^", "X", p.description.ll$id.full)
p.description.ll$id <- sub("-", ".", p.description.ll$id)
colnames(p.description.ll)[4] <- "variable"
p.description.ll <- p.description.ll[-1]

p.description.ll <- p.description.ll[!duplicated(p.description.ll[4]), ]

p.data.lon.lat <- inner_join(p.data.lon.lat, p.description.ll)

p.data.lon.lat <- p.data.lon.lat[!duplicated(p.data.lon.lat), ]

p.data.lon.lat$s.improvement <- 0
p.data.lon.lat$s.improvement <- as.numeric(p.data.lon.lat$s.improvement)

not.imp.li <- which(p.data.lon.lat$improvement != 1)
imp.li <- which(p.data.lon.lat$improvement == 1)


p.monthly.agg.unpiv$s.improvement <- 0




# relocate the s.improvement
p.data.lon.lat.1 <- p.data.lon.lat %>% relocate(s.improvement, .after = improvement)
p.monthly.agg.unpiv <- p.data.lon.lat.1
p.monthly.agg.unpiv <- p.monthly.agg.unpiv[-ncol(p.monthly.agg.unpiv)]

# relocate the s.improvement
p.monthly.agg.unpiv <- p.monthly.agg.unpiv %>% relocate(s.improvement, .after = improvement)


paris_streetnet <- dodgr_streetnet("Paris, France")

p.monthly.agg.unpiv$variable <- sub(".SC", "", p.monthly.agg.unpiv$variable)

p.bike_usage <- p.monthly.agg.unpiv %>%
  janitor::clean_names() %>%
  dplyr::rename_with(~ gsub("group_1", "date", .x)) %>%
  dplyr::rename_with(~ gsub("^variable_x", "counter_", .x)) %>%
  dplyr::rename(channel_id = "variable") %>%
  mutate(channel_id = gsub("^X", "", channel_id))
head(p.bike_usage)


p.counters_with_physical_improvement <- c(100044506)
p.counters_with_paint_improvement <- c()
p.counters_with_shared_lane_improvement <- c()

p.bike_usage <- p.bike_usage %>%
  mutate(
    improvement_physical = ifelse(channel_id %in% p.counters_with_physical_improvement & improvement == 1, 1, 0),
    improvement_paint = ifelse(channel_id %in% p.counters_with_paint_improvement & improvement == 1, 1, 0),
    improvement_shared_lane = ifelse(channel_id %in% p.counters_with_shared_lane_improvement & improvement == 1, 1, 0)
  )

# p.bike_usage <- p.bike_usage %>% filter(! channel_id %in% c(100041488))

p.counters_sf <- p.description.ll %>%
  select(variable, lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
p.counters_sf$variable <- sub("X", "", p.counters_sf$variable)
p.counters_sf$variable <- sub(".SC", "", p.counters_sf$variable)

p.counters_sf <- p.counters_sf %>% filter(variable %in% p.bike_usage$channel_id)
p.counters_sf
colnames(p.counters_sf)[1] <- "channel_id"

paris_streetnet <- dodgr_streetnet(pts = st_coordinates(p.counters_sf))





# Match points to graph
p.graph <- weight_streetnet(paris_streetnet, wt_profile = "bicycle")
p.verts <- dodgr_vertices(p.graph)
p.pts <- match_points_to_graph(p.graph, p.counters_sf)

# Shortest distance between counters
p.d <- dodgr_dists(p.graph,
  from = p.graph[p.pts, ]$from_id,
  to = p.graph[p.pts, ]$from_id
)

#  Bandwidth selection for basic GWR
p.dd <- p.d / 1000


# counters_spatial_weights <- GWmodel::gw.weight(dd,bw=11.37453,kernel="exponential",adaptive=FALSE)
p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = max(p.dd), kernel = "gaussian", adaptive = FALSE)
# p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = max(p.dd), kernel = "bisquare", adaptive = FALSE)
# p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = max(p.dd), kernel = "tricube", adaptive = FALSE)
# p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = max(p.dd), kernel = "boxcar", adaptive = FALSE)

# Fit Spatial Panel Linear Models
p.results <- matrix(0, 16, 2)
p.charlist <- c("2019-02-01", "2019-03-01", "2019-04-01", "2019-05-01", "2019-06-01", "2019-07-01", "2019-08-01", "2019-09-01", "2019-10-01", "2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01")
# charlist <- c("2018-09-01","2018-10-01","2018-11-01","2018-12-01","2019-01-01","2019-02-01","2019-03-01","2019-04-01","2019-05-01","2019-06-01","2019-07-01","2019-08-01","2019-09-01","2019-10-01","2019-11-01","2019-12-01","2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01")

for (i in p.charlist)
{
  i <- c("2024-02-01")
  p.bike_usage_panel_df <- p.bike_usage %>%
    filter(date < i) %>%
    select(channel_id, date, value, starts_with("improvement")) %>%
    plm::pdata.frame(index = c("channel_id", "date"))

  p.sdid_ind <-
    splm::spml(
      value ~ improvement_physical + factor(date),
      data = p.bike_usage_panel_df,
      listw = spdep::mat2listw(p.counters_spatial_weights),
      model = "within",
      spatial.error = "b",
      Hess = FALSE,
      lag = TRUE,
      effect = "individual"
    )

  summary(p.sdid_ind)
  p.results[which(p.charlist == i), 1:2] <- p.sdid_ind$coefficients[3:4]
}

summary(p.sdid_ind)


##############################################################################
# checking bandwidth

i <- c("2024-02-01")
l.bike_usage_panel_df <- l.bike_usage %>%
  filter(date < i) %>%
  select(channel_id, date, value, starts_with("improvement")) %>%
  plm::pdata.frame(index = c("channel_id", "date"))

i <- c("2024-02-01")
p.bike_usage_panel_df <- p.bike_usage %>%
  filter(date < i) %>%
  select(channel_id, date, value, starts_with("improvement")) %>%
  plm::pdata.frame(index = c("channel_id", "date"))

##########
# compute_linear_weight <- function(d, cutoff) {
#   stopifnot(cutoff > 0)
#   for (row in seq_len(nrow(d))) {
#     d_row_cutoff <- ifelse(d[row, ] <= cutoff, d[row, ], NA)
#     d_row_cutoff <- (cutoff - d_row_cutoff) / cutoff
#     d_row_cutoff <- ifelse(is.na(d_row_cutoff), 0, d_row_cutoff)
#     d[row, ] <- d_row_cutoff
#   }
#   d
# }


##########
# checking different bandwith
bandwidth.results.l.sigma2 <- matrix(0, 5, 5)
bandwidth.results.l.rmse <- matrix(0, 5, 5)

bandwidth.results.p.sigma2 <- matrix(0, 5, 5)
bandwidth.results.p.rmse <- matrix(0, 5, 5)

for (bandwidth in 1:4) {
  bandwidthcounter <- bandwidth
  listt <- c(1, 2.5, 5, 10)
  bandwidth <- listt[bandwidth]

  for (bandwidth.t in 1:3)
  {
    if (bandwidth.t == 1) {
      l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = bandwidth, kernel = "exponential", adaptive = FALSE)
    }
    if (bandwidth.t == 2) {
      l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = bandwidth, kernel = "gaussian", adaptive = FALSE)
    }
    if (bandwidth.t == 3) {
      listt <- c(5, 10, 15, 20)
      linearbandwidth <- listt[bandwidthcounter]
      l.counters_spatial_weights <- compute_linear_weight(l.dd, cutoff = linearbandwidth)
    }



    if (bandwidth.t == 1) {
      p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = bandwidth, kernel = "exponential", adaptive = FALSE)
    }
    if (bandwidth.t == 2) {
      p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = bandwidth, kernel = "gaussian", adaptive = FALSE)
    }
    if (bandwidth.t == 3) {
      listt <- c(5, 10, 15, 20)
      linearbandwidth <- listt[bandwidthcounter]
      p.counters_spatial_weights <- compute_linear_weight(p.dd, cutoff = linearbandwidth)
    }




    diag(l.counters_spatial_weights) <- 0
    diag(p.counters_spatial_weights) <- 0

    # l.counters_spatial_weights <-l.counters_spatial_weights/rowSums(l.counters_spatial_weights)
    # Fit Spatial Panel Linear Models

    l.sdid_ind <-
      splm::spml(
        value ~ improvement_physical + factor(date),
        data = l.bike_usage_panel_df,
        listw = spdep::mat2listw(l.counters_spatial_weights),
        model = "within",
        spatial.error = "b",
        Hess = FALSE,
        lag = TRUE,
        effect = "individual"
      )

    summary(l.sdid_ind)

    actual <- as.numeric(l.sdid_ind$model[, 1])

    # reference  :  file:///C:/Users/z5217441/Downloads/Spatial_Panel_Data_Models.pdf
    rss <- sum((l.sdid_ind$residuals)^2) ## residual sum of squares
    tss <- sum((actual - mean(actual))^2) ## total sum of squares
    rsq.l <- 1 - rss / tss
    rmse.l <- sqrt(mean((l.sdid_ind$residuals)^2))







    p.sdid_ind <-
      splm::spml(
        value ~ improvement_physical + factor(date),
        data = p.bike_usage_panel_df,
        listw = spdep::mat2listw(p.counters_spatial_weights),
        model = "within",
        spatial.error = "b",
        Hess = FALSE,
        lag = TRUE,
        effect = "individual"
      )

    summary(p.sdid_ind)



    actual <- as.numeric(p.sdid_ind$model[, 1])

    rss <- sum((p.sdid_ind$residuals)^2) ## residual sum of squares
    tss <- sum((actual - mean(actual))^2) ## total sum of squares
    rsq.p <- 1 - rss / tss
    rmse.p <- sqrt(mean((p.sdid_ind$residuals)^2))



    bandwidth.results.l.sigma2[bandwidth.t, bandwidthcounter] <- l.sdid_ind$sigma2
    bandwidth.results.l.rmse[bandwidth.t, bandwidthcounter] <- rmse.l


    bandwidth.results.p.sigma2[bandwidth.t, bandwidthcounter] <- p.sdid_ind$sigma2
    bandwidth.results.p.rmse[bandwidth.t, bandwidthcounter] <- rmse.p
  }
}



# building the model without the spatial impact
sdid.withoutspatialterm <-
  plm::plm(
    value ~ improvement_physical + factor(date),
    data = p.bike_usage_panel_df,
    model = "within",
    spatial.error = "b",
    Hess = FALSE,
    lag = TRUE,
    effect = "twoways"
  )


summary(sdid.withoutspatialterm)

actual <- as.numeric(sdid.withoutspatialterm$model[, 1])
rmse.without <- sqrt(mean((sdid.withoutspatialterm$residuals)^2))
sigma2.without <- var(sdid.withoutspatialterm$residuals)
#########################################################################################
# checking time to get significant out come



significance.results <- matrix(0, 20, 2)

l.charlist <- c("2019-01-01", "2019-02-01", "2019-03-01", "2019-04-01", "2019-05-01", "2019-06-01", "2019-07-01", "2019-08-01", "2019-09-01", "2019-10-01", "2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01")
p.charlist <- c("2019-09-01", "2019-10-01", "2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-01", "2020-08-01", "2020-09-01", "2020-10-01", "2020-11-01")

for (i in 1:15)
{
  l.i <- l.charlist[i]
  p.i <- p.charlist[i]

  l.bike_usage_panel_df <- l.bike_usage %>%
    filter(date <= l.i) %>%
    select(channel_id, date, value, starts_with("improvement")) %>%
    plm::pdata.frame(index = c("channel_id", "date"))

  p.bike_usage_panel_df <- p.bike_usage %>%
    filter(date <= p.i) %>%
    select(channel_id, date, value, starts_with("improvement")) %>%
    plm::pdata.frame(index = c("channel_id", "date"))




  l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = 10, kernel = "gaussian", adaptive = FALSE)
  p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = 10, kernel = "gaussian", adaptive = FALSE)


  diag(l.counters_spatial_weights) <- 0
  diag(p.counters_spatial_weights) <- 0

  # l.counters_spatial_weights <-l.counters_spatial_weights/rowSums(l.counters_spatial_weights)
  # Fit Spatial Panel Linear Models

  l.sdid_ind <-
    splm::spml(
      value ~ improvement_physical + factor(date),
      data = l.bike_usage_panel_df,
      listw = spdep::mat2listw(l.counters_spatial_weights),
      model = "within",
      spatial.error = "b",
      Hess = FALSE,
      lag = TRUE,
      effect = "individual"
    )

  print(summary(l.sdid_ind))




  p.sdid_ind <-
    splm::spml(
      value ~ improvement_physical + factor(date),
      data = p.bike_usage_panel_df,
      listw = spdep::mat2listw(p.counters_spatial_weights),
      model = "within",
      spatial.error = "b",
      Hess = FALSE,
      lag = TRUE,
      effect = "individual"
    )

  print(summary(p.sdid_ind))


  significance.results[i, 1] <- l.sdid_ind$coefficients[3]
  significance.results[i, 2] <- p.sdid_ind$coefficients[3]
}






##########################################################################################
# analysis:


p.plot.c <- p.counters_sf[6, ]
p.map <- mapview::mapview(p.counters_sf, col.regions = "red", xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE, label = TRUE) +
  mapview::mapview(p.plot.c, col.regions = "green", xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE, label = TRUE)
p.map

l.plot.c <- l.counters_sf[16, ]
l.map <- mapview::mapview(l.counters_sf, col.regions = "red", xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE, label = TRUE) +
  mapview::mapview(l.plot.c, col.regions = "green", xcol = "lon", ycol = "lat", crs = 4269, grid = FALSE, label = TRUE)
l.map
plot(p.counters_spatial_weights)



##########
# plot corraltion structures
par(mfrow = c(4, 3))
for (bandwidth in 1:4) {
  bandwidthcounter <- bandwidth
  listt <- c(1, 2.5, 5, 10)
  bandwidth <- listt[bandwidth]



  for (bandwidth.t in 1:3)
  {
    if (bandwidth.t == 3) {
      l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = bandwidth, kernel = "exponential", adaptive = FALSE)
    }
    if (bandwidth.t == 2) {
      l.counters_spatial_weights <- GWmodel::gw.weight(l.dd, bw = bandwidth, kernel = "gaussian", adaptive = FALSE)
    }
    if (bandwidth.t == 1) {
      listt <- c(5, 10, 15, 25)
      linearbandwidth <- listt[bandwidthcounter]
      l.counters_spatial_weights <- compute_linear_weight(l.dd, cutoff = linearbandwidth)
    }



    if (bandwidth.t == 3) {
      p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = bandwidth, kernel = "exponential", adaptive = FALSE)
    }
    if (bandwidth.t == 2) {
      p.counters_spatial_weights <- GWmodel::gw.weight(p.dd, bw = bandwidth, kernel = "gaussian", adaptive = FALSE)
    }
    if (bandwidth.t == 1) {
      listt <- c(5, 10, 15, 25)
      linearbandwidth <- listt[bandwidthcounter]
      p.counters_spatial_weights <- compute_linear_weight(p.dd, cutoff = linearbandwidth)
    }


    if (bandwidth.t == 1) {
      plot(p.dd, p.counters_spatial_weights, type = "p", xlab = "Distance (km)", ylab = "Spatial Weights", main = c(paste("Linear, cut-off=", linearbandwidth)), cex.lab = 1.5)
    }

    if (bandwidth.t == 2) {
      plot(p.dd, p.counters_spatial_weights, type = "p", xlab = "Distance (km)", ylab = "Spatial Weights", main = c(paste("Gaussian, bandwidth=", bandwidth)), cex.lab = 1.5)
    }

    if (bandwidth.t == 3) {
      plot(p.dd, p.counters_spatial_weights, type = "p", xlab = "Distance (km)", ylab = "Spatial Weights", main = c(paste("Exponential, bandwidth=", bandwidth)), cex.lab = 1.5)
    }
  }
}
##########################################################################################
## plot usage trends

## Visualisation


library(ggplot2)


l.improvement <- l.bike_usage %>%
  filter(improvement == 1) %>%
  arrange(channel_id, date) %>%
  group_by(channel_id) %>%
  slice_head(n = 1) %>%
  select(date, channel_id, value)

ggplot() +
  # untreated
  geom_line(
    data = l.bike_usage %>% filter(!channel_id %in% l.counters_with_physical_improvement),
    aes(x = date, y = value, group = channel_id), alpha = 0.2
  ) +

  # treated
  geom_line(
    data = l.bike_usage %>% filter(channel_id %in% l.counters_with_physical_improvement),
    aes(x = date, y = value, group = channel_id, color = channel_id)
  ) +
  geom_point(data = l.improvement, aes(x = date, y = value, color = channel_id, size = 1.5)) +
  # geom_vline(data = improvement, aes(xintercept = date), alpha = 0.2) +
  guides(
    alpha = "none",
    # color = "none",
    size = "none"
  ) +
  theme_bw() +
  scale_x_date(date_break = "1 month") +
  theme(axis.text.x = element_text(angle = 90))






p.improvement <- p.bike_usage %>%
  filter(improvement == 1) %>%
  arrange(channel_id, date) %>%
  group_by(channel_id) %>%
  slice_head(n = 1) %>%
  select(date, channel_id, value)

ggplot() +
  # untreated
  geom_line(
    data = p.bike_usage %>% filter(!channel_id %in% p.counters_with_physical_improvement),
    aes(x = date, y = value, group = channel_id), alpha = 0.2
  ) +

  # treated
  geom_line(
    data = p.bike_usage %>% filter(channel_id %in% p.counters_with_physical_improvement),
    aes(x = date, y = value, group = channel_id, color = channel_id)
  ) +
  geom_point(data = p.improvement, aes(x = date, y = value, color = channel_id, size = 1.5)) +
  # geom_vline(data = improvement, aes(xintercept = date), alpha = 0.2) +
  guides(
    alpha = "none",
    # color = "none",
    size = "none"
  ) +
  theme_bw() +
  scale_x_date(date_break = "1 month") +
  theme(axis.text.x = element_text(angle = 90))
##########################################################################################
# 3D plot of the results
## LYON
coords <- st_coordinates(l.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- l.counters_spatial_weights
z[1:85, ] <- 0
diag(z) <- l.counters_spatial_weights[15, ] * l.sdid_ind$coefficients[1]


# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-15, ]
plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


## PARIS
coords <- st_coordinates(p.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- p.counters_spatial_weights
z[1:14, ] <- 0
diag(z) <- p.counters_spatial_weights[6, ] * p.sdid_ind$coefficients[1]



# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-6, ]
plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)





##########################################################################################
# case studies with multiple improvements 3D plot of the results with linear compounding
## LYON
coords <- st_coordinates(l.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- l.counters_spatial_weights
z[1:85, ] <- 0
diag(z) <- (l.counters_spatial_weights[12, ] * l.sdid_ind$coefficients[1]) + (l.counters_spatial_weights[18, ] * l.sdid_ind$coefficients[1])


# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-18, ]
data <- data[-12, ]
plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


## PARIS
coords <- st_coordinates(p.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- p.counters_spatial_weights
z[1:14, ] <- 0
diag(z) <- (p.counters_spatial_weights[12, ] * p.sdid_ind$coefficients[1]) + (p.counters_spatial_weights[9, ] * p.sdid_ind$coefficients[1])



# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-12, ]
data <- data[-9, ]

plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.04, -0.03, -0.02, -0.01),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  colors = colorRamp(c("red", "yellow", "green"))
)




##########################################################################################
# case studies with multiple improvements 3D plot of the results with non-linear compounding
## LYON
coords <- st_coordinates(l.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- l.counters_spatial_weights
z[1:85, ] <- 0
z1 <- z2 <- z
diag(z1) <- (l.counters_spatial_weights[12, ] * l.sdid_ind$coefficients[1])
diag(z2) <- (l.counters_spatial_weights[18, ] * l.sdid_ind$coefficients[1])
combined_vector <- (-1 * diag(z1))^1 + (-1 * diag(z2))^0.9

diag(z) <- ((combined_vector) * -1)
# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-18, ]
data <- data[-12, ]
plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  intensity = data$z,
  color = c(-0.004, -0.003, -0.002, -0.001),
  colors = colorRamp(c("red", "yellow", "green"))
)


## PARIS
coords <- st_coordinates(p.counters_sf)
x <- coords[, "X"]
y <- coords[, "Y"]
z <- p.counters_spatial_weights
z[1:14, ] <- 0
z1 <- z2 <- z
diag(z1) <- (p.counters_spatial_weights[12, ] * p.sdid_ind$coefficients[1])
diag(z2) <- (p.counters_spatial_weights[9, ] * p.sdid_ind$coefficients[1])
combined_vector <- (-1 * diag(z1))^1 + (-1 * diag(z2))^0.9
diag(z) <- ((combined_vector) * -1)

# mesh
z <- diag(z)
data <- data.frame(x = x, y = y, z = z)
data <- data[-12, ]
data <- data[-9, ]

plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "mesh3d",
  intensity = data$z,
  color = c(-0.04, -0.03, -0.02, -0.01),
  colors = colorRamp(c("red", "yellow", "green"))
)


plot_ly() %>% add_trace(
  data = data, x = data$x, y = data$y, z = data$z, type = "contour",
  colors = colorRamp(c("red", "yellow", "green"))
)
