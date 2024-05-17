##############################################################################################################################################
## Example of GPS processing used in Riaz et al 2023_Global Ecology and Conservation
##############################################################################################################################################

## Read in relevant packages and set data path


library(ggplot2)
library(tidyverse)
library(aniMotum)
library(hms)
library(viridis)
library(scales)
library(readr)
library(rgdal)
library(grid)
library(raster)
library(stringr)
library(lubridate)
library(sp)
library(diveMove)
library(sf)
library(data.table)
library(parallel)
library(future)
library(future.apply)
library(adehabitatLT)
library(ggpubr)
library(terrainr)
library(ggOceanMaps)
library(ggspatial)

setwd("//saeri-file/users$/jriaz/Documents/Github/")

data_path <- ("//saeri-file/users$/jriaz/Documents/Github/")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Read in DF and some basic time and date formatting
## Coded this way to allow reading in raw data for multiple individuals at once based on pattern in naming

tbl_locs <- list.files(file.path(data_path), full.names = TRUE,  
                       pattern = "*-1-FastGPS.csv", recursive=T)  
tbl_locs

myfiles = do.call(rbind, lapply(tbl_locs, function(x) read.csv(x, stringsAsFactors = FALSE, skip = 3, header = TRUE)))

names(myfiles)

tbl_locs<-data.frame(depID = myfiles$Name, lat = myfiles$Latitude, lon=myfiles$Longitude, sats = myfiles$Satellites, haulout = myfiles$Hauled.Out, day = myfiles$Day, time = myfiles$Time, date = myfiles$InitTime)

gps<-subset(tbl_locs, tbl_locs$sats > 4)

BirdGPSData <- gps

BirdGPSData$Location <- "Bird Island"

BirdGPSData$Year <- "2019"

BirdGPSData$Sex <- "Female"

BirdGPSData$date<-as.POSIXct(paste(BirdGPSData$day, BirdGPSData$time),format = "%d-%b-%Y %H:%M:%S", tz="GMT")#format times

BirdGPSData$Type = "GPS" 

BirdGPSData$id = BirdGPSData$depID

BirdGPSData$lc = "G"


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Initial visualisation

Basemap <- basemap(limits = c(-70, -58, -56, -47), bathy.style = "rcb" , rotate = TRUE, grid.col = NA) + 
  theme(legend.position = "right")  +
  theme(legend.background = element_rect(fill = "white", colour = "white"),
        legend.key = element_rect(fill = "white")) 
Basemap

## Select 1 here to download basemap

Basemap + ggspatial::geom_spatial_point(data=BirdGPSData, aes(x=lon,y=lat,color=depID))


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#Pre-filtering of data: adapted from Ropert-Couldert et al 2020
# The next steps below are to clean the GPS data
# Step-by-step annotations provides

newd <- BirdGPSData %>%
  dplyr::group_by(id) %>%
  dplyr::arrange(date)


## Step 1 - Remove any duplicates or near duplicates that occurred within 10 s
##############################################################################################################################################

d1new <- newd %>%
  do(distinct(., date, .keep_all = TRUE)) %>%
  do(mutate(
    .,
    dup = difftime(date, lag(date), units = "secs") < 10 &
      Type == "GPS")) %>% 
  do(arrange(., order(date))) %>%
  filter(.,!dup) %>%
  dplyr::select(-dup)

d1.repnew <- nrow(newd) - nrow(d1new)
cat(sprintf("%d duplicate time &/or near duplicate location/time records removed\n", d1.repnew))


## Step 2 - Remove deployments with less than 20 observations
##############################################################################################################################################

min_obsnew <- 20 

d3new <- d1new %>%
  do(filter(., n() >= min_obsnew))
d3.rep.new <- n_groups(d1new) - n_groups(d3new)
cat(sprintf(
  paste(
    "%d individuals with fewer than",
    min_obsnew,
    "observations removed\n"
  ),
  d3.rep.new
))  



## Step 3 - Remove deployments that last less than 1 day
##############################################################################################################################################

min_daysnew <-1

d4new <- d3new %>%
  do(filter(
    .,
    difftime(max(date), min(date), units = "days") >= min_daysnew
  ))

d4.rep.new <- n_groups(d3new) - n_groups(d4new)
cat(sprintf(
  paste(
    "%d individuals with fewer than",
    min_daysnew,
    "deployment days removed\n"
  ),
  d4.rep.new
)) 



## Step 4 - Remove extreme travel rate locations (> 4 m/s)
##############################################################################################################################################

vmax <- 4

d5new <- d4new %>%
  do(mutate(., keep = grpSpeedFilter(cbind(date, lon, lat), speed.thr = vmax)))
d5.rep.new <- nrow(d4new) - sum(d5new$keep)
cat(sprintf(
  paste(
    "\n%d records with travel rates >",
    vmax,
    "m/s will be ignored by SSM filter\n"
  ),
  d5.rep.new
)) 

d5new <- d5new %>%
  dplyr::filter(keep == "TRUE")


##############################################################################################################################################
## Summary of the pre-filtering stage

indnew <- c(NA, NA, d3.rep.new, d4.rep.new, NA)
ignnew <- c(NA, NA, NA, NA, d5.rep.new)


repnew <- tibble(
  test = c(
    "duplicated_dates",
    paste("obs_lessthan_", min_obsnew, sep = ""),
    paste("deploy_lessthan_", min_daysnew, "_days", sep =
            ""),
    paste("obs_greaterthan_", vmax, "_ms", sep =
            "")
  ),
  obs_removed = c(d1.repnew, NA, NA, NA),
  individuals_removed = c(NA, d3.rep.new, d4.rep.new, NA),
  obs_flagged_to_ignore = c(NA, NA, NA, d5.rep.new)
)

repnew

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##Format according to aniMotum requirements

GPS_FilteredTracks <- d5new %>%
  dplyr :: select(id, date, lc, lon, lat) %>%
  ungroup() %>%
  dplyr::arrange(date, id) %>%
  dplyr::distinct()

unique(GPS_FilteredTracks$id)


ggplot()+ geom_point(data=GPS_FilteredTracks, aes(x=lon,y=lat,color=id)) +
  scale_colour_viridis_d() + theme(legend.position = "NULL")

########################################################## START OF SSM FITTING SCRIPT #######################################################
##############################################################################################################################################
##############################################################################################################################################
########################################################## START OF SSM FITTING SCRIPT #######################################################

# Run tracks through SSM. These are are GPS tracks so I've run with a correlated random walk with 15 minute time step
# Code written/structured to loop through individuals in batch 
# Grab the Predicted and the Fitted. Predicted tracks had to be rerouted around land barriers - with some individuals showing unrealistic movements over land

d3 <- data.frame(GPS_FilteredTracks) 

split_d <- split(d3, d3$id)[unique(d3$id)]


## Step 1 - Fit predicted and rerouted and do some checks
##############################################################################################################################################

plan(multisession)
ts <- Sys.time()
fit_all <- future_lapply(split_d, function(z) {
  fits <- fit_ssm(z, model="crw", 
                  control = ssm_control(optim=c("nlminb")), 
                  time.step = 0.25)
  Reroute <- route_path(fits, what = "predicted", map_scale = 10, buffer = 500) ## Reroute around land barriers
})
plan(sequential)
difftime(Sys.time(),ts)


fit_all <- bind_rows(fit_all)

plot(fit_all, what = "r") ## Check fit

redi <- osar(fit_all) ## Takes long time

plot(redi, "qq") / plot(redi, "acf")

unique(fit_all$id)

Predictedframe <- grab(fit_all, "predicted", as_sf=FALSE) 

Rerouteframe <- grab(fit_all, "rerouted", as_sf=FALSE) 

Predictedframe$Type <- "Predicted"

Rerouteframe$Type <- "Rerouted"


## Plot predicted vs rerouted to compare
Basemap + 
  ggspatial::geom_spatial_point(data = Predictedframe, aes(x = lon, y = lat, group = id), colour = "black") +
  ggspatial::geom_spatial_point(data = Rerouteframe, aes(x = lon, y = lat, group = id), colour = "red") 


## Step 2 - Fit fitted
##############################################################################################################################################

plan(multisession)
ts <- Sys.time()
fit_all_obs <- future_lapply(split_d, function(z) {
  ## time step of 24 = 1 per day
  fits <- fit_ssm(z, model="crw", 
                  control = ssm_control(optim=c("nlminb")), 
                  time.step = 0.25)
  return(grab(fits, "fitted", as_sf=FALSE))
})
plan(sequential)
difftime(Sys.time(),ts)

fit_all_obs <- bind_rows(fit_all_obs)

fit_all_obs$Type <- "Fitted"

FittedFrame <- fit_all_obs


Basemap + 
   ggspatial::geom_spatial_point(data = FittedFrame, aes(x = lon, y = lat, group = id), colour = "orange") 
  # ggspatial::geom_spatial_point(data = Predictedframe, aes(x = lon, y = lat, group = id), colour = "black") +
  # ggspatial::geom_spatial_point(data = Rerouteframe, aes(x = lon, y = lat, group = id), colour = "red") 



##############################################################################################################################################
##############################################################################################################################################
# Loop through each ID to plot Fitted vs Predicted/Rerouted to inspect and see if there are any dodgy tracks with the CRW.

Fit_Predicted_GPS <- FittedFrame %>%
  full_join(Rerouteframe)

for(i in unique(Fit_Predicted_GPS$id)){
  subID <- subset(Fit_Predicted_GPS, id == i)
  gi <- Basemap +
    ggspatial::geom_spatial_point(data = subID, aes(x = lon, y = lat, colour = Type), size = 1) + 
    scale_colour_manual(values = c("orange", "red")) + 
    xlab("Longitude") + ylab("Latitude") +
    # coord_sf(xlim = c(min(subID$lon), max(subID$lon)), ylim = c(min(subID$lat), max(subID$lat)), expand = TRUE) +
    theme_bw(10) + theme(strip.background = element_rect(fill="gray85"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black")) +
    theme(legend.position= "right") +
    ggtitle(paste0(i))
  ggsave(filename = sprintf('%s.png', i), plot = gi, width = 9, height = 9)
}


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

## Save the data

write_csv(Fit_Predicted_GPS, "GPS_SSM_output.csv")


##############################################################################################################################################
## THE END
##############################################################################################################################################


