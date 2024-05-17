##############################################################################################################################################
## Example of dive processing code used in Riaz et al 2023_Global Ecology and Conservation
##############################################################################################################################################

## Read in relevant packages and set data path

library(diveMove)
library(readr)
library(tidyverse)

rm(list=ls())

setwd("P:/Projects/DPLUS168_pinniped bycatch/07. Data/Bird Island_Fur Seal Dive Data/2019 data/")

datapath <- setwd("P:/Projects/DPLUS168_pinniped bycatch/07. Data/Bird Island_Fur Seal Dive Data/2019 data/")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Read in TDR data for 1 individual and do some formatting so it's good to go

TDR_ID <- "16A0600"

TDR1_data <- read_csv("Example_ID.csv", show_col_types = FALSE) 

TDR1_data <- TDR1_data %>%
  dplyr::select(Time, Depth)
  
TDR1_data <- TDR1_data[complete.cases(TDR1_data$Depth), ]

TDR1_data$Time <- as.POSIXct(strptime(TDR1_data$Time,"%H:%M:%S %d-%b-%Y", tz = "GMT"))

TDR1_data$Time <- as.POSIXct(TDR1_data$Time, format = "%d-%m-%Y %H:%M:%S", tz = "GMT")

ggplot() + geom_line(data = TDR1_data, aes(x= Time, y = Depth)) +
  scale_y_reverse()

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Start the processing using divemove

tdrX <- createTDR(time=TDR1_data$Time,
                  depth=TDR1_data$Depth,
                  dtime=1, file=datapath)

plotTDR(tdrX) ## Assess ZOC required for each ID

dcalib <- calibrateDepth(tdrX,
                         zoc.method="offset",
                         offset= 0.5
                         ,dive.thr = 3)

plotTDR(dcalib)

plotTDR(dcalib, diveNo=100:102) ##Subset to assess performance

Divemetrics <- diveStats(dcalib) ## Calculate dive statistics

Divemetrics$id <- TDR_ID
Divemetrics$Year <- "2019"
Divemetrics$Species <- "SAFS"

Divemetrics$Divenum <- 1:nrow(Divemetrics) ## Tally number of dives
Divemetrics$dive.tot <- max(Divemetrics$Divenum) ## Total number of dives in dataframe
Divemetrics$DiveEnd <- Divemetrics$begdesc + seconds(Divemetrics$divetim) # Dive end timestamp


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## Eyeball processed data

ggplot() + geom_line(data = Divemetrics, aes(x= begdesc, y = maxdep)) +
  scale_y_reverse() +
  geom_hline(yintercept = 0, colour  = "black") +
  theme_classic() +
  labs(x = "Date", y = "Depth (m)")

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

## Save the data

write.csv(Divemetrics, (paste("Dive_Processed_2019_",TDR_ID,".csv", sep = "")), row.names=FALSE)


##############################################################################################################################################
## THE END
##############################################################################################################################################
