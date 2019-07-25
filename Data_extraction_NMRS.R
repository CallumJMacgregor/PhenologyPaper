####################################################
####   Script for preparing data for analysis   ####
####################################################

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","rnrfa","rgdal","blighty","dplyr","raster","RColorBrewer","lme4","ggplot2","lubridate","mgcv","gridExtra","reshape2")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# we have moth data from two sources - NMRS and RIS
# we are going to use these data to estimate changes in a range of variables between two time windows:

# used in the original Mair et al. paper...

# distribution (percentage change in occupied hectads per year)
# abundance (percentage change in abundance index per year)

# dispersal ability? very tricky to replicate for moths... so maybe just use butterflies for this paper?
# habitat availability? suggest talk to Phil about using SSI for this... he has less faith in SSI for moths so again a case for just using butterflies?
#  - could either use national SSI (easier to extract)
#  - or try to replicate Mair more exactly by focussing it on colonised hectads!


# added in for this paper...

# phenology (number of days change in first-generation emergence date per year)
# northern range margin (km advance per year)


## selecting data to use
# we have different criteria for inclusion in the study for each data type 
# (NMRS = hectad-level annual presence/absence; RIS = site-level daily count data)

# we need to apply these to select out species to use from the hectad data, and species/site combinations to use from the site data

### NMRS data ####
# in a separate script we've calculated the metrics for the NMRS data
# now let's do it on the hectad data for moths, from the NMRS

## read in the data
NMRS_raw <- read.csv("../../Data/NMRS/NMRS Macro-moths 10km 1950-2014/Report - NMRS Macro-moths 10km 1950-2014.csv", header = TRUE)

summary(NMRS_raw)


# tweak the colnames to a standard

colnames(NMRS_raw) <- c("SCI_NAME","COMMON_NAME","HECTAD","YEAR")

NMRS_raw$fYEAR <- as.factor(NMRS_raw$YEAR)

# attach a column to indicate presence

NMRS_raw$PRESENCE <- 1

## next, we want to attach eastings and northings to this data so that it can be plotted
# as there are > 1000 repeats of each hectad let's generate a non-redundant list

NMRS_hectads <- ddply(NMRS_raw, .(HECTAD), summarise, COUNT = sum(PRESENCE))

# at the moment the dataset includes a whole bunch of records from sites in Northern Ireland, and some in the Channel Islands.
# that's nice but I don't want them - they will cause problems elsewhere (esp. when analysing range-margin stuff)

# the question is, how to pick them out and get rid?
# first let's take first two characters of every grid ref: these are the 100km square

NMRS_hectads$SQUARE <- as.factor(substr(NMRS_hectads$HECTAD, 1,2))

# now, all sites in GB have two letters, sites in NI have one letter and one number
# we can use this to separate them

NMRS_hectadsGB <- NMRS_hectads[which(!grepl("[[:digit:]]", NMRS_hectads$SQUARE)),]

# to get rid of Channel Islands, we need to get rid of squares WA and WV

NMRS_hectadsGB <- NMRS_hectadsGB[which(!(NMRS_hectadsGB$SQUARE %in% c("WA","WV"))), ]


# now generate a vector containing these site numbers only

# first we need to drop unused levels (i.e. those from NI)
NMRS_hectadsGB <- droplevels(NMRS_hectadsGB)


# now we are ready to parse the easting/northings

latlon <- osg_parse(NMRS_hectadsGB$HECTAD)

NMRS_hectadsGB$EASTING <- latlon[[1]]
NMRS_hectadsGB$NORTHING <- latlon[[2]]

### we now have eastings and northings! let's do a brief visual check that all's well by plotting them on a map...

# to do this we actually also need these with fewer digits - specifically, in units of 1 km
NMRS_hectadsGB$east.km <- NMRS_hectadsGB$EASTING/1000
NMRS_hectadsGB$north.km <- NMRS_hectadsGB$NORTHING/1000

# plot base map to test
blighty()

# plot all records to test
# we need to offset all points slightly for them to be centred in the middle of hectads, rather than the bottom-left corner
points(5+NMRS_hectadsGB$east.km,5+NMRS_hectadsGB$north.km,pch=16) # this looks a bit messy but looks fine in a zoom window


# great - but this is a total mess!
# we also need to merge these eastings and northings back into the main data
NMRS_hectadsGB <- NMRS_hectadsGB[,-2]


# we can now use this to trim down the other data to only include these sites
NMRS_GB <- merge(NMRS_raw,NMRS_hectadsGB)



# let's plot an example species to get a distribution

NMRS_GB.bb <- NMRS_GB[which(NMRS_GB$COMMON_NAME=="Burnished Brass"), ]
blighty()
points((5+NMRS_GB.bb$east.km),(5+NMRS_GB.bb$north.km),pch=16)
title(main = "Burnished Brass")


# let's also try looping this across the range of years of interest (1986-1995 and 2001-2010) to see whether any expansion is visible

# set up the period of interest
years <- 1995:2014

for (x in years){
  year <- NMRS_GB.bb[which(NMRS_GB.bb$YEAR==x), ]
  png(paste0("Data/NMRS/Burnished Brass/Range/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points((5+year$east.km),(5+year$north.km),pch=16)
  title(main = x)
  dev.off()
}


## now restrict the data down to fit data-quality filters


### I want to select only 'heavily-recorded' hectads to take forwards
# this follows Hickling et al 2005, 'common squares' (see Hassall & Thompson 2010)
# the idea of this is that we only use hectads where we are more confident that we're looking at a true absence rather than a lack of recording
# i.e. we check every hectad individually and only use those where a high enough proportion of the regional fauna
# has been recorded at least once in each time window




## first, let's create a non-redundant list of recording in each hectad in each year, dropping the species information

NMRS_rec.hectad.years <- ddply(NMRS_GB, .(HECTAD, YEAR), numcolwise(mean))

# for inspection, let's plot the total hectad-level recording across all years

years.extended <- 1980:2014

for (n in years.extended){
  year <- NMRS_rec.hectad.years[which(NMRS_rec.hectad.years$YEAR==n), ]
  png(paste0("Data/NMRS/HectadRecording/",n,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points(year$east.km,year$north.km,pch=16)
  title(main = n)
  dev.off()
}

# as for butterflies, we can see that recording was pretty sparse at some points
# and in some places, especially in the North and in Lincolnshire (tut tut)
# but gets fairly robust especially in most recent years

# let's actually count the number of recorded hectads and the number of species*hectad records in each year 

recording <- data.frame(YEAR = numeric(),
                        HECTADS = numeric(),
                        RECORDS = numeric())

for (n in years.extended){
  year.hec <- NMRS_rec.hectad.years[which(NMRS_rec.hectad.years$YEAR==n), ]
  HECTADS <- nrow(year.hec)
  
  year.rec <- NMRS_GB[which(NMRS_GB$YEAR==n), ]
  RECORDS <- nrow(year.rec)
  
  YEAR <- n
  
  out <- data.frame(cbind(YEAR,HECTADS,RECORDS))
  recording <- rbind(recording,out)
  
  
}

summary(recording)


# now plot the number of hectads and records over time
plot(HECTADS ~ YEAR, data = recording)
plot(RECORDS ~ YEAR, data = recording)


# it's not so clear-cut as for butterflies - the amount of recording is still increasing year on year
total.hecs <- max(recording$HECTADS)

recording$PERC <- recording$HECTADS*100/total.hecs

plot(PERC ~ YEAR, data = recording)

# after around 1995, we can say that over 60% of the maximum number of hectads have been recorded in all years for both datasets
# so we can justify using 1995 as the start point of a 20-year window

# therefore, our final target window is:

interval <- 1995:2014


# now, let's do some data cleaning
# we want to make sure that each hectad used was recorded at least once in both the first and second half of each interval
# we could even ask to only use hectads recorded in every single year, but things might get a bit sparse with this approach

# let's divide the target interval into two halves
inter1 <- 1995:2004
inter2 <- 2005:2014


# first follow Hickling et al directly by calculating which hectads were recorded at least once in each interval

# first select out the intervals
NMRS_rec.hectad.years.inter1 <- NMRS_rec.hectad.years[which(NMRS_rec.hectad.years$YEAR %in% inter1), ]
NMRS_rec.hectad.years.inter2 <- NMRS_rec.hectad.years[which(NMRS_rec.hectad.years$YEAR %in% inter2), ]

# then get rid of the year information to leave a non-redundant list of hectads recorded in each interval

NMRS_rec.hectads.inter1 <- ddply(NMRS_rec.hectad.years.inter1, .(HECTAD), numcolwise(mean))
NMRS_rec.hectads.inter2 <- ddply(NMRS_rec.hectad.years.inter2, .(HECTAD), numcolwise(mean))


# plot these on a map

blighty()
points(NMRS_rec.hectads.inter2$east.km, NMRS_rec.hectads.inter2$north.km, pch=16, col="red")
points(NMRS_rec.hectads.inter1$east.km, NMRS_rec.hectads.inter1$north.km, pch=16)


# we can see that almost every hectad was recorded in interval 2, and most in interval 1
# though parts of the Highlands, and again Lincolnshire, are slightly sparse

## now we want to recombine these datasets in a way that ditches anything not present in both

NMRS_rec.hectads.inter1.vec <- as.character(NMRS_rec.hectads.inter1$HECTAD)
NMRS_rec.hectads.inter2.vec <- as.character(NMRS_rec.hectads.inter2$HECTAD)

NMRS_rec.hectads.full <- append(NMRS_rec.hectads.inter1.vec,NMRS_rec.hectads.inter2.vec)

NMRS_rec.hectads.full <- data.frame(cbind(NMRS_rec.hectads.full,1))
colnames(NMRS_rec.hectads.full) <- c("HECTAD","INTERVALS")
NMRS_rec.hectads.full$HECTAD <- as.factor(as.character(NMRS_rec.hectads.full$HECTAD))
NMRS_rec.hectads.full$INTERVALS <- as.numeric(as.character(NMRS_rec.hectads.full$INTERVALS))


NMRS_rec.hectads.good <- ddply(NMRS_rec.hectads.full, .(HECTAD), numcolwise(sum))

NMRS_rec.hectads.good <- NMRS_rec.hectads.good[which(NMRS_rec.hectads.good$INTERVALS > 1), ]

# this leaves us with a final list of 2353 'good hectads'

## now we can use this list of 'good' hectads to select out data from the main set to use
# first create a vector of 'good' hectads

# drop unused levels
NMRS_rec.hectads.good <- droplevels(NMRS_rec.hectads.good)

# generate the vector
NMRS_good.hectads <- levels(NMRS_rec.hectads.good$HECTAD)

# pick out the good records
NMRS_GB_good <- NMRS_GB[which(NMRS_GB$HECTAD %in% NMRS_good.hectads), ]

# it's promising that we lose only just over 1% of records

## we can also further restrict to the hectads that have records of a minimum percentage of the regional species richness
# this is going to take quite a complex loop to construct:
# for every hectad in the NMRS_GB_good dataframe, we need to calculate the observed species richness
# then we need to calculate pairwise distances to all other recorded hectads, pick out the nearest 100, 
# and calculate the regional species richness of these 100
# and finally turn this into a percentage of regional SR that is recorded in the hectad

# start with the list of recorded hectads, re-attach the eastings and northings
NMRS_good <- ddply(NMRS_rec.hectad.years, .(HECTAD,EASTING,NORTHING), summarise,
                  COUNT = sum(PRESENCE))

# select out only the good hectads
NMRS_good <- NMRS_good[which(NMRS_good$HECTAD %in% NMRS_good.hectads), ]

# now start to construct the loop:
hec.records <- data.frame(YEAR = numeric(),
                          HECTAD = factor(),
                          EASTING = numeric(),
                          NORTHING = numeric(),
                          hec.SR = numeric(),
                          reg.SR = numeric(),
                          hec.perc = numeric())
i <- 1

for (x in NMRS_good.hectads){
  if (round(i, -1) == i){
    print(i)
  }
  hec <- NMRS_good[which(NMRS_good$HECTAD == x), ]    # pick out details of focal hectad
  candidates <- NMRS_good[which(NMRS_good$HECTAD != x), ]   # pick out details of all others
  
  hec_east <- hec[[1,2]]     # easting of target
  hec_north <- hec[[1,3]]    # northing of target
  
  candidates$east_diff <- candidates$EASTING - hec_east             # longitudinal difference
  candidates$north_diff <- candidates$NORTHING - hec_north          # latitudinal difference
  
  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
  
  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
  
  closest <- candidates[1:100,] # select out closest 100
  
  # now we want to calculate species richness from the hectad and from the region in each year
  for (n in years){
    year.recs <- NMRS_GB_good[which(NMRS_GB_good$YEAR == n), ]
    
    hec.recs <- year.recs[which(year.recs$HECTAD == x), ] # pull out hectad records
    hec.SR <- nlevels(droplevels(hec.recs$COMMON_NAME))       # calculate species richness
    
    reg.recs <- year.recs[which(year.recs$HECTAD %in% droplevels(closest$HECTAD)), ]  # pull out region records
    reg.recs <- rbind(reg.recs,hec.recs)    # add in hectad records (they're part of the regional richness too!)
    reg.SR <- nlevels(droplevels(reg.recs$COMMON_NAME))
    
    hec.perc <- hec.SR*100/reg.SR
    
    out <- cbind(x,n,hec_east,hec_north,hec.SR,reg.SR,hec.perc)
    hec.records <- rbind(hec.records,out)
  }    
  i <- i+1
}

colnames(hec.records) <- c("HECTAD","YEAR","EASTING","NORTHING","HECTAD.SR","REGION.SR","PERCENT.RECORDED")
summary(hec.records)

# make things that should be numeric, numeric
hec.records$YEAR <- as.numeric(as.character(hec.records$YEAR))
hec.records$EASTING <- as.numeric(as.character(hec.records$EASTING))
hec.records$NORTHING <- as.numeric(as.character(hec.records$NORTHING))
hec.records$HECTAD.SR <- as.numeric(as.character(hec.records$HECTAD.SR))
hec.records$REGION.SR <- as.numeric(as.character(hec.records$REGION.SR))
hec.records$PERCENT.RECORDED <- as.numeric(as.character(hec.records$PERCENT.RECORDED))

summary(hec.records)

# most of the hectad*year combinations are "well recorded" but only slightly under a quarter are "heavily recorded"

# this data takes a while to generate so let's back it up
write.table(hec.records, "Data/NMRS/HectadRecording/RecordingLevels.txt", row.names = F)

# now we want to generate a single row for each hectad with the minimum recording level in any year
# we'll use this to assign hectads to different recording levels
hec.records$HECTAD <- as.character(hec.records$HECTAD)


hec.records.overall <- ddply(hec.records, .(HECTAD,EASTING,NORTHING), summarise,
                             MIN.PERC = min(PERCENT.RECORDED),
                             MEAN.PERC = mean(PERCENT.RECORDED),
                             MED.PERC = median(PERCENT.RECORDED),
                             LQ.PERC = quantile(PERCENT.RECORDED)[2])


summary(hec.records.overall)


# now we have a LOT of hectads which aren't even recorded in every year showing up,
# and a large proportion of hectads fail to reach either the well-recorded or heavily-recorded thresholds when using the 'in every year' (i.e. min) criterion

# using the median criterion (i.e. 'in a median year'), we still have most cells at least well-recorded





# so let's try plotting the geographic spread of recording density a couple of ways
# first, label each hectad according to its maximum recording level
hec.records.overall$RECORDING <- as.factor(ifelse(hec.records.overall$LQ.PERC >= 25, "Heavily recorded",
                                                  ifelse(hec.records.overall$LQ.PERC >= 10, "Well recorded",
                                                         ifelse(hec.records.overall$LQ.PERC >0, "Recorded",
                                                                "Not recorded"))))

summary(hec.records.overall$RECORDING)

blighty()
points(5+(hec.records.overall$EASTING/1000),5+(hec.records.overall$NORTHING/1000),pch=16,col=c("Black","Grey90","Grey70","Red")[hec.records.overall$RECORDING])
legend("topleft", col=c("Grey90","Grey70","Red","Black"), pch=16,
       legend=c("Not recorded","Recorded","Well recorded","Heavily recorded"))

# again, many less well-recorded cells are in the uplands, but there are pockets elsewhere

# next, do a gradient for percentage

cols <- brewer.pal(11, "RdYlBu")
pal <- colorRampPalette(rev(cols))

hec.records.overall$ORDER <- findInterval(hec.records.overall$LQ.PERC, sort(hec.records.overall$LQ.PERC))

blighty()

points(5+(hec.records.overall$EASTING/1000),5+(hec.records.overall$NORTHING/1000),pch=16,
       col = pal(nrow(hec.records.overall))[hec.records.overall$ORDER])
legend("topleft", col = pal(2), pch=16,
       legend=c("Less recording","More recording"))


# save these plots
png(paste0("Data/NMRS/HectadRecording/RecordingLevels.png"), width = 1600, height = 800, units = "px", bg = "white")
par( mfrow = c(1,2), oma = c(0,0,2,0))

blighty()
points(5+(hec.records.overall$EASTING/1000),5+(hec.records.overall$NORTHING/1000),pch=16,col=c("Black","Grey90","Grey70","Red")[hec.records.overall$RECORDING])
legend("topright", col=c("Grey90","Grey70","Red","Black"),pch=16,cex=2,
       legend=c("Not recorded","Recorded","Well recorded","Heavily recorded"))

blighty()
points(5+(hec.records.overall$EASTING/1000),5+(hec.records.overall$NORTHING/1000),pch=16,
       col = pal(nrow(hec.records.overall))[hec.records.overall$ORDER])
legend("topleft", col = pal(2), pch=16,cex=2,
       legend=c("Less recording","More recording"))

dev.off()       


# so - there is a clear bias towards the South
# and also towards some northern cities (Manchester, Edinburgh)


# we can also plot out the hectads so we can visualise where their region lies...
# but this takes a LONG time so once I've done it I'll blank it out!
#for (x in NMRS_good.hectads){
#  hec <- NMRS_good[which(NMRS_good$HECTAD == x), ]    # pick out details of focal hectad
#  candidates <- NMRS_good[which(NMRS_good$HECTAD != x), ]   # pick out details of all others
#  
#  hec_east <- hec[[1,2]]     # easting of target
#  hec_north <- hec[[1,3]]    # northing of target
#  
#  candidates$east_diff <- candidates$EASTING - hec_east             # longitudinal difference
#  candidates$north_diff <- candidates$NORTHING - hec_north          # latitudinal difference
#  
#  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
#  
#  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
#  
#  closest <- candidates[1:100,] # select out closest 100
#  
#  # plot out the hectads
#  png(paste0("Data/NMRS/HectadRecording/HectadRegions/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
#  
#  blighty()
#  points(5+(candidates$EASTING/1000), 5+(candidates$NORTHING/1000),pch=16)
#  points(5+(closest$EASTING/1000), 5+(closest$NORTHING/1000),pch=16,col="red")
#  points(5+(hec_east/1000),5+(hec_north/1000),pch=16,col="royalblue")
#  title(main = x)
#  
#  dev.off()
#}

# so the main dataframe (NMRS_GB_good) is effectively our dataframe for 'recorded'
# let's generate an additional one for each of 'well recorded' and 'heavily recorded'
hec.records.overall$HECTAD <- as.factor(hec.records.overall$HECTAD)

rec.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC > 0), ]
well.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC >= 10), ]
heavy.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC >= 25), ]

NMRS_GB_rec <- NMRS_GB_good[which(NMRS_GB_good$HECTAD %in% droplevels(rec.hecs$HECTAD)), ]
NMRS_GB_well <- NMRS_GB_good[which(NMRS_GB_good$HECTAD %in% droplevels(well.hecs$HECTAD)), ]
NMRS_GB_heavy <- NMRS_GB_good[which(NMRS_GB_good$HECTAD %in% droplevels(heavy.hecs$HECTAD)), ]


### selecting species ####
## so, we have our good hectads: now we need to select which species to use from those that remain in the data

## we'll need to repeat this for all three recording levels. Let's first do recorded.

### recorded hectads ####
# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
NMRS_GB_good.inters <- NMRS_GB_rec[which(NMRS_GB_rec$YEAR %in% years), ]

# get rid of the year information, keeping only which hectads each species was recorded in:
NMRS_GB_good.hecs <- ddply(NMRS_GB_good.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                          PRESENCE = 1)

# and now count the number of hectads per species
species.hecs <- ddply(NMRS_GB_good.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                      HECTADS = sum(PRESENCE))

print(species.hecs)

# we need a list of species with fewer than 20 hectads of records:
under.recorded.species <- species.hecs[which(species.hecs$HECTADS < 20), ]

print(under.recorded.species)

under.recorded <- levels(droplevels(under.recorded.species$COMMON_NAME))

# this chops out a bunch of things - most either very rare (e.g. Barberry Carpet, Fisher's Estuarine)
# or migrant (e.g. Spurge Hawk-moth)

# we also want to manually create a list of migrants (urgh!)
# I'm going to write out all species into a csv, manually add a column for migrant or not, and read it back in
# this will additionally let me get to the bottom of another rule - no species aggs!

write.table(species.hecs,"Data/NMRS/allspecies.csv", row.names = F, sep = ",")

migrants.frame <- read.csv("Data/NMRS/migrants.csv", header = T)
summary(migrants.frame)

migrants.only <- migrants.frame[which(migrants.frame$MIGRANT_ADVENTIVE==T), ]
migrants <- levels(droplevels(migrants.only$COMMON_NAME))

aggregates.only <- migrants.frame[which(migrants.frame$AGGREGATE==T), ]
aggregates <- levels(droplevels(aggregates.only$COMMON_NAME))

spec.to.cut <- NULL

for (x in c(under.recorded,migrants,aggregates)){
  if (!(x %in% spec.to.cut)){
    spec.to.cut <- append(spec.to.cut,x)
  }
}


# next we want to focus on southerly species: i.e. those that reach a northern range margin in the UK, are likely climate-limited here,
# and have the capacity to expand northwards

# we won't completely drop the non-southerly species, but we will create sub-datasets containing *only* southerly species at the end,
# meaning that we need to identify now which species are southerly and which aren't

# the rule for this is that the species' northern range margin in each period can't be within 100 km of John O'Groats
# so obviously for this, we need to estimate the northern range margins!

# these are the mean northing of the 10 most northerly occupied hectads in the interval

# and now calculate the margins for each species
margins <- data.frame(SCI_NAME = factor(),
                      NORTHING = numeric(),
                      COMMON_NAME = factor())

species.list <- levels(droplevels(species.hecs$COMMON_NAME))
species.list <- species.list[which(!(species.list %in% spec.to.cut))]

for (x in species.list){
  sp.dat <- NMRS_GB_good.hecs[which(NMRS_GB_good.hecs$COMMON_NAME==x), ]
  northern10 <- top_n(sp.dat, 10, sp.dat$NORTHING)
  species.rm <- ddply(northern10, .(COMMON_NAME), numcolwise(mean))
  species.rm <- species.rm[,c(1,3)]
  colnames(species.rm) <- c("COMMON_NAME","NORTHING")
  
  species.rm$SCI_NAME <- sp.dat[[1,2]]
  
  
  margins <- rbind(margins,species.rm)
}

summary(margins)


## now we can calculate which species to drop:
# the hectad containing John O'Groats is ND37, so the northing is 970000
# 100 km south of here would be 870000 - around Lossiemouth, just north of Inverness
# we want to check whether each species' NRM is south of here in BOTH intervals

margins$EXCLUDE <- as.factor(ifelse(margins$NORTHING < 870000, "Southerly","Ubiquitous"))
summary(margins$EXCLUDE)

head(margins)

# here we lose quite a range of species with distributions that extend into Scotland
# some are truly ubiquitous -
# e.g. Brimstone Moth
# some have patchy distributions but happen to be present in northern Scotland -
# e.g. Common Heath
# and some are northerly-distributed
# e.g. Argent & Sable

# none of these are climate-limited (southerly-distributed) so we're fine to exclude them in a sub-dataset
# however we *don't* want to add them to the list of species to cut, because we will also want the full data for them

margins.drop <- margins[which(margins$EXCLUDE == "Ubiquitous"), ]

margins.spec.drop <- levels(droplevels(margins.drop$COMMON_NAME))


## at this point we also want to check the elevational distributions of the remaining species
# montane species may go uphill rather than north, which is much harder to track...
# so we want to exclude anything for which the mean elevation of occupied hectads is >200m

# we have a raster containing DEM data at hectad level for the entire British Isles (except a few cells in Shetland) 
# start by reading this in

DEM <- raster("Data/dem.mean.asc")
plot(DEM)

# first I want to get a basic distribution for each species (i.e. every cell in which it's ever been recorded)
# then for each species, I want to extract the elevation of all of these cells
# and take a mean

# it's probably easiest to do this within a loop

elevations <- data.frame(SCI_NAME = factor(),
                         COMMON_NAME = factor(),
                         HECTADS = numeric(),
                         ELEVATION = numeric(),
                         ELEVATION.WEIGHTED = numeric())

for (x in species.list){
  spec <- NMRS_GB[which(NMRS_GB$COMMON_NAME == x),]
  SCI_NAME <- as.character(spec[[1,2]])
  COMMON_NAME <- x
  
  spec.sum <- ddply(spec, .(SCI_NAME,COMMON_NAME,HECTAD,EASTING,NORTHING), summarise,
                    YEARS = sum(PRESENCE), PRESENCE = mean(PRESENCE))
  
  # at present the grid ref points to the corner of the cell 
  # which is an issue when extracting using raster::extract, so we add a tiny amount to point to somewhere just inside each cell
  spec.sum$EASTING5 <- spec.sum$EASTING + 5             
  spec.sum$NORTHING5 <- spec.sum$NORTHING + 5            
  
  spec.sum$ELEVATION <- extract(DEM, spec.sum[,8:9])
  
  # occasionally this returns a small number of NAs.
  # These are always coastal squares which are mostly water, so it's not unreasonable to just set them to 0
  
  spec.sum$ELEVATION[is.na(spec.sum$ELEVATION)] <- 0
  
  
  HECTADS <- nrow(spec.sum)
  ELEVATION <- mean(spec.sum$ELEVATION)
  ELEVATION.WEIGHTED <- weighted.mean(spec.sum$ELEVATION,spec.sum$YEARS)
  
  # export the final frame
  
  out <- data.frame(cbind(SCI_NAME,COMMON_NAME,HECTADS,ELEVATION,ELEVATION.WEIGHTED))
  elevations <- rbind(elevations,out)
}


elevations$HECTADS <- as.numeric(as.character(elevations$HECTADS))
elevations$ELEVATION <- as.numeric(as.character(elevations$ELEVATION))
elevations$ELEVATION.WEIGHTED <- as.numeric(as.character(elevations$ELEVATION.WEIGHTED))

summary(elevations)

# we have our list of mean elevations at hectad level
# now we want to drop out anything where the mean elevation is >200m

elevs.drop <- as.character(elevations$COMMON_NAME[which(elevations$ELEVATION > 200)])

# and add these to the list of species to cut
for (x in elevs.drop){
  if(!(x %in% spec.to.cut)){
    spec.to.cut <- append(spec.to.cut,x)
  }
}

# finally, in an early iteration of this analysis, I found that high up among the expanding species are almost all of the Clearwings (Sesiidae)
# these species are rarely recorded except by means of pheromone lures - most of which have been developed/improved fairly recently
# therefore I think these increases are due to increased targeted recording that is independent of general trends in recording effort
# I'm going to exclude all these species at this point:

pheromones.spec <- elevations[which(grepl("Clearwing",elevations$COMMON_NAME)), ]

pheromones <- levels(droplevels(pheromones.spec$COMMON_NAME))

# there are other species with pheromone lures but they are either already excluded from the data at this point (e.g. Emperor moth - ubiquitous)
# or the lures are not commercially available (e.g. Barred Tooth-striped) and therefore will not have caused a major upturn in recording

for (x in pheromones){
  if(!(x %in% spec.to.cut)){
    spec.to.cut <- append(spec.to.cut,x)
  }
}

# finally, we only want species that have been recorded in a minimum proportion of the 20 years
# count how many years each species has a record (in any hectad) for 

NMRS_GB_good.years <- ddply(NMRS_GB_good.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

NMRS_GB_good.years.count <- ddply(NMRS_GB_good.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(NMRS_GB_good.years.count)

# now retain only species recorded in all 20 years
insufficient.years.spec <- NMRS_GB_good.years.count[which(NMRS_GB_good.years.count$YEARS < 20), ]

insufficient.years <- levels(droplevels(insufficient.years.spec$COMMON_NAME))

for (x in insufficient.years){
  if(!(x %in% spec.to.cut)){
    spec.to.cut <- append(spec.to.cut,x)
  }
}

# we finally have a final list of species to cut - so we can derive a list of species to use!

spec.to.use <- NULL

for (x in species.list){
  if(!(x %in% spec.to.cut)){
    spec.to.use <- append(spec.to.use,x)
  }
}

### variable extraction ####

# for each species left, we want to extract the distribution change and northern range margin advance from this data
# and the abundance change and phenology advance from the RIS data

# since we have the data read in already, let's start with the distribution change
# we need to calculate the number of occupied hectads (from those chosen for retention based on recording level) in each time period
# and the % change per year

## select out data from years of interest
hectads.years <- NMRS_GB_good[which(NMRS_GB_good$YEAR %in% years), ]

# calculate how many hectads were recorded in each year
# first count how many species recorded in each hectad per year
hectads.rec.years <- ddply(hectads.years, .(HECTAD,YEAR,EASTING,NORTHING), summarise, SPECIES = sum(PRESENCE), PRESENCE = 1)

# then how many hectads in each year
hectads.rec <- ddply(hectads.rec.years, .(YEAR), summarise, HECTADS = sum(PRESENCE))
colnames(hectads.rec)[2] <- "ANNUAL.TOTAL"

# and now pick out species of interest
hectads.years <- hectads.years[which(hectads.years$COMMON_NAME %in% spec.to.use),]

# and count the number of hectads in each year for each species
hectads.species <- ddply(hectads.years, .(SCI_NAME,COMMON_NAME,YEAR), summarise, HECTADS = sum(PRESENCE))

# now merge in the annual total number of hectads
hectads.annualtotals <- merge(hectads.species,hectads.rec)

# and calculate the percentage of all recorded hectads in which each species was recorded each year
hectads.annualtotals$PERC.REC <- hectads.annualtotals$HECTADS*100/hectads.annualtotals$ANNUAL.TOTAL


## we also want a scaled percentage (i.e. PERC compared to PERC in 1995)
# this gives us an estimate of distribution size change relative to the size of the starting distribution

# extract out the 1995 estimates
hectads.1995 <- hectads.annualtotals[which(hectads.annualtotals$YEAR == 1995), c(2:3,6)]
colnames(hectads.1995)[3] <- "START.PERC.REC"

# merge them back in 
hectads.annualtotals <- merge(hectads.annualtotals,hectads.1995)

# and calculate the scaled distribution change
hectads.annualtotals$PERC.CHANGE <- hectads.annualtotals$PERC.REC*100/hectads.annualtotals$START.PERC.REC

## now we need change over time for each species - this is a linear regression of PERC against YEAR
# let's set it up inside a loop to export slope and P-value, as well as a plot

distrib.slopes <- data.frame(SCI_NAME = factor(),
                             COMMON_NAME = factor(),
                             DISTRIB.SLOPE = numeric(),
                             DISTRIB.SE = numeric(),
                             DISTRIB.CHI = numeric(),
                             DISTRIB.P = numeric())

for (x in spec.to.use){
  spec <- hectads.annualtotals[which(hectads.annualtotals$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(PERC.CHANGE ~ YEAR, data = spec)
  
  DISTRIB.SLOPE <- round(summary(test1)$coefficients[2,1], 2)
  DISTRIB.SE <- round(summary(test1)$coefficients[2,2], 2)
  DISTRIB.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  DISTRIB.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/Distribution/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PERC.CHANGE ~ spec$YEAR, 
       ylim = c(0,max(300,max(spec$PERC.CHANGE))),
       xlab = "YEAR", ylab = "Annual change in % of recorded hectads occupied (scaled to 1995 level)")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (DISTRIB.P == 0){
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p = ",DISTRIB.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, DISTRIB.SLOPE, DISTRIB.SE, DISTRIB.CHI, DISTRIB.P))
  distrib.slopes <- rbind(distrib.slopes, out)
  
}

distrib.slopes$DISTRIB.SLOPE <- as.numeric(as.character(distrib.slopes$DISTRIB.SLOPE))
distrib.slopes$DISTRIB.SE <- as.numeric(as.character(distrib.slopes$DISTRIB.SE))
distrib.slopes$DISTRIB.CHI <- as.numeric(as.character(distrib.slopes$DISTRIB.CHI))
distrib.slopes$DISTRIB.P <- as.numeric(as.character(distrib.slopes$DISTRIB.P))

summary(distrib.slopes)

# this puts V-moth, The Four-Spotted and Striped Lychnis among the biggest losers, which makes sense
# and known expanding species like Cypress Carpet and Channel Islands Pug among the biggest winners


## and now northern range margins
# we want to select out the ten northernmost hectads for each species in each interval
# retaining ALL hectads that are tied for tenth (if there are scattered cells above a hard margin, this pulls towards the hard margin)
# and take the mean latitude of these cells


# so, for each year, first select out the northernmost hectads for each species and take their mean
# seed an output table

margins <- data.frame(SCI_NAME = factor(),
                      COMMON_NAME = factor(),
                      YEAR = numeric(),
                      NORTHING = numeric())


for (x in spec.to.use){
  spec <- hectads.years[which(hectads.years$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,2]])
  
  for (n in years){
    YEAR <- n
    spec.year <- spec[which(spec$YEAR==n), ]
    northern10 <- top_n(spec.year, 10, spec.year$NORTHING)
    spec.year.rm <- ddply(northern10, .(SCI_NAME,COMMON_NAME,YEAR), summarise, NORTHING = mean(NORTHING))
    
    NORTHING <- as.character(spec.year.rm[[1,4]])
    
    out <- data.frame(cbind(SCI_NAME,COMMON_NAME,YEAR,NORTHING))
    margins <- rbind(margins,out)
  }
}

margins$YEAR <- as.numeric(as.character(margins$YEAR))
margins$NORTHING <- as.numeric(as.character(margins$NORTHING))

summary(margins)

cols.pal <- colorRampPalette(colors=c("#000000", "#f0f0f0", "#FF0000"))
cols.vec <- cols.pal(20)

cols <- data.frame(cbind(years,cols.vec))
colnames(cols) <- c("YEAR","COL")

# plot this on those maps as lines over the ranges
for (x in spec.to.use){
  spec <- hectads.years[which(hectads.years$COMMON_NAME==x), ]
  spec.rm <- margins[which(margins$COMMON_NAME==x), ]
  spec.rm <- merge(spec.rm,cols)
  png(paste0("Data/Derived/RangeMargins/Margins/Moths/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points((5+spec$EASTING/1000),(5+spec$NORTHING/1000),pch=16, col="black")
  title(main = x)
  for (n in 1995:2014){
    abline(h=(spec.rm[which(spec.rm$YEAR==n),4]/1000), col=as.character(cols[which(cols$YEAR==n),2]), lwd = 2)
  }
  legend("bottomright", col=cols.pal(3), lwd=2,
         legend=c("1995","2005","2014"))
  
  dev.off()
} 


# let's try calculating change in range margins, then:
# here we want the absolute change in km per year

margin.slopes <- data.frame(SCI_NAME = factor(),
                            COMMON_NAME = factor(),
                            MARGIN.SLOPE = numeric(),
                            MARGIN.SE = numeric(),
                            MARGIN.CHI = numeric(),
                            MARGIN.P = numeric())

for (x in spec.to.use){
  spec <- margins[which(margins$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(NORTHING ~ YEAR, data = spec)
  
  MARGIN.SLOPE <- round((summary(test1)$coefficients[2,1])/1000, 2)
  MARGIN.SE <- round((summary(test1)$coefficients[2,2])/1000, 2)
  MARGIN.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  MARGIN.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/RangeMargins/Change/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$NORTHING ~ spec$YEAR, 
       ylim = c(70000,1200000),
       xlab = "YEAR", ylab = "Annual Northing of Northern Range Margin")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (MARGIN.P == 0){
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p = ",MARGIN.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, MARGIN.SLOPE, MARGIN.SE, MARGIN.CHI, MARGIN.P))
  margin.slopes <- rbind(margin.slopes, out)
  
}

margin.slopes$MARGIN.SLOPE <- as.numeric(as.character(margin.slopes$MARGIN.SLOPE))
margin.slopes$MARGIN.SE <- as.numeric(as.character(margin.slopes$MARGIN.SE))
margin.slopes$MARGIN.CHI <- as.numeric(as.character(margin.slopes$MARGIN.CHI))
margin.slopes$MARGIN.P <- as.numeric(as.character(margin.slopes$MARGIN.P))

summary(margin.slopes)



# we have a couple of retractions e.g. Pale Mottled Willow where northerly populations have gone extinct
# and some big expansions e.g. Narrow-bordered Bee Hawkmoth
# merge these variables into a single data frame to output

derived.NMRS <- merge(margin.slopes,distrib.slopes)

# write it out
write.csv(derived.NMRS, file = "Data/NMRS/Derived/derived.csv", row.names = F)





### well-recorded hectads ####
# we now want to repeat all of the above using only data from the well-recorded hectads

# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
NMRS_GB_well.inters <- NMRS_GB_well[which(NMRS_GB_well$YEAR %in% years), ]


# get rid of the year information, keeping only which hectads each species was recorded in:
NMRS_GB_well.hecs <- ddply(NMRS_GB_well.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                           PRESENCE = 1)

# and now count the number of hectads per species
species.hecs.well <- ddply(NMRS_GB_well.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                      HECTADS = sum(PRESENCE))

print(species.hecs.well)

# we need a list of species with fewer than 20 hectads of records:
under.recorded.species.well <- species.hecs.well[which(species.hecs.well$HECTADS < 20), ]

print(under.recorded.species.well)

under.recorded.well <- levels(droplevels(under.recorded.species.well$COMMON_NAME))

# this chops out a bunch of things - most either very rare (e.g. Barberry Carpet, Fisher's Estuarine)
# or migrant (e.g. Spurge Hawk-moth)

# we also want to manually remove migrants and aggregates as before

spec.to.cut.well <- NULL

for (x in c(under.recorded.well,migrants,aggregates,pheromones)){
  if (!(x %in% spec.to.cut.well)){
    spec.to.cut.well <- append(spec.to.cut.well,x)
  }
}


# next we want to focus on southerly species: i.e. those that reach a northern range margin in the UK, are likely climate-limited here,
# and have the capacity to expand northwards

# the rule for this is that the species' northern range margin across the entire period can't be within 100 km of John O'Groats
# so obviously for this, we need to estimate the northern range margins!

# these are the mean northing of the 10 most northerly occupied hectads in the interval

# calculate the margins for each species
margins.well <- data.frame(SCI_NAME = factor(),
                           NORTHING = numeric(),
                           COMMON_NAME = factor())

species.list.well <- levels(droplevels(species.hecs.well$COMMON_NAME))
species.list.well <- species.list.well[which(!(species.list.well %in% spec.to.cut.well))]

for (x in species.list.well){
  sp.dat <- NMRS_GB_well.hecs[which(NMRS_GB_well.hecs$COMMON_NAME==x), ]
  northern10 <- top_n(sp.dat, 10, sp.dat$NORTHING)
  species.rm <- ddply(northern10, .(COMMON_NAME), numcolwise(mean))
  species.rm <- species.rm[,c(1,3)]
  colnames(species.rm) <- c("COMMON_NAME","NORTHING")
  
  species.rm$SCI_NAME <- sp.dat[[1,2]]
  
  
  margins.well <- rbind(margins.well,species.rm)
}

summary(margins.well)


## now we can calculate which species to drop:
# the hectad containing John O'Groats is ND37, so the northing is 970000
# 100 km south of here would be 870000 - around Lossiemouth, just north of Inverness
# we want to check whether each species' NRM is south of here in BOTH intervals

margins.well$EXCLUDE <- as.factor(ifelse(margins.well$NORTHING < 870000, "Southerly","Ubiquitous"))
summary(margins.well$EXCLUDE)

print(margins.well)

# here we lose quite a range of species with distributions that extend into Scotland
# some are truly ubiquitous -
# e.g. Common Blue, Green-veined White and Meadow Brown
# some have patchy distributions but happen to be present in northern Scotland -
# e.g. Grayling, Small Blue
# and some are northerly-distributed
# e.g. Large Heath, Scotch Argus

# none of these are climate-limited (southerly-distributed) so we're fine to get rid of all of them 
# however we *will* want the data on them, so don't append them to the spec.to.cut list

margins.drop.well <- margins.well[which(margins.well$EXCLUDE == "Ubiquitous"), ]

margins.spec.drop.well <- levels(droplevels(margins.drop.well$COMMON_NAME))


## elevational species assessment used ALL recorded hectads, so doesn't need repeating...
# we have our list of mean elevations at hectad level
# now we want to drop out anything where the mean elevation is >200m
# and add these to the list of species to cut
for (x in elevs.drop){
  if(!(x %in% spec.to.cut.well)){
    spec.to.cut.well <- append(spec.to.cut.well,x)
  }
}


# finally, we only want species that have been recorded in a minimum proportion of the 20 years
# count how many years each species has a record (in any hectad) for 

NMRS_GB_well.years <- ddply(NMRS_GB_well.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

NMRS_GB_well.years.count <- ddply(NMRS_GB_well.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(NMRS_GB_well.years.count)

# now retain only species recorded in at least 15/20 years
insufficient.years.spec.well <- NMRS_GB_well.years.count[which(NMRS_GB_well.years.count$YEARS < 20), ]

insufficient.years.well <- levels(droplevels(insufficient.years.spec.well$COMMON_NAME))

for (x in insufficient.years.well){
  if(!(x %in% spec.to.cut.well)){
    spec.to.cut.well <- append(spec.to.cut.well,x)
  }
}




# we finally have a final list of species to cut - so we can derive a list of species to use!

spec.to.use.well <- NULL

for (x in species.list.well){
  if(!(x %in% spec.to.cut.well)){
    spec.to.use.well <- append(spec.to.use.well,x)
  }
}


### variable extraction ####

# for each species left, we want to extract the distribution change and northern range margin advance from this data
# and the abundance change and phenology advance from the UKBMS data

# since we have the data read in already, let's start with the distribution change
# we need to calculate the percentage of recorded hectads (from those chosen for retention based on recording level) that were occupied in each year
# and plot a slope through time to establish whether this has (significantly) changed

## select out data from years of interest
hectads.years.well <- NMRS_GB_well[which(NMRS_GB_well$YEAR %in% years), ]

# calculate how many hectads were recorded in each year
# first count how many species recorded in each hectad per year
hectads.rec.years.well <- ddply(hectads.years.well, .(HECTAD,YEAR,EASTING,NORTHING), summarise, SPECIES = sum(PRESENCE), PRESENCE = 1)

# then how many hectads in each year
hectads.rec.well <- ddply(hectads.rec.years.well, .(YEAR), summarise, HECTADS = sum(PRESENCE))
colnames(hectads.rec.well)[2] <- "ANNUAL.TOTAL"


# and now pick out species of interest
hectads.years.well <- hectads.years.well[which(hectads.years.well$COMMON_NAME %in% spec.to.use.well),]

# and count the number of hectads in each year for each species
hectads.species.well <- ddply(hectads.years.well, .(SCI_NAME,COMMON_NAME,YEAR), summarise, HECTADS = sum(PRESENCE))

# now merge in the annual total number of hectads
hectads.annualtotals.well <- merge(hectads.species.well,hectads.rec.well)

# and calculate the percentage of all recorded hectads in which each species was recorded each year
hectads.annualtotals.well$PERC.REC <- hectads.annualtotals.well$HECTADS*100/hectads.annualtotals.well$ANNUAL.TOTAL

## we also want a scaled percentage (i.e. PERC compared to PERC in 1995)
# this gives us an estimate of distribution size change relative to the size of the starting distribution

# extract out the 1995 estimates
hectads.1995.well <- hectads.annualtotals.well[which(hectads.annualtotals.well$YEAR == 1995), c(2:3,6)]
colnames(hectads.1995.well)[3] <- "START.PERC.REC"

# merge them back in 
hectads.annualtotals.well <- merge(hectads.annualtotals.well,hectads.1995.well)

# and calculate the scaled distribution change
hectads.annualtotals.well$PERC.CHANGE <- hectads.annualtotals.well$PERC.REC*100/hectads.annualtotals.well$START.PERC.REC

## now we need change over time for each species - this is a linear regression of PERC against YEAR
# let's set it up inside a loop to export slope and P-value, as well as a plot

distrib.slopes.well <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  DISTRIB.SLOPE = numeric(),
                                  DISTRIB.SE = numeric(),
                                  DISTRIB.CHI = numeric(),
                                  DISTRIB.P = numeric())

for (x in spec.to.use.well){
  spec <- hectads.annualtotals.well[which(hectads.annualtotals.well$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(PERC.CHANGE ~ YEAR, data = spec)
  
  DISTRIB.SLOPE <- round(summary(test1)$coefficients[2,1], 2)
  DISTRIB.SE <- round(summary(test1)$coefficients[2,2], 2)
  DISTRIB.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  DISTRIB.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/Distribution/Moths/WellRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PERC.CHANGE ~ spec$YEAR, 
       ylim = c(0,300),
       xlab = "YEAR", ylab = "Annual change in % of recorded hectads occupied (scaled to 1995 level)")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (DISTRIB.P == 0){
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p = ",DISTRIB.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, DISTRIB.SLOPE, DISTRIB.SE, DISTRIB.CHI, DISTRIB.P))
  distrib.slopes.well <- rbind(distrib.slopes.well, out)
  
}

distrib.slopes.well$DISTRIB.SLOPE <- as.numeric(as.character(distrib.slopes.well$DISTRIB.SLOPE))
distrib.slopes.well$DISTRIB.SE <- as.numeric(as.character(distrib.slopes.well$DISTRIB.SE))
distrib.slopes.well$DISTRIB.CHI <- as.numeric(as.character(distrib.slopes.well$DISTRIB.CHI))
distrib.slopes.well$DISTRIB.P <- as.numeric(as.character(distrib.slopes.well$DISTRIB.P))

summary(distrib.slopes.well)



## and now northern range margins
# we want to select out the ten northernmost hectads for each species in each year
# retaining ALL hectads that are tied for tenth (if there are scattered cells above a hard margin, this pulls towards the hard margin)
# and take the mean latitude of these cells


# so, for each year, first select out the northernmost hectads for each species and take their mean
# seed an output table

margins.well <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           YEAR = numeric(),
                           NORTHING = numeric())


for (x in spec.to.use.well){
  spec <- hectads.years.well[which(hectads.years.well$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,2]])
  
  for (n in years){
    YEAR <- n
    spec.year <- spec[which(spec$YEAR==n), ]
    northern10 <- top_n(spec.year, 10, spec.year$NORTHING)
    spec.year.rm <- ddply(northern10, .(SCI_NAME,COMMON_NAME,YEAR), summarise, NORTHING = mean(NORTHING))
    
    NORTHING <- as.character(spec.year.rm[[1,4]])
    
    out <- data.frame(cbind(SCI_NAME,COMMON_NAME,YEAR,NORTHING))
    margins.well <- rbind(margins.well,out)
  }
}

margins.well$YEAR <- as.numeric(as.character(margins.well$YEAR))
margins.well$NORTHING <- as.numeric(as.character(margins.well$NORTHING))

summary(margins.well)


# plot this on those maps as lines over the ranges
for (x in spec.to.use.well){
  spec <- hectads.years.well[which(hectads.years.well$COMMON_NAME==x), ]
  spec.rm <- margins.well[which(margins.well$COMMON_NAME==x), ]
  spec.rm <- merge(spec.rm,cols)
  png(paste0("Data/Derived/RangeMargins/Margins/Moths/WellRecorded/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points((5+spec$EASTING/1000),(5+spec$NORTHING/1000),pch=16, col="black")
  title(main = x)
  for (n in 1995:2014){
    abline(h=(spec.rm[which(spec.rm$YEAR==n),4]/1000), col=as.character(cols[which(cols$YEAR==n),2]), lwd = 2)
  }
  legend("bottomright", col=cols.pal(3), lwd=2,
         legend=c("1995","2005","2014"))
  
  dev.off()
} 



# let's try calculating change in range margins, then:
# here we want the absolute change in km per year

margin.slopes.well <- data.frame(SCI_NAME = factor(),
                                 COMMON_NAME = factor(),
                                 MARGIN.SLOPE = numeric(),
                                 MARGIN.SE = numeric(),
                                 MARGIN.CHI = numeric(),
                                 MARGIN.P = numeric())

for (x in spec.to.use.well){
  spec <- margins.well[which(margins.well$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(NORTHING ~ YEAR, data = spec)
  
  MARGIN.SLOPE <- round((summary(test1)$coefficients[2,1])/1000, 2)
  MARGIN.SE <- round((summary(test1)$coefficients[2,2])/1000, 2)
  MARGIN.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  MARGIN.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/RangeMargins/Change/Moths/WellRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$NORTHING ~ spec$YEAR, 
       ylim = c(70000,1200000),
       xlab = "YEAR", ylab = "Annual Northing of Northern Range Margin")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (MARGIN.P == 0){
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p = ",MARGIN.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, MARGIN.SLOPE, MARGIN.SE, MARGIN.CHI, MARGIN.P))
  margin.slopes.well <- rbind(margin.slopes.well, out)
  
}

margin.slopes.well$MARGIN.SLOPE <- as.numeric(as.character(margin.slopes.well$MARGIN.SLOPE))
margin.slopes.well$MARGIN.SE <- as.numeric(as.character(margin.slopes.well$MARGIN.SE))
margin.slopes.well$MARGIN.CHI <- as.numeric(as.character(margin.slopes.well$MARGIN.CHI))
margin.slopes.well$MARGIN.P <- as.numeric(as.character(margin.slopes.well$MARGIN.P))

summary(margin.slopes.well)



# we have a couple of retractions e.g. Grizzled Skipper where northerly populations have gone extinct
# and some big expansions of which Comma is unsurprisingly the biggest

# merge these variables into a single data frame to output

derived.NMRS.well <- merge(margin.slopes.well,distrib.slopes.well)

# write it out
write.csv(derived.NMRS.well, file = "Data/NMRS/Derived/derived_wellrecorded.csv", row.names = F)





### heavily recorded hectads ####
# we now want to repeat all of the above using only data from the heavy-recorded hectads

# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
NMRS_GB_heavy.inters <- NMRS_GB_heavy[which(NMRS_GB_heavy$YEAR %in% years), ]


# get rid of the year information, keeping only which hectads each species was recorded in:
NMRS_GB_heavy.hecs <- ddply(NMRS_GB_heavy.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                           PRESENCE = 1)

# and now count the number of hectads per species
species.hecs.heavy <- ddply(NMRS_GB_heavy.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                           HECTADS = sum(PRESENCE))

print(species.hecs.heavy)

# we need a list of species with fewer than 20 hectads of records:
under.recorded.species.heavy <- species.hecs.heavy[which(species.hecs.heavy$HECTADS < 20), ]

print(under.recorded.species.heavy)

under.recorded.heavy <- levels(droplevels(under.recorded.species.heavy$COMMON_NAME))

# this chops out a bunch of things - most either very rare (e.g. Barberry Carpet, Fisher's Estuarine)
# or migrant (e.g. Spurge Hawk-moth)

# we also want to manually remove migrants and aggregates as before

spec.to.cut.heavy <- NULL

for (x in c(under.recorded.heavy,migrants,aggregates,pheromones)){
  if (!(x %in% spec.to.cut.heavy)){
    spec.to.cut.heavy <- append(spec.to.cut.heavy,x)
  }
}


# next we want to focus on southerly species: i.e. those that reach a northern range margin in the UK, are likely climate-limited here,
# and have the capacity to expand northwards

# the rule for this is that the species' northern range margin across the entire period can't be within 100 km of John O'Groats
# so obviously for this, we need to estimate the northern range margins!

# these are the mean northing of the 10 most northerly occupied hectads in the interval

# calculate the margins for each species
margins.heavy <- data.frame(SCI_NAME = factor(),
                           NORTHING = numeric(),
                           COMMON_NAME = factor())

species.list.heavy <- levels(droplevels(species.hecs.heavy$COMMON_NAME))
species.list.heavy <- species.list.heavy[which(!(species.list.heavy %in% spec.to.cut.heavy))]

for (x in species.list.heavy){
  sp.dat <- NMRS_GB_heavy.hecs[which(NMRS_GB_heavy.hecs$COMMON_NAME==x), ]
  northern10 <- top_n(sp.dat, 10, sp.dat$NORTHING)
  species.rm <- ddply(northern10, .(COMMON_NAME), numcolwise(mean))
  species.rm <- species.rm[,c(1,3)]
  colnames(species.rm) <- c("COMMON_NAME","NORTHING")
  
  species.rm$SCI_NAME <- sp.dat[[1,2]]
  
  
  margins.heavy <- rbind(margins.heavy,species.rm)
}

summary(margins.heavy)


## now we can calculate which species to drop:
# the hectad containing John O'Groats is ND37, so the northing is 970000
# 100 km south of here would be 870000 - around Lossiemouth, just north of Inverness
# we want to check whether each species' NRM is south of here in BOTH intervals

margins.heavy$EXCLUDE <- as.factor(ifelse(margins.heavy$NORTHING < 870000, "Southerly","Ubiquitous"))
summary(margins.heavy$EXCLUDE)

print(margins.heavy)

# here we lose quite a range of species with distributions that extend into Scotland
# however we *will* want the data on them, so don't append them to the spec.to.cut list

margins.drop.heavy <- margins.heavy[which(margins.heavy$EXCLUDE == "Ubiquitous"), ]

margins.spec.drop.heavy <- levels(droplevels(margins.drop.heavy$COMMON_NAME))


## elevational species assessment used ALL recorded hectads, so doesn't need repeating...
# we have our list of mean elevations at hectad level
# now we want to drop out anything where the mean elevation is >200m
# and add these to the list of species to cut
for (x in elevs.drop){
  if(!(x %in% spec.to.cut.heavy)){
    spec.to.cut.heavy <- append(spec.to.cut.heavy,x)
  }
}


# finally, we only want species that have been recorded in a minimum proportion of the 20 years
# count how many years each species has a record (in any hectad) for 

NMRS_GB_heavy.years <- ddply(NMRS_GB_heavy.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

NMRS_GB_heavy.years.count <- ddply(NMRS_GB_heavy.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(NMRS_GB_heavy.years.count)

# now retain only species recorded in at least 15/20 years
insufficient.years.spec.heavy <- NMRS_GB_heavy.years.count[which(NMRS_GB_heavy.years.count$YEARS < 20), ]

insufficient.years.heavy <- levels(droplevels(insufficient.years.spec.heavy$COMMON_NAME))

for (x in insufficient.years.heavy){
  if(!(x %in% spec.to.cut.heavy)){
    spec.to.cut.heavy <- append(spec.to.cut.heavy,x)
  }
}




# we finally have a final list of species to cut - so we can derive a list of species to use!

spec.to.use.heavy <- NULL

for (x in species.list.heavy){
  if(!(x %in% spec.to.cut.heavy)){
    spec.to.use.heavy <- append(spec.to.use.heavy,x)
  }
}


### variable extraction ####

# for each species left, we want to extract the distribution change and northern range margin advance from this data
# and the abundance change and phenology advance from the UKBMS data

# since we have the data read in already, let's start with the distribution change
# we need to calculate the percentage of recorded hectads (from those chosen for retention based on recording level) that were occupied in each year
# and plot a slope through time to establish whether this has (significantly) changed

## select out data from years of interest
hectads.years.heavy <- NMRS_GB_heavy[which(NMRS_GB_heavy$YEAR %in% years), ]

# calculate how many hectads were recorded in each year
# first count how many species recorded in each hectad per year
hectads.rec.years.heavy <- ddply(hectads.years.heavy, .(HECTAD,YEAR,EASTING,NORTHING), summarise, SPECIES = sum(PRESENCE), PRESENCE = 1)

# then how many hectads in each year
hectads.rec.heavy <- ddply(hectads.rec.years.heavy, .(YEAR), summarise, HECTADS = sum(PRESENCE))
colnames(hectads.rec.heavy)[2] <- "ANNUAL.TOTAL"


# and now pick out species of interest
hectads.years.heavy <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME %in% spec.to.use.heavy),]

# and count the number of hectads in each year for each species
hectads.species.heavy <- ddply(hectads.years.heavy, .(SCI_NAME,COMMON_NAME,YEAR), summarise, HECTADS = sum(PRESENCE))

# now merge in the annual total number of hectads
hectads.annualtotals.heavy <- merge(hectads.species.heavy,hectads.rec.heavy)

# and calculate the percentage of all recorded hectads in which each species was recorded each year
hectads.annualtotals.heavy$PERC.REC <- hectads.annualtotals.heavy$HECTADS*100/hectads.annualtotals.heavy$ANNUAL.TOTAL

### at this point we want to quickly calculate the total number of records (i.e. hectads x years x species)
sum(hectads.annualtotals.heavy$HECTADS)

## and generate species-level summaries
NMRS.summary.heavy <- data.frame(SCI_NAME = factor(),
                                 COMMON_NAME = factor(),
                                 HECTAD.RECORDS = numeric(),    # i.e. total number of records summed across all years
                                 RECORDED.HECTADS = numeric())  # i.e. total number of hectads recorded in any year

for (x in spec.to.use.heavy){
  spec.anntot <- hectads.annualtotals.heavy[which(hectads.annualtotals.heavy$COMMON_NAME==x), ]
  
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec.anntot[[1,2]])
  
  HECTAD.RECORDS <- sum(spec.anntot$HECTADS)
  
  spec <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME==x), ]
  spec.coll <- ddply(spec, .(HECTAD,SCI_NAME,COMMON_NAME), summarise,
                     YEARS = sum(PRESENCE))
  
  RECORDED.HECTADS <- nrow(spec.coll)
  
  out <- data.frame(cbind(SCI_NAME,COMMON_NAME,HECTAD.RECORDS,RECORDED.HECTADS))
  NMRS.summary.heavy <- rbind(NMRS.summary.heavy,out)
}


summary(NMRS.summary.heavy)


## we also want a scaled percentage (i.e. PERC compared to PERC in 1995)
# this gives us an estimate of distribution size change relative to the size of the starting distribution

# extract out the 1995 estimates
hectads.1995.heavy <- hectads.annualtotals.heavy[which(hectads.annualtotals.heavy$YEAR == 1995), c(2:3,6)]
colnames(hectads.1995.heavy)[3] <- "START.PERC.REC"

# merge them back in 
hectads.annualtotals.heavy <- merge(hectads.annualtotals.heavy,hectads.1995.heavy)

# and calculate the scaled distribution change
hectads.annualtotals.heavy$PERC.CHANGE <- hectads.annualtotals.heavy$PERC.REC*100/hectads.annualtotals.heavy$START.PERC.REC

## now we need change over time for each species - this is a linear regression of PERC against YEAR
# let's set it up inside a loop to export slope and P-value, as heavy as a plot

distrib.slopes.heavy <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  DISTRIB.SLOPE = numeric(),
                                  DISTRIB.SE = numeric(),
                                  DISTRIB.CHI = numeric(),
                                  DISTRIB.P = numeric())

for (x in spec.to.use.heavy){
  spec <- hectads.annualtotals.heavy[which(hectads.annualtotals.heavy$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(PERC.CHANGE ~ YEAR, data = spec)
  
  DISTRIB.SLOPE <- round(summary(test1)$coefficients[2,1], 2)
  DISTRIB.SE <- round(summary(test1)$coefficients[2,2], 2)
  DISTRIB.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  DISTRIB.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/Distribution/Moths/HeavyRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PERC.CHANGE ~ spec$YEAR, 
       ylim = c(0,300),
       xlab = "YEAR", ylab = "Annual change in % of recorded hectads occupied (scaled to 1995 level)")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (DISTRIB.P == 0){
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",DISTRIB.SLOPE," +/- ",DISTRIB.SE,"; chisq = ",DISTRIB.CHI,"; p = ",DISTRIB.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, DISTRIB.SLOPE, DISTRIB.SE, DISTRIB.CHI, DISTRIB.P))
  distrib.slopes.heavy <- rbind(distrib.slopes.heavy, out)
  
}

distrib.slopes.heavy$DISTRIB.SLOPE <- as.numeric(as.character(distrib.slopes.heavy$DISTRIB.SLOPE))
distrib.slopes.heavy$DISTRIB.SE <- as.numeric(as.character(distrib.slopes.heavy$DISTRIB.SE))
distrib.slopes.heavy$DISTRIB.CHI <- as.numeric(as.character(distrib.slopes.heavy$DISTRIB.CHI))
distrib.slopes.heavy$DISTRIB.P <- as.numeric(as.character(distrib.slopes.heavy$DISTRIB.P))

summary(distrib.slopes.heavy)



## and now northern range margins
# we want to select out the ten northernmost hectads for each species in each year
# retaining ALL hectads that are tied for tenth (if there are scattered cells above a hard margin, this pulls towards the hard margin)
# and take the mean latitude of these cells


# so, for each year, first select out the northernmost hectads for each species and take their mean
# seed an output table

margins.heavy <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           YEAR = numeric(),
                           NORTHING = numeric())


for (x in spec.to.use.heavy){
  spec <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,2]])
  
  for (n in years){
    YEAR <- n
    spec.year <- spec[which(spec$YEAR==n), ]
    northern10 <- top_n(spec.year, 10, spec.year$NORTHING)
    spec.year.rm <- ddply(northern10, .(SCI_NAME,COMMON_NAME,YEAR), summarise, NORTHING = mean(NORTHING))
    
    NORTHING <- as.character(spec.year.rm[[1,4]])
    
    out <- data.frame(cbind(SCI_NAME,COMMON_NAME,YEAR,NORTHING))
    margins.heavy <- rbind(margins.heavy,out)
  }
}

margins.heavy$YEAR <- as.numeric(as.character(margins.heavy$YEAR))
margins.heavy$NORTHING <- as.numeric(as.character(margins.heavy$NORTHING))

summary(margins.heavy)


# plot this on those maps as lines over the ranges
for (x in spec.to.use.heavy){
  spec <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME==x), ]
  spec.rm <- margins.heavy[which(margins.heavy$COMMON_NAME==x), ]
  spec.rm <- merge(spec.rm,cols)
  png(paste0("Data/Derived/RangeMargins/Margins/Moths/HeavyRecorded/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points((5+spec$EASTING/1000),(5+spec$NORTHING/1000),pch=16, col="black")
  title(main = x)
  for (n in 1995:2014){
    abline(h=(spec.rm[which(spec.rm$YEAR==n),4]/1000), col=as.character(cols[which(cols$YEAR==n),2]), lwd = 2)
  }
  legend("bottomright", col=cols.pal(3), lwd=2,
         legend=c("1995","2005","2014"))
  
  dev.off()
} 



# let's try calculating change in range margins, then:
# here we want the absolute change in km per year

margin.slopes.heavy <- data.frame(SCI_NAME = factor(),
                                 COMMON_NAME = factor(),
                                 MARGIN.SLOPE = numeric(),
                                 MARGIN.SE = numeric(),
                                 MARGIN.CHI = numeric(),
                                 MARGIN.P = numeric())

for (x in spec.to.use.heavy){
  spec <- margins.heavy[which(margins.heavy$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec[[1,1]])
  
  test1 <- lm(NORTHING ~ YEAR, data = spec)
  
  MARGIN.SLOPE <- round((summary(test1)$coefficients[2,1])/1000, 2)
  MARGIN.SE <- round((summary(test1)$coefficients[2,2])/1000, 2)
  MARGIN.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
  MARGIN.P <- round(drop1(test1, test = "F")[2,6], 4)
  
  png(paste0("Data/Derived/RangeMargins/Change/Moths/HeavyRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$NORTHING ~ spec$YEAR, 
       ylim = c(70000,1200000),
       xlab = "YEAR", ylab = "Annual Northing of Northern Range Margin")  
  abline(test1)
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (MARGIN.P == 0){
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",MARGIN.SLOPE," +/- ",MARGIN.SE,"; chisq = ",MARGIN.CHI,"; p = ",MARGIN.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, MARGIN.SLOPE, MARGIN.SE, MARGIN.CHI, MARGIN.P))
  margin.slopes.heavy <- rbind(margin.slopes.heavy, out)
  
}

margin.slopes.heavy$MARGIN.SLOPE <- as.numeric(as.character(margin.slopes.heavy$MARGIN.SLOPE))
margin.slopes.heavy$MARGIN.SE <- as.numeric(as.character(margin.slopes.heavy$MARGIN.SE))
margin.slopes.heavy$MARGIN.CHI <- as.numeric(as.character(margin.slopes.heavy$MARGIN.CHI))
margin.slopes.heavy$MARGIN.P <- as.numeric(as.character(margin.slopes.heavy$MARGIN.P))

summary(margin.slopes.heavy)



# merge these variables into a single data frame to output

derived.NMRS.heavy.pre <- merge(NMRS.summary.heavy,margin.slopes.heavy)
derived.NMRS.heavy <- merge(derived.NMRS.heavy.pre,distrib.slopes.heavy)

# write it out
write.csv(derived.NMRS.heavy, file = "Data/NMRS/Derived/derived_heavyrecorded.csv", row.names = F)






### RIS data ####
# now let's pull in the RIS data and see how much of that passes its quality filters - and what species remain when we combine both datasets

# read in the data
RIS_raw <- read.csv("../../Data/RIS/GAMphenology_refined.csv", header = T)

summary(RIS_raw)

# we have already applied a couple of the QC filters in generating this file -
# we only have rows of data here where a transect was recorded 10 times in a year,
# and the species of interest was recorded on 3 of those
# we also have "Fail" entered under PEAKDAY wherever the GAM either doesn't start the year at zero, 
# or fails to peak and decrease before the end of the year
# and we have quite a few instances of "Fail" under the voltinism ratio, which we don't want to chuck

# we can dispose of the 'fails' now, and make the response variable numeric again
RIS_raw <- RIS_raw[which(RIS_raw$PEAKDAY != "Fail"), ]
RIS_raw$PEAKDAY <- as.numeric(as.character(RIS_raw$PEAKDAY)) 

# the remaining QC filters are:
# species present in every year (most years, if too harsh?) during the study period 1995-2014
# species present at some point in  1980-1990 (so no colonisation effects by first window)

# this creates two time windows, with slightly different rules for the first
# we need to separately check whether each population (i.e. site x species combination) passes the filters in each time window
# and then select out all the populations which pass the filters in both windows to use

## prequel time window - species present by 1981
# set up the time window
prequel <- 1980:1990

# select out the data
RIS_prequel <- RIS_raw[which(RIS_raw$YEAR %in% prequel), ]

# collapse down to a list of populations, counting how many years each was recorded in
RIS_prequel_popns <- ddply(RIS_prequel, .(siteXspecies,SCI_NAME,SPECIES,COMMON_NAME,SITE), summarise, YEARS.PREQUEL = length(PEAKDAY))


## interval: present in all years 1995-2014

# select out the data
years <- 1995:2014
RIS_years <- RIS_raw[which(RIS_raw$YEAR %in% years), ]

# collapse down, counting years
RIS_years_popns <- ddply(RIS_years, .(siteXspecies,SCI_NAME,SPECIES,COMMON_NAME,SITE), summarise, YEARS = length(PEAKDAY))


# try a few different threshold levels -
# recorded in *every* year
RIS_years_t20 <- RIS_years_popns[which(RIS_years_popns$YEARS>=20), ]

# recorded in at least 18/20 years
RIS_years_t18 <- RIS_years_popns[which(RIS_years_popns$YEARS>=18), ]

# recorded in at least 15/20 years
RIS_years_t15 <- RIS_years_popns[which(RIS_years_popns$YEARS>=15), ]

# recorded in at least 10/20 years
RIS_years_t10 <- RIS_years_popns[which(RIS_years_popns$YEARS>=10), ]




## try merging these together to produce candidate datasets at each threshold

RIS_t20 <- merge(RIS_years_t20, RIS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
RIS_t20_species <- ddply(RIS_t20, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

RIS_t20_good <- RIS_t20_species[which(RIS_t20_species$POPULATIONS >= 3), ]

RIS_t20_use <- RIS_t20_good[which(RIS_t20_good$COMMON_NAME %in% spec.to.use), ]

# setting the threshold at 20 (i.e. must have been recorded in *every* year) leaves very few populations behind,
# and most of those are populations of ubiquitous species that won't be included in the NRM analysis
# this nonetheless leaves us with a fairly respectable 41 species using this approach

# ideally we want to retain as many of the good NMRS species as possible
# therefore we'll need to use a less strict threshold - let's try some of the other options

# t18
RIS_t18 <- merge(RIS_years_t18,RIS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
RIS_t18_species <- ddply(RIS_t18, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

RIS_t18_good <- RIS_t18_species[which(RIS_t18_species$POPULATIONS >= 3), ]

RIS_t18_use <- RIS_t18_good[which(RIS_t18_good$COMMON_NAME %in% spec.to.use), ]

# this immediately increase the number of species to 69!
# let's delve even further...

# t15
RIS_t15 <- merge(RIS_years_t15,RIS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
RIS_t15_species <- ddply(RIS_t15, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

RIS_t15_good <- RIS_t15_species[which(RIS_t15_species$POPULATIONS >= 3), ]

RIS_t15_use <- RIS_t15_good[which(RIS_t15_good$COMMON_NAME %in% spec.to.use), ]

# now up to 104 species...

# t10
RIS_t10 <- merge(RIS_years_t10,RIS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
RIS_t10_species <- ddply(RIS_t10, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

RIS_t10_good <- RIS_t10_species[which(RIS_t10_species$POPULATIONS >= 3), ]

RIS_t10_use <- RIS_t10_good[which(RIS_t10_good$COMMON_NAME %in% spec.to.use), ]

# even dropping it this far adds even more species to the dataset - up to 149



# summarise which species to use/exclude at all possible levels

moths <- data.frame(COMMON_NAME = levels(droplevels(NMRS_raw$COMMON_NAME)))

moths$NMRS.QC <- ifelse(moths$COMMON_NAME %in% spec.to.use, T,F)
moths$RIS.QC.t20 <- ifelse(moths$COMMON_NAME %in% RIS_t20_good$COMMON_NAME, T,F)
moths$RIS.QC.t18 <- ifelse(moths$COMMON_NAME %in% RIS_t18_good$COMMON_NAME, T,F)
moths$RIS.QC.t15 <- ifelse(moths$COMMON_NAME %in% RIS_t15_good$COMMON_NAME, T,F)
moths$RIS.QC.t10 <- ifelse(moths$COMMON_NAME %in% RIS_t10_good$COMMON_NAME, T,F)


moths$USE <- as.factor(ifelse(moths$NMRS.QC == T & moths$RIS.QC.t20 == T, "Use (threshold 20)",
                                    ifelse(moths$NMRS.QC == T & moths$RIS.QC.t20 == F & moths$RIS.QC.t18 == T, "Use? (threshold 18)",
                                           ifelse(moths$NMRS.QC == T & moths$RIS.QC.t18 == F & moths$RIS.QC.t15 == T, "Use? (threshold 15)",
                                                  ifelse(moths$NMRS.QC == T & moths$RIS.QC.t15 == F & moths$RIS.QC.t10 == T, "Exclude? (threshold 10)",
                                                         ifelse(moths$NMRS.QC == F & moths$COMMON_NAME %in% migrants, "Exclude (migrant)",
                                                                ifelse(moths$NMRS.QC == F & moths$COMMON_NAME %in% under.recorded, "Exclude (small distribution)",
                                                                       ifelse(moths$NMRS.QC == F & moths$COMMON_NAME %in% elevs.drop, "Exclude (upland species)",
                                                                              "Exclude (no RIS data)"))))))))

moths$NORTHERLY <- ifelse(moths$COMMON_NAME %in% margins.spec.drop, T, F)



# write out this list
write.csv(moths, "Data/mothspeciestouse.csv", row.names = F)



### extract data for good RIS populations ####

# we want to use populations for species labelled "Use (threshold 20)", "Use (threshold 18)" and "Use (threshold 15)"

moths.to.use.spec <- moths[which(moths$USE %in% c("Use (threshold 20)","Use? (threshold 18)","Use? (threshold 15)")), ]
moths.to.use <- levels(droplevels(moths.to.use.spec$COMMON_NAME))

# now pull out the good sites for these species
# these are held in RIS_t15 - the lowest acceptable threshold for inclusion

popns.to.use <- RIS_t15[which(RIS_t15$COMMON_NAME %in% moths.to.use), ]

# for each of these populations we want to extract the raw data
# we need a single column for population to do this, which we can create by combining species & site names

popns.to.use$POPULATION <- paste(popns.to.use$COMMON_NAME, popns.to.use$SITE, sep=".")

# now repeat this for the raw data, which is in RIS_raw

RIS_years$POPULATION <- paste(RIS_years$COMMON_NAME, RIS_years$SITE, sep=".")

# and use this to extract only the data from good populations

RIS_popns <- RIS_years[which(RIS_years$POPULATION %in% popns.to.use$POPULATION), ]

## now we want a final cleaning step - we have taken forward species with at least 3 populations, where each population is recorded in at least 15 years
# but we additionally want to take forward only species that were recorded *somewhere* in all 20 years

RIS_spec.years <- ddply(RIS_popns, .(YEAR,SCI_NAME,COMMON_NAME), summarise, POPULATIONS = length(PEAKDAY))

RIS_specs <- ddply(RIS_spec.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = length(POPULATIONS))

RIS_specs.good <- RIS_specs[which(RIS_specs$YEARS==20),]
RIS.moths <- levels(droplevels(RIS_specs.good$COMMON_NAME))

# we only lose 4 species by this step: Angle Shades, Common Wainscot, Small Blood-vein, and Shaded Broad-bar
# we can absorb those losses in order to have the best quality of data

RIS_popns <- RIS_popns[which(RIS_popns$COMMON_NAME %in% RIS.moths), ]

# now we have our final dataset, we can start producing outputs
# first, though, we want to extract the average abundance per count (rather than raw annual abundance)
# this is to partly address the possible issue of heavily-monitored transects appearing to have higher abundance than occasionally-monitored ones

RIS_popns$MEAN.ABUND <- RIS_popns$COUNT/RIS_popns$RECS

# let's generate some summary stats about how many RIS populations and how many population*year records go into the final dataset
# first, the overall stat of how many sites in total are included
length(levels(droplevels(RIS_popns$SITE)))

# and now some species-level stats
summary.stats <- data.frame(SCI_NAME = factor(),
                            COMMON_NAME = factor(),
                            SITES = numeric(),
                            SITE.YEARS = numeric(),
                            RECORDS = numeric(),
                            INDIVIDUALS = numeric())

for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  SCI_NAME <- as.character(spec[[1,3]])
  spec$YEARS <- 1
  by.year <- ddply(spec, .(SCI_NAME,COMMON_NAME,SITE), summarise,
                   SITE.YEARS = sum(YEARS),
                   RECORDS = sum(OBS),
                   INDIVIDUALS = sum(COUNT))
  SITES <- nrow(by.year)
  SITE.YEARS <- sum(by.year$SITE.YEARS)
  RECORDS <- sum(by.year$RECORDS)
  INDIVIDUALS <- sum(by.year$INDIVIDUALS)
  
  out <- data.frame(cbind(SCI_NAME,x,SITES,SITE.YEARS,RECORDS,INDIVIDUALS))
  summary.stats <- rbind(summary.stats,out)
}

colnames(summary.stats) <- c("SCI_NAME","COMMON_NAME","SITES","SITE.YEARS","RECORDS","INDIVIDUALS")

summary.stats$SITES <- as.numeric(as.character(summary.stats$SITES))
summary.stats$SITE.YEARS <- as.numeric(as.character(summary.stats$SITE.YEARS))
summary.stats$RECORDS <- as.numeric(as.character(summary.stats$RECORDS))
summary.stats$INDIVIDUALS <- as.numeric(as.character(summary.stats$INDIVIDUALS))

summary(summary.stats)


# now we can calculate the trends through the emergence dates over time
# this is a slight tweak on what we've done above - rather than getting a single value in each year and putting a simple trend through it,
# let's fit a model to the values for all sites, with site as a random effect
# for this, we want to invert the slope so that this is a measure of phenological *advancement* as opposed to phenological *change*

pheno.slopes <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           PHENO.SLOPE = numeric(),
                           PHENO.SE = numeric(),
                           PHENO.CHI = numeric(),
                           PHENO.P = numeric())


for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Phenology/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PEAKDAY ~ spec$YEAR,
       ylim = c(0,365),
       xlab = "YEAR", ylab = "Annual julian day of peak first-generation emergence")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (PHENO.P == 0){
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p = ",PHENO.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, PHENO.SLOPE, PHENO.SE, PHENO.CHI, PHENO.P))
  pheno.slopes <- rbind(pheno.slopes, out)
  
}

pheno.slopes$PHENO.SLOPE <- as.numeric(as.character(pheno.slopes$PHENO.SLOPE))
pheno.slopes$PHENO.SE <- as.numeric(as.character(pheno.slopes$PHENO.SE))
pheno.slopes$PHENO.CHI <- as.numeric(as.character(pheno.slopes$PHENO.CHI))
pheno.slopes$PHENO.P <- as.numeric(as.character(pheno.slopes$PHENO.P))

summary(pheno.slopes)

## once again, we have a pattern where most species are getting earlier - three-quarters of a day per year for Setaceous Hebrew Character
# and a few are doing the opposite (notably Grey Pine Carpet)


### now we want to replicate this process for abundance data from the same set of RIS populations
# (we will have good abundance )
# this data is already read in - 
summary(RIS_popns$MEAN.ABUND)

abund.slopes <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           ABUND.EXP.SLOPE = numeric(),
                           ABUND.SLOPE = numeric(),
                           ABUND.SE = numeric(),
                           ABUND.CHI = numeric(),
                           ABUND.P = numeric())


for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = spec)
  
  ABUND.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
  ABUND.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
  ABUND.SE <- round(summary(test1)$coefficients[2,2], 4)
  ABUND.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  ABUND.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Abundance/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$MEAN.ABUND) ~ spec$YEAR,
       ylim = c(0,10),
       xlab = "YEAR", ylab = "Log (mean no. recorded individuals per trap)")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (ABUND.P == 0){
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"); chisq = ",ABUND.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"; chisq = ",ABUND.CHI,"; p = ",ABUND.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, ABUND.EXP.SLOPE, ABUND.SLOPE, ABUND.SE, ABUND.CHI, ABUND.P))
  abund.slopes <- rbind(abund.slopes, out)
  
}

abund.slopes$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes$ABUND.EXP.SLOPE))
abund.slopes$ABUND.SLOPE <- as.numeric(as.character(abund.slopes$ABUND.SLOPE))
abund.slopes$ABUND.SE <- as.numeric(as.character(abund.slopes$ABUND.SE))
abund.slopes$ABUND.CHI <- as.numeric(as.character(abund.slopes$ABUND.CHI))
abund.slopes$ABUND.P <- as.numeric(as.character(abund.slopes$ABUND.P))

summary(abund.slopes)

## we have a range from species in sharp decline to strongly-increasing species

## finally, do the same for the voltinism ratio
# here we have a much-reduced set of data because (a) there were more failures (where the GAM failed to return to 0 by the end of the year),
# and (b) there are loads of zeroes (which are univoltine emergences)
# we clearly want to chop out (a), and we also want rid of (b) because we are specifically interested in the relationship between gen 1 and gen 2+
# for some univoltine species this will leave little or no data, but that's fine as we're only going to analyse multivoltines anyway for this one
# however, we want to trim out things with much too little data as they just won't fit a line properly - whilst also retaining all true multivoltines

summary(RIS_popns$VOLTINISM.RATIO) # at the moment the failures are in and it's being treated as a factor!

RIS_volti_pre <- RIS_popns[which(RIS_popns$VOLTINISM.RATIO != "Fail"), ]
RIS_volti_pre$VOLTINISM.RATIO <- as.numeric(as.character(RIS_volti_pre$VOLTINISM.RATIO))
RIS_volti <- RIS_volti_pre[which(RIS_volti_pre$VOLTINISM.RATIO != 0), ]

summary(RIS_volti)

# read in a list of what voltinism things have

voltinism <- read.csv("Data/voltinism.csv", header = T)


# check amount of data going in to this analysis per species, compared to the total amount
volti.spec <- data.frame(cbind(summary(RIS_volti$COMMON_NAME),summary(RIS_popns$COMMON_NAME)))
colnames(volti.spec) <- c("VOLTI_RECS","ALL_RECS")
volti.spec$ratio <- volti.spec$VOLTI_RECS/volti.spec$ALL_RECS
volti.spec$COMMON_NAME <- rownames(volti.spec)

volti.spec <- merge(voltinism[,c(3,5)],volti.spec)

# we can see there's not much attrition for true multivoltines but very heavy attrition for true univoltines
# it's not as clear as for butterflies, but we can get most of the multivoltines and only a few univoltines with a cutoff at 10%

volti.spec <- volti.spec[which(volti.spec$ratio > 0.1), ]

RIS.moths.volti <- as.character(volti.spec$COMMON_NAME)

# prep the analysis
# the data are very clearly overdispersed (log-normal)...
hist(RIS_volti$VOLTINISM.RATIO)
hist(log(RIS_volti$VOLTINISM.RATIO))


volti.slopes <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           VOLTI.EXP.SLOPE = numeric(),
                           VOLTI.SLOPE = numeric(),
                           VOLTI.SE = numeric(),
                           VOLTI.CHI = numeric(),
                           VOLTI.P = numeric())


for (x in RIS.moths.volti){
  spec <- RIS_volti[which(RIS_volti$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(log(VOLTINISM.RATIO) ~ YEAR + (1|SITE), data = spec)
  
  VOLTI.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
  VOLTI.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
  VOLTI.SE <- round(summary(test1)$coefficients[2,2], 4)
  VOLTI.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  VOLTI.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Voltinism/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$VOLTINISM.RATIO) ~ spec$YEAR,
       ylim = c(0,10),
       xlab = "YEAR", ylab = "Log (voltinism ratio)")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (ABUND.P == 0){
    mtext(paste0("effect: ",VOLTI.EXP.SLOPE,"(",VOLTI.SLOPE," +/- ",VOLTI.SE,"); chisq = ",VOLTI.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",VOLTI.EXP.SLOPE,"(",VOLTI.SLOPE," +/- ",VOLTI.SE,"; chisq = ",VOLTI.CHI,"; p = ",VOLTI.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, VOLTI.EXP.SLOPE, VOLTI.SLOPE, VOLTI.SE, VOLTI.CHI, VOLTI.P))
  volti.slopes <- rbind(volti.slopes, out)
  
}

volti.slopes$VOLTI.EXP.SLOPE <- as.numeric(as.character(volti.slopes$VOLTI.EXP.SLOPE))
volti.slopes$VOLTI.SLOPE <- as.numeric(as.character(volti.slopes$VOLTI.SLOPE))
volti.slopes$VOLTI.SE <- as.numeric(as.character(volti.slopes$VOLTI.SE))
volti.slopes$VOLTI.CHI <- as.numeric(as.character(volti.slopes$VOLTI.CHI))
volti.slopes$VOLTI.P <- as.numeric(as.character(volti.slopes$VOLTI.P))

summary(volti.slopes)


# perhaps surprisingly, most species have got *less backloaded* rather than more?


# merge all this together

RIS_out_pre <- merge(pheno.slopes,abund.slopes)
RIS_out_sum <- merge(summary.stats,RIS_out_pre)
RIS_out <- merge(RIS_out_sum,volti.slopes, all=T)


# and merge this into the NMRS derived datasets
# some species might have different sci names between the two datasets
# the NMRS one is more up-to-date so ditch the sci name column from the RIS dataset before merging
RIS_NMRS <- merge(RIS_out[,-1], derived.NMRS)
RIS_NMRS_well <- merge(RIS_out[,-1], derived.NMRS.well)
RIS_NMRS_heavy <- merge(RIS_out[,-1], derived.NMRS.heavy)


# and write out all versions
write.csv(RIS_NMRS, "Data/Derived/RIS_NMRS.csv", row.names = F)
write.csv(RIS_NMRS_well, "Data/Derived/RIS_NMRS_well.csv", row.names = F)
write.csv(RIS_NMRS_heavy, "Data/Derived/RIS_NMRS_heavy.csv", row.names = F)



### intra-specific ####
# we want to extract a site-level trend in abundance and phenology for every combination of site * species

# first let's do phenology
pheno.slopes.spec <- data.frame(SCI_NAME = factor(),
                                COMMON_NAME = factor(),
                                SITE = factor(),
                                PHENO.SLOPE = numeric(),
                                PHENO.SE = numeric(),
                                PHENO.CHI = numeric(),
                                PHENO.P = numeric())

for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  sites.vec <- levels(droplevels(spec$SITE))
  
  for (y in sites.vec){
    site <- spec[which(spec$SITE==y), ]
    SITE <- y
    
    test1 <- lm(PEAKDAY~ YEAR, data = site)
    
    PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
    PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
    PHENO.CHI <- round(drop1(test1, test = "F")[2,5], 2) 
    PHENO.P <- round(drop1(test1, test = "F")[2,6], 4)
    
    ifelse(!dir.exists(paste0("Data/Derived/Intraspecific/Phenology/Moths/",x,"/")), 
           dir.create(paste0("Data/Derived/Intraspecific/Phenology/Moths/",x,"/"), recursive = T),
           FALSE)
    
    png(paste0("Data/Derived/Intraspecific/Phenology/Moths/",x,"/",y,".png"), width = 800, height = 800, units = "px", bg = "white")
    
    plot(site$PEAKDAY ~ site$YEAR,
         ylim = c(0,365),
         xlab = "YEAR", ylab = "Annual julian day of peak first-generation emergence")
    
    abline(test1)
    title(main = paste(x,SCI_NAME,y, sep = " - "))
    if (PHENO.P == 0){
      mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p < 0.0001"), line = 0)
    } else {
      mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p = ",PHENO.P), line = 0)
    }
    
    dev.off()
    
    out <- data.frame(cbind(SCI_NAME, COMMON_NAME, SITE, PHENO.SLOPE, PHENO.SE, PHENO.CHI, PHENO.P))
    pheno.slopes.spec <- rbind(pheno.slopes.spec, out)
    
  }
  
}



pheno.slopes.spec$PHENO.SLOPE <- as.numeric(as.character(pheno.slopes.spec$PHENO.SLOPE))
pheno.slopes.spec$PHENO.SE <- as.numeric(as.character(pheno.slopes.spec$PHENO.SE))
pheno.slopes.spec$PHENO.CHI <- as.numeric(as.character(pheno.slopes.spec$PHENO.CHI))
pheno.slopes.spec$PHENO.P <- as.numeric(as.character(pheno.slopes.spec$PHENO.P))

summary(pheno.slopes.spec)


### now we want to replicate this process for abundance data
abund.slopes.spec <- data.frame(SCI_NAME = factor(),
                                COMMON_NAME = factor(),
                                SITE = factor(),
                                ABUND.EXP.SLOPE = numeric(),
                                ABUND.SLOPE = numeric(),
                                ABUND.SE = numeric(),
                                ABUND.CHI = numeric(),
                                ABUND.P = numeric())

for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  sites.vec <- levels(droplevels(spec$SITE))
  
  for (y in sites.vec){
    site <- spec[which(spec$SITE==y), ]
    SITE <- y
    
    test1 <- lm(log(MEAN.ABUND) ~ YEAR, data = site)
    
    ABUND.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
    ABUND.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
    ABUND.SE <- round(summary(test1)$coefficients[2,2], 4)
    ABUND.CHI <- round(drop1(test1, test = "Chi")[2,4], 2) 
    ABUND.P <- round(drop1(test1, test = "Chi")[2,5], 4)
    
    ifelse(!dir.exists(paste0("Data/Derived/Intraspecific/Abundance/Moths/",x,"/")), 
           dir.create(paste0("Data/Derived/Intraspecific/Abundance/Moths/",x,"/"), recursive = T),
           FALSE)
   
     png(paste0("Data/Derived/Intraspecific/Abundance/Moths/",x,"/",y,".png"), width = 800, height = 800, units = "px", bg = "white")
    
    plot(log(site$MEAN.ABUND) ~ site$YEAR,
         ylim = c(0,10),
         xlab = "YEAR", ylab = "Log (mean no. recorded individuals per trap)")
    
    abline(test1)
    title(main = paste(x,SCI_NAME,y, sep = " - "))
    if (ABUND.P == 0){
      mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"); chisq = ",ABUND.CHI,"; p < 0.0001"), line = 0)
    } else {
      mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"; chisq = ",ABUND.CHI,"; p = ",ABUND.P), line = 0)
    }
    
    dev.off()
    
    out <- data.frame(cbind(SCI_NAME, COMMON_NAME, SITE, ABUND.EXP.SLOPE, ABUND.SLOPE, ABUND.SE, ABUND.CHI, ABUND.P))
    abund.slopes.spec <- rbind(abund.slopes.spec, out)
    
  }
  
}


abund.slopes.spec$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes.spec$ABUND.EXP.SLOPE))
abund.slopes.spec$ABUND.SLOPE <- as.numeric(as.character(abund.slopes.spec$ABUND.SLOPE))
abund.slopes.spec$ABUND.SE <- as.numeric(as.character(abund.slopes.spec$ABUND.SE))
abund.slopes.spec$ABUND.CHI <- as.numeric(as.character(abund.slopes.spec$ABUND.CHI))
abund.slopes.spec$ABUND.P <- as.numeric(as.character(abund.slopes.spec$ABUND.P))

summary(abund.slopes.spec)


## merge these data together and write them out
intraspecific <- merge(pheno.slopes.spec, abund.slopes.spec)

write.csv(intraspecific, "Data/Derived/RIS_intraspec.csv", row.names = F)


#### temporal ####

# we now want to divide up the dataset into two 12-year periods, with a slight overlap (1995-2006 & 2003-2014)
# and see whether the pattern changes between them

# first divide up the dataset
inter1 <- 1995:2006
inter2 <- 2003:2014

RIS_inter1 <- RIS_popns[which(RIS_popns$YEAR %in% inter1), ]
RIS_inter2 <- RIS_popns[which(RIS_popns$YEAR %in% inter2), ]


# let's generate some summary stats about how many RIS populations and how many population*year records go into the final dataset
# interval 1
summary.stats.inter1 <- data.frame(SCI_NAME = factor(),
                                   COMMON_NAME = factor(),
                                   SITES = numeric(),
                                   SITE.YEARS = numeric(),
                                   RECORDS = numeric(),
                                   INDIVIDUALS = numeric())

for (x in RIS.moths){
  spec <- RIS_inter1[which(RIS_inter1$COMMON_NAME==x), ]
  SCI_NAME <- as.character(spec[[1,3]])
  spec$YEARS <- 1
  by.year <- ddply(spec, .(SCI_NAME,COMMON_NAME,SITE), summarise,
                   SITE.YEARS = sum(YEARS),
                   RECORDS = sum(OBS),
                   INDIVIDUALS = sum(COUNT))
  SITES <- nrow(by.year)
  SITE.YEARS <- sum(by.year$SITE.YEARS)
  RECORDS <- sum(by.year$RECORDS)
  INDIVIDUALS <- sum(by.year$INDIVIDUALS)
  
  out <- data.frame(cbind(SCI_NAME,x,SITES,SITE.YEARS,RECORDS,INDIVIDUALS))
  summary.stats.inter1 <- rbind(summary.stats.inter1,out)
}

colnames(summary.stats.inter1) <- c("SCI_NAME","COMMON_NAME","SITES","SITE.YEARS","RECORDS","INDIVIDUALS")

summary.stats.inter1$SITES <- as.numeric(as.character(summary.stats.inter1$SITES))
summary.stats.inter1$SITE.YEARS <- as.numeric(as.character(summary.stats.inter1$SITE.YEARS))
summary.stats.inter1$RECORDS <- as.numeric(as.character(summary.stats.inter1$RECORDS))
summary.stats.inter1$INDIVIDUALS <- as.numeric(as.character(summary.stats.inter1$INDIVIDUALS))

summary(summary.stats.inter1)


# interval 2
summary.stats.inter2 <- data.frame(SCI_NAME = factor(),
                                   COMMON_NAME = factor(),
                                   SITES = numeric(),
                                   SITE.YEARS = numeric(),
                                   RECORDS = numeric(),
                                   INDIVIDUALS = numeric())

for (x in RIS.moths){
  spec <- RIS_inter2[which(RIS_inter2$COMMON_NAME==x), ]
  SCI_NAME <- as.character(spec[[1,3]])
  spec$YEARS <- 1
  by.year <- ddply(spec, .(SCI_NAME,COMMON_NAME,SITE), summarise,
                   SITE.YEARS = sum(YEARS),
                   RECORDS = sum(OBS),
                   INDIVIDUALS = sum(COUNT))
  SITES <- nrow(by.year)
  SITE.YEARS <- sum(by.year$SITE.YEARS)
  RECORDS <- sum(by.year$RECORDS)
  INDIVIDUALS <- sum(by.year$INDIVIDUALS)
  
  out <- data.frame(cbind(SCI_NAME,x,SITES,SITE.YEARS,RECORDS,INDIVIDUALS))
  summary.stats.inter2 <- rbind(summary.stats.inter2,out)
}

colnames(summary.stats.inter2) <- c("SCI_NAME","COMMON_NAME","SITES","SITE.YEARS","RECORDS","INDIVIDUALS")

summary.stats.inter2$SITES <- as.numeric(as.character(summary.stats.inter2$SITES))
summary.stats.inter2$SITE.YEARS <- as.numeric(as.character(summary.stats.inter2$SITE.YEARS))
summary.stats.inter2$RECORDS <- as.numeric(as.character(summary.stats.inter2$RECORDS))
summary.stats.inter2$INDIVIDUALS <- as.numeric(as.character(summary.stats.inter2$INDIVIDUALS))

summary(summary.stats.inter2)


# phenology, interval 1

pheno.slopes.inter1 <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  PHENO.SLOPE = numeric(),
                                  PHENO.SE = numeric(),
                                  PHENO.CHI = numeric(),
                                  PHENO.P = numeric())


for (x in RIS.moths){
  spec <- RIS_inter1[which(RIS_inter1$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval1/Phenology/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PEAKDAY ~ spec$YEAR,
       ylim = c(0,365),
       xlab = "YEAR", ylab = "Annual julian day of peak first-generation emergence")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (PHENO.P == 0){
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p = ",PHENO.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, PHENO.SLOPE, PHENO.SE, PHENO.CHI, PHENO.P))
  pheno.slopes.inter1 <- rbind(pheno.slopes.inter1, out)
  
}

pheno.slopes.inter1$PHENO.SLOPE <- as.numeric(as.character(pheno.slopes.inter1$PHENO.SLOPE))
pheno.slopes.inter1$PHENO.SE <- as.numeric(as.character(pheno.slopes.inter1$PHENO.SE))
pheno.slopes.inter1$PHENO.CHI <- as.numeric(as.character(pheno.slopes.inter1$PHENO.CHI))
pheno.slopes.inter1$PHENO.P <- as.numeric(as.character(pheno.slopes.inter1$PHENO.P))

summary(pheno.slopes.inter1)

# phenology, interval 2

pheno.slopes.inter2 <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  PHENO.SLOPE = numeric(),
                                  PHENO.SE = numeric(),
                                  PHENO.CHI = numeric(),
                                  PHENO.P = numeric())


for (x in RIS.moths){
  spec <- RIS_inter2[which(RIS_inter2$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval2/Phenology/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PEAKDAY ~ spec$YEAR,
       ylim = c(0,365),
       xlab = "YEAR", ylab = "Annual julian day of peak first-generation emergence")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (PHENO.P == 0){
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p = ",PHENO.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, PHENO.SLOPE, PHENO.SE, PHENO.CHI, PHENO.P))
  pheno.slopes.inter2 <- rbind(pheno.slopes.inter2, out)
  
}

pheno.slopes.inter2$PHENO.SLOPE <- as.numeric(as.character(pheno.slopes.inter2$PHENO.SLOPE))
pheno.slopes.inter2$PHENO.SE <- as.numeric(as.character(pheno.slopes.inter2$PHENO.SE))
pheno.slopes.inter2$PHENO.CHI <- as.numeric(as.character(pheno.slopes.inter2$PHENO.CHI))
pheno.slopes.inter2$PHENO.P <- as.numeric(as.character(pheno.slopes.inter2$PHENO.P))

summary(pheno.slopes.inter2)

##




### now we want to replicate this process for abundance data from the same set of RIS populations
# interval 1
abund.slopes.inter1 <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  ABUND.EXP.SLOPE = numeric(),
                                  ABUND.SLOPE = numeric(),
                                  ABUND.SE = numeric(),
                                  ABUND.CHI = numeric(),
                                  ABUND.P = numeric())


for (x in RIS.moths){
  spec <- RIS_inter1[which(RIS_inter1$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = spec)
  
  ABUND.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
  ABUND.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
  ABUND.SE <- round(summary(test1)$coefficients[2,2], 4)
  ABUND.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  ABUND.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval1/Abundance/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$COUNT) ~ spec$YEAR,
       ylim = c(0,10),
       xlab = "YEAR", ylab = "Log (mean no. recorded individuals per trap)")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (ABUND.P == 0){
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"); chisq = ",ABUND.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"; chisq = ",ABUND.CHI,"; p = ",ABUND.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, ABUND.EXP.SLOPE, ABUND.SLOPE, ABUND.SE, ABUND.CHI, ABUND.P))
  abund.slopes.inter1 <- rbind(abund.slopes.inter1, out)
  
}

abund.slopes.inter1$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes.inter1$ABUND.EXP.SLOPE))
abund.slopes.inter1$ABUND.SLOPE <- as.numeric(as.character(abund.slopes.inter1$ABUND.SLOPE))
abund.slopes.inter1$ABUND.SE <- as.numeric(as.character(abund.slopes.inter1$ABUND.SE))
abund.slopes.inter1$ABUND.CHI <- as.numeric(as.character(abund.slopes.inter1$ABUND.CHI))
abund.slopes.inter1$ABUND.P <- as.numeric(as.character(abund.slopes.inter1$ABUND.P))

summary(abund.slopes.inter1)



# interval 2
abund.slopes.inter2 <- data.frame(SCI_NAME = factor(),
                                  COMMON_NAME = factor(),
                                  ABUND.EXP.SLOPE = numeric(),
                                  ABUND.SLOPE = numeric(),
                                  ABUND.SE = numeric(),
                                  ABUND.CHI = numeric(),
                                  ABUND.P = numeric())


for (x in RIS.moths){
  spec <- RIS_inter2[which(RIS_inter2$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = spec)
  
  ABUND.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
  ABUND.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
  ABUND.SE <- round(summary(test1)$coefficients[2,2], 4)
  ABUND.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  ABUND.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval2/Abundance/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$COUNT) ~ spec$YEAR,
       ylim = c(0,10),
       xlab = "YEAR", ylab = "Log (mean no. recorded individuals per trap)")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (ABUND.P == 0){
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"); chisq = ",ABUND.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"; chisq = ",ABUND.CHI,"; p = ",ABUND.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, ABUND.EXP.SLOPE, ABUND.SLOPE, ABUND.SE, ABUND.CHI, ABUND.P))
  abund.slopes.inter2 <- rbind(abund.slopes.inter2, out)
  
}

abund.slopes.inter2$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes.inter2$ABUND.EXP.SLOPE))
abund.slopes.inter2$ABUND.SLOPE <- as.numeric(as.character(abund.slopes.inter2$ABUND.SLOPE))
abund.slopes.inter2$ABUND.SE <- as.numeric(as.character(abund.slopes.inter2$ABUND.SE))
abund.slopes.inter2$ABUND.CHI <- as.numeric(as.character(abund.slopes.inter2$ABUND.CHI))
abund.slopes.inter2$ABUND.P <- as.numeric(as.character(abund.slopes.inter2$ABUND.P))

summary(abund.slopes.inter2)



## now merge these together
base.inter1 <- merge(summary.stats.inter1,pheno.slopes.inter1)
inter1.derived <- merge(base.inter1, abund.slopes.inter1)

base.inter2 <- merge(summary.stats.inter2,pheno.slopes.inter2)
inter2.derived <- merge(base.inter2, abund.slopes.inter2)

## label the interval
inter1.derived$INTERVAL <- 1
inter2.derived$INTERVAL <- 2

## and merge again
inters.derived <- rbind(inter1.derived,inter2.derived)

## and export
write.csv(inters.derived, "Data/Derived/RIS_intervals.csv", row.names = F)


### site details ####
# I want to assess the annual climatic conditions at each of our sites, using GDD (growing degree days)
# so I need to put together a list of the sites that we're actually using, and attach grid references and target dates (for the GDD)

# first gather a list of sites
summary(RIS_popns)

RIS_sites <- ddply(RIS_popns, .(SITE), summarise,
                   SPECIES = nlevels(droplevels(SCI_NAME)))


# and generate a dataframe containing all combinations of site, year and three target dates: meteorological spring, summer and autumn
# these are Julian day 60, 152, 244

GDD_targets <- expand.grid(SITE = factor(levels(droplevels(RIS_sites$SITE))),
                           YEAR = c(1995:2014),
                           JULIAN_DAY = c(60,152,244))

summary(GDD_targets)

# we want to quickly adjust these for leap years

GDD_targets$JULIAN_DAY <- ifelse(GDD_targets$YEAR/4 == round(GDD_targets$YEAR/4), 1+GDD_targets$JULIAN_DAY, GDD_targets$JULIAN_DAY)



# now read in a dataframe containing details of all the sites
site_details <- read.csv("../../Data/RIS/sites_full.csv", header = T)

summary(site_details)

# merge them together
GDD_target_sites <- merge(GDD_targets,site_details, all.x = T)



# we can now write this out

write.csv(GDD_target_sites, "Data/RIS_sites.csv", row.names = F)


# for the purpose of plotting a figure, we want a version of the RIS_sites dataset with voltinism included
# merge the voltinism into the RIS_popns

RIS_popns_volti <- merge(RIS_popns, voltinism[,c(1,3,5)])

RIS_sites_volti <- ddply(RIS_popns_volti, .(SITE,VOLTINISM), summarise,
                         SPECIES = nlevels(droplevels(SCI_NAME)), .drop=FALSE)

RIS_sites_volti <- RIS_sites_volti[which(RIS_sites_volti$SITE %in% levels(droplevels(RIS_popns_volti$SITE))), ]

RIS_sites_uv <- ddply(RIS_sites_volti[which(RIS_sites_volti$VOLTINISM == 1), ], .(SITE), summarise,
                      UNIVOLTINES = sum(SPECIES))

RIS_sites_mv <- ddply(RIS_sites_volti[which(RIS_sites_volti$VOLTINISM == 2), ], .(SITE), summarise,
                      MULTIVOLTINES = sum(SPECIES))

RIS_sites_count <- merge(RIS_sites_uv, RIS_sites_mv)

RIS_sites_count$SPECIES <- RIS_sites_count$UNIVOLTINES + RIS_sites_count$MULTIVOLTINES


# now merge in the site details
RIS_site_richness <- merge(RIS_sites_count, site_details)

# and write it out
write.csv(RIS_site_richness, "Data/RIS_sitecounts.csv", row.names = F)



### yearly analysis ####

# as well as these trends over time, I also want to extract a range of annual site-level variables for analysis;
# we want this to include phenology, overall abundance, and voltinism, all of which are in the raw data, so
# essentially this just means exporting the trimmed version of the dataset containing only good sites (RIS_popns)


summary(RIS_popns)

write.csv(RIS_popns, file = "Data/RIS_annual.csv", row.names = F)






### immediate response ####
# I want to assess whether an early emergence in year t leads to an increase in abundance in t+1
# this will require generation of another derived dataset - for each site * species * year combination,
# we need the emergence date and the abundance change to the next year

combos <- levels(droplevels(RIS_popns$siteXspecies))

indices <- data.frame(COMMON_NAME = factor(),
                      SCI_NAME = factor(),
                      SITE = factor(),
                      siteXspecies = factor(),
                      YEAR = numeric(),
                      PEAKDAY.INDEX = numeric(),
                      ABUNDANCE = numeric(),
                      INDEX = numeric())

previous.species <- "NULL"

for (x in combos){
  combo <- RIS_popns[which(RIS_popns$siteXspecies==x), ]
  siteXspecies <- x
  COMMON_NAME <- as.character(combo$COMMON_NAME[[1]])
  SCI_NAME <- as.character(combo$SCI_NAME[[1]])
  
  if(!(SCI_NAME == previous.species)){
    print(SCI_NAME)
  }
  
  SITE <- as.character(combo$SITE[[1]])
  
  PEAKDAY.MDN <- median(combo$PEAKDAY)
  
  for (n in 1995:2013){
    YEAR <- n
    NEXT <- n+1
    
    ABUNDANCE <- combo[which(combo$YEAR==YEAR),'MEAN.ABUND']
    nextabund <- combo[which(combo$YEAR==NEXT),'MEAN.ABUND']
    
    INDEX <- (nextabund - ABUNDANCE)/ABUNDANCE
    
    if (length(ABUNDANCE)==0){
      ABUNDANCE <- NA
    }
    
    if (length(INDEX)==0){
      INDEX <- NA
    }
    
    
    PEAKDAY <- combo[which(combo$YEAR==YEAR),'PEAKDAY']
    PEAKDAY.INDEX <- PEAKDAY.MDN - PEAKDAY
    
    if (length(PEAKDAY.INDEX)==0){
      PEAKDAY.INDEX <- NA
    }
    
    out <- data.frame(cbind(COMMON_NAME,SCI_NAME,SITE,siteXspecies,YEAR,PEAKDAY.INDEX,ABUNDANCE,INDEX))
    indices <- rbind(indices,out)
    
    previous.species <- SCI_NAME
    
  }
  
}

indices$YEAR <- as.numeric(as.character(indices$YEAR))
indices$PEAKDAY.INDEX <- as.numeric(as.character(indices$PEAKDAY.INDEX))
indices$ABUNDANCE <- as.numeric(as.character(indices$ABUNDANCE))
indices$INDEX <- as.numeric(as.character(indices$INDEX))

summary(indices)

# great! let's write out this data now
write.csv(indices, "Data/Derived/RIS_indices.csv", row.names = F)



#### we want to produce the same thing but for NRM (obviously at a species-level rather than a site-level),
# with a *national* average of abundance and phenology for each year tied to it

## we already have annual NRM estimates held in margins.heavy

summary(margins.heavy)

NRM.indices <- data.frame(COMMON_NAME = factor(),
                          SCI_NAME = factor(),
                          YEAR = numeric(),
                          PEAKDAY.MEAN = numeric(),
                          ABUNDANCE.MEAN = numeric(),
                          NRM = numeric(),
                          NRM.INDEX = numeric())


previous.species <- "NULL"

for (x in RIS.moths){
  spec <- RIS_popns[which(RIS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- as.character(spec$COMMON_NAME[[1]])
  SCI_NAME <- as.character(spec$SCI_NAME[[1]])
  
  NRM.spec <- margins.heavy[which(margins.heavy$COMMON_NAME==x), ]
  
  if(!(COMMON_NAME == previous.species)){
    print(COMMON_NAME)
  }
  
  
  for (n in 1995:2013){
    YEAR <- n
    NEXT <- n+1
    
    spec.year <- spec[which(spec$YEAR==YEAR), ]
    
    ABUNDANCE.MEAN <- mean(spec.year$MEAN.ABUND)
    PEAKDAY.MEAN <- mean(spec.year$PEAKDAY)
    
    
    if (length(ABUNDANCE)==0){
      ABUNDANCE <- NA
    }
    
    if (length(PEAKDAY.MEAN)==0){
      PEAKDAY.MEAN <- NA
    }
    
    NRM <- NRM.spec[which(NRM.spec$YEAR == YEAR),'NORTHING']
    NRM.next <- NRM.spec[which(NRM.spec$YEAR == NEXT),'NORTHING']
    
    NRM.INDEX <- (NRM.next - NRM)/1000
    
    if (length(NRM.INDEX)==0){
      NRM.INDEX <- NA
    }
    
    out <- data.frame(cbind(COMMON_NAME,SCI_NAME,YEAR,PEAKDAY.MEAN,ABUNDANCE.MEAN,NRM,NRM.INDEX))
    NRM.indices <- rbind(NRM.indices,out)
    
    previous.species <- COMMON_NAME
    
  }
  
}

NRM.indices$YEAR <- as.numeric(as.character(NRM.indices$YEAR))
NRM.indices$PEAKDAY.MEAN <- as.numeric(as.character(NRM.indices$PEAKDAY.MEAN))
NRM.indices$ABUNDANCE.MEAN <- as.numeric(as.character(NRM.indices$ABUNDANCE.MEAN))
NRM.indices$NRM <- as.numeric(as.character(NRM.indices$NRM))
NRM.indices$NRM.INDEX <- as.numeric(as.character(NRM.indices$NRM.INDEX))

summary(NRM.indices)

# great! let's write out this data now
write.csv(indices, "Data/Derived/NMRS_NRMindices.csv", row.names = F)




### post reviewers' comments ####

# one of the reviewers wants us to do a supplementary analysis using the full length of the available butterfly data
# i.e. for every population, rather than restricting to 1995-2014, calculate the trends from start to finish of recording
# but there aren't many butterflies and the RIS is almost as long as the BMS, so let's do it for the lot


# we still have a list of all the populations to use, in 'popns.to.use'

# so now we return to the completely raw data, with all years included, to extract the relevant data


# now repeat this for the raw data, which is in RIS_raw

RIS_raw$POPULATION <- paste(RIS_raw$COMMON_NAME, RIS_raw$SITE, sep=".")

# and use this to extract only the data from good populations

RIS_allyears <- RIS_raw[which(RIS_raw$POPULATION %in% popns.to.use$POPULATION), ]

summary(RIS_allyears)

## now extract abundance and phenology trends for each of these populations, as before

RIS_allyears$MEAN.ABUND <- RIS_allyears$COUNT/RIS_allyears$RECS


# let's generate some summary stats about how many RIS populations and how many population*year records go into the final dataset
# first, the overall stat of how many sites in total are included (should be 110 as before)
length(levels(droplevels(RIS_allyears$SITE)))

# and now some species-level stats

summary.stats.ay <- data.frame(SCI_NAME = factor(),
                               COMMON_NAME = factor(),
                               SITES = numeric(),
                               SITE.YEARS = numeric(),
                               RECORDS = numeric(),
                               INDIVIDUALS = numeric(),
                               FIRST.YEAR = numeric(),
                               LAST.YEAR = numeric(),
                               DURATION = numeric())

for (x in RIS.moths){
  spec <- RIS_allyears[which(RIS_allyears$COMMON_NAME==x), ]
  SCI_NAME <- as.character(spec[[1,3]])
  spec$YEARS <- 1
  by.year <- ddply(spec, .(SCI_NAME,COMMON_NAME,SITE), summarise,
                   SITE.YEARS = sum(YEARS),
                   RECORDS = sum(OBS),
                   INDIVIDUALS = sum(COUNT),
                   FIRST.YEAR = min(YEAR),
                   LAST.YEAR = max(YEAR))
  SITES <- nrow(by.year)
  SITE.YEARS <- sum(by.year$SITE.YEARS)
  RECORDS <- sum(by.year$RECORDS)
  INDIVIDUALS <- sum(by.year$INDIVIDUALS)
  FIRST.YEAR <- min(by.year$FIRST.YEAR)
  LAST.YEAR <- max(by.year$LAST.YEAR)
  DURATION <- LAST.YEAR - FIRST.YEAR
  
  out <- data.frame(cbind(SCI_NAME,x,SITES,SITE.YEARS,RECORDS,INDIVIDUALS,FIRST.YEAR,LAST.YEAR,DURATION))
  summary.stats.ay <- rbind(summary.stats.ay,out)
}

colnames(summary.stats.ay) <- c("SCI_NAME","COMMON_NAME","SITES","SITE.YEARS","RECORDS","INDIVIDUALS","FIRST.YEAR","LAST.YEAR","DURATION")

summary.stats.ay$SITES <- as.numeric(as.character(summary.stats.ay$SITES))
summary.stats.ay$SITE.YEARS <- as.numeric(as.character(summary.stats.ay$SITE.YEARS))
summary.stats.ay$RECORDS <- as.numeric(as.character(summary.stats.ay$RECORDS))
summary.stats.ay$INDIVIDUALS <- as.numeric(as.character(summary.stats.ay$INDIVIDUALS))
summary.stats.ay$FIRST.YEAR <- as.numeric(as.character(summary.stats.ay$FIRST.YEAR))
summary.stats.ay$LAST.YEAR <- as.numeric(as.character(summary.stats.ay$LAST.YEAR))
summary.stats.ay$DURATION <- as.numeric(as.character(summary.stats.ay$DURATION))

summary(summary.stats.ay)


# now we can calculate the trends through the emergence dates over time
# this is a slight tweak on what we've done above - rather than getting a single value in each year and putting a simple trend through it,
# let's fit a model to the values for all sites, with site as a random effect
# for this, we want to invert the slope so that this is a measure of phenological *advancement* as opposed to phenological *change*

pheno.slopes.ay <- data.frame(SCI_NAME = factor(),
                              COMMON_NAME = factor(),
                              PHENO.SLOPE = numeric(),
                              PHENO.SE = numeric(),
                              PHENO.CHI = numeric(),
                              PHENO.P = numeric())


for (x in RIS.moths){
  spec <- RIS_allyears[which(RIS_allyears$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Phenology/FullRange/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(spec$PEAKDAY ~ spec$YEAR,
       ylim = c(0,365),
       xlab = "YEAR", ylab = "Annual julian day of peak first-generation emergence")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (PHENO.P == 0){
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",PHENO.SLOPE," +/- ",PHENO.SE,"; chisq = ",PHENO.CHI,"; p = ",PHENO.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, PHENO.SLOPE, PHENO.SE, PHENO.CHI, PHENO.P))
  pheno.slopes.ay <- rbind(pheno.slopes.ay, out)
  
}

pheno.slopes.ay$PHENO.SLOPE <- as.numeric(as.character(pheno.slopes.ay$PHENO.SLOPE))
pheno.slopes.ay$PHENO.SE <- as.numeric(as.character(pheno.slopes.ay$PHENO.SE))
pheno.slopes.ay$PHENO.CHI <- as.numeric(as.character(pheno.slopes.ay$PHENO.CHI))
pheno.slopes.ay$PHENO.P <- as.numeric(as.character(pheno.slopes.ay$PHENO.P))

summary(pheno.slopes.ay)



### now we want to replicate this process for abundance data from the same set of RIS populations
# (we will have good abundance )
# this data is already read in - 
summary(RIS_allyears$MEAN.ABUND)

abund.slopes.ay <- data.frame(SCI_NAME = factor(),
                              COMMON_NAME = factor(),
                              ABUND.EXP.SLOPE = numeric(),
                              ABUND.SLOPE = numeric(),
                              ABUND.SE = numeric(),
                              ABUND.CHI = numeric(),
                              ABUND.P = numeric())


for (x in RIS.moths){
  spec <- RIS_allyears[which(RIS_allyears$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = spec)
  
  ABUND.EXP.SLOPE <- round(exp(summary(test1)$coefficients[2,1]), 4)
  ABUND.SLOPE <- round(summary(test1)$coefficients[2,1], 4)
  ABUND.SE <- round(summary(test1)$coefficients[2,2], 4)
  ABUND.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  ABUND.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Abundance/FullRange/Moths/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$MEAN.ABUND) ~ spec$YEAR,
       ylim = c(0,10),
       xlab = "YEAR", ylab = "Log (mean no. recorded individuals per transect)")
  
  abline(coef[1],coef[2])
  title(main = paste(x,SCI_NAME, sep = " - "))
  if (ABUND.P == 0){
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"); chisq = ",ABUND.CHI,"; p < 0.0001"), line = 0)
  } else {
    mtext(paste0("effect: ",ABUND.EXP.SLOPE,"(",ABUND.SLOPE," +/- ",ABUND.SE,"; chisq = ",ABUND.CHI,"; p = ",ABUND.P), line = 0)
  }
  
  dev.off()
  
  out <- data.frame(cbind(SCI_NAME, COMMON_NAME, ABUND.EXP.SLOPE, ABUND.SLOPE, ABUND.SE, ABUND.CHI, ABUND.P))
  abund.slopes.ay <- rbind(abund.slopes.ay, out)
  
}

abund.slopes.ay$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes.ay$ABUND.EXP.SLOPE))
abund.slopes.ay$ABUND.SLOPE <- as.numeric(as.character(abund.slopes.ay$ABUND.SLOPE))
abund.slopes.ay$ABUND.SE <- as.numeric(as.character(abund.slopes.ay$ABUND.SE))
abund.slopes.ay$ABUND.CHI <- as.numeric(as.character(abund.slopes.ay$ABUND.CHI))
abund.slopes.ay$ABUND.P <- as.numeric(as.character(abund.slopes.ay$ABUND.P))

summary(abund.slopes.ay)


# merge all this together

RIS_out_pre.ay <- merge(pheno.slopes.ay,abund.slopes.ay)
RIS_out.ay <- merge(summary.stats.ay,RIS_out_pre.ay)


# and write out all versions
write.csv(RIS_out.ay, "Data/Derived/RIS_ay.csv", row.names = F)


