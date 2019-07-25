####################################################
####   Script for preparing data for analysis   ####
####################################################

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","rnrfa","rgdal","blighty","dplyr","raster","RColorBrewer","lme4","ggplot2","lubridate","mgcv","gridExtra","reshape2","effects")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# we have data from four sources - BNM, NMRS, UKBMS and RIS
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
# (BNM/NMRS = hectad-level annual presence/absence; UKBMS/RIS = site-level daily count data)

# we need to apply these to select out species to use from the hectad data, and species/site combinations to use from the site data

### BNM data ####
# first let's do it on the hectad data for butterflies

## read in the data
BNM_raw <- read.csv("../../Data/BNM/Butterflies/Report - BNM Butterflies 10km 1950-2014.csv", header = TRUE)

summary(BNM_raw)


# tweak the colnames to a standard

colnames(BNM_raw) <- c("SCI_NAME","COMMON_NAME","HECTAD","YEAR")

BNM_raw$fYEAR <- as.factor(BNM_raw$YEAR)

# attach a column to indicate presence

BNM_raw$PRESENCE <- 1

## next, we want to attach eastings and northings to this data so that it can be plotted
# as there are > 1000 repeats of each hectad let's generate a non-redundant list

BNM_hectads <- ddply(BNM_raw, .(HECTAD), summarise, COUNT = sum(PRESENCE))

# at the moment the dataset includes a whole bunch of records from sites in Northern Ireland, and some in the Channel Islands.
# that's nice but I don't want them - they will cause problems elsewhere (esp. when analysing range-margin stuff)

# the question is, how to pick them out and get rid?
# first let's take first two characters of every grid ref: these are the 100km square

BNM_hectads$SQUARE <- as.factor(substr(BNM_hectads$HECTAD, 1,2))

# now, all sites in GB have two letters, sites in NI have one letter and one number
# we can use this to separate them

BNM_hectadsGB <- BNM_hectads[which(!grepl("[[:digit:]]", BNM_hectads$SQUARE)),]

# to get rid of Channel Islands, we need to get rid of squares WA and WV

BNM_hectadsGB <- BNM_hectadsGB[which(!(BNM_hectadsGB$SQUARE %in% c("WA","WV"))), ]


# now generate a vector containing these site numbers only

# first we need to drop unused levels (i.e. those from NI)
BNM_hectadsGB <- droplevels(BNM_hectadsGB)


# now we are ready to parse the easting/northings

latlon <- osg_parse(BNM_hectadsGB$HECTAD)

BNM_hectadsGB$EASTING <- latlon[[1]]
BNM_hectadsGB$NORTHING <- latlon[[2]]

### we now have eastings and northings! let's do a brief visual check that all's well by plotting them on a map...

# to do this we actually also need these with fewer digits - specifically, in units of 1 km
BNM_hectadsGB$east.km <- BNM_hectadsGB$EASTING/1000
BNM_hectadsGB$north.km <- BNM_hectadsGB$NORTHING/1000

# plot base map to test
blighty()

# plot all records to test
# we need to offset all points slightly for them to be centred in the middle of hectads, rather than the bottom-left corner
points(5+BNM_hectadsGB$east.km,5+BNM_hectadsGB$north.km,pch=16) # this looks a bit messy but looks fine in a zoom window


# great - but this is a total mess!
# we also need to merge these eastings and northings back into the main data
BNM_hectadsGB <- BNM_hectadsGB[,-2]


# we can now use this to trim down the other data to only include these sites
BNM_GB <- merge(BNM_raw,BNM_hectadsGB)



# let's plot an example species to get a distribution

BNM_GB.db <- BNM_GB[which(BNM_GB$COMMON_NAME=="Duke of Burgundy"), ]
blighty()
points((5+BNM_GB.db$east.km),(5+BNM_GB.db$north.km),pch=16)
title(main = "Duke of Burgundy")


# let's also try looping this across the range of years of interest (1986-1995 and 2001-2010) to see whether the expansion is visible

# set up the period of interest
years <- 1995:2014

for (x in years){
  year <- BNM_GB.db[which(BNM_GB.db$YEAR==x), ]
  png(paste0("Data/BNM/DukeOfBurgundy/Range/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points((5+year$east.km),(5+year$north.km),pch=16)
  title(main = x)
  dev.off()
}


## now restrict the data down to fit data-quality filters


### I want to select only 'well-recorded' hectads to take forwards
# this follows Hickling et al 2005, 'common squares' (see Hassall & Thompson 2010)
# the idea of this is that we only use hectads where we are more confident that we're looking at a true absence rather than a lack of recording
# i.e. we check every hectad individually and only use those where a high enough proportion of the regional fauna
# has been recorded at least once in each time window


## first, let's create a non-redundant list of recording in each hectad in each year, dropping the species information

BNM_rec.hectad.years <- ddply(BNM_GB, .(HECTAD, YEAR), numcolwise(mean))

# for inspection, let's plot the total hectad-level recording across all years

years.extended <- 1980:2014

for (n in years.extended){
  year <- BNM_rec.hectad.years[which(BNM_rec.hectad.years$YEAR==n), ]
  png(paste0("Data/BNM/HectadRecording/",n,".png"),width = 612, height = 627, units = "px", bg = "white")
  blighty()
  points(year$east.km,year$north.km,pch=16)
  title(main = n)
  dev.off()
}

# we can see that recording was pretty sparse at some points and in some places, especially in the North
# but gets fairly robust especially in most recent years

# let's actually count the number of recorded hectads and the number of species*hectad records in each year 

recording <- data.frame(YEAR = numeric(),
                        HECTADS = numeric(),
                        RECORDS = numeric())

for (n in years.extended){
  year.hec <- BNM_rec.hectad.years[which(BNM_rec.hectad.years$YEAR==n), ]
  HECTADS <- nrow(year.hec)
  
  year.rec <- BNM_GB[which(BNM_GB$YEAR==n), ]
  RECORDS <- nrow(year.rec)
  
  YEAR <- n
  
  out <- data.frame(cbind(YEAR,HECTADS,RECORDS))
  recording <- rbind(recording,out)
  
  
}

summary(recording)


# now plot the number of hectads and records over time
plot(HECTADS ~ YEAR, data = recording)
plot(RECORDS ~ YEAR, data = recording)

total.hecs <- max(recording$HECTADS)

recording$PERC <- recording$HECTADS*100/total.hecs

plot(PERC ~ YEAR, data = recording)



# for butterflies, there is a really obvious jump in recording in 1995, and thereafter it's at a pretty stable level

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
BNM_rec.hectad.years.inter1 <- BNM_rec.hectad.years[which(BNM_rec.hectad.years$YEAR %in% inter1), ]
BNM_rec.hectad.years.inter2 <- BNM_rec.hectad.years[which(BNM_rec.hectad.years$YEAR %in% inter2), ]

# then get rid of the year information to leave a non-redundant list of hectads recorded in each interval

BNM_rec.hectads.inter1 <- ddply(BNM_rec.hectad.years.inter1, .(HECTAD), numcolwise(mean))
BNM_rec.hectads.inter2 <- ddply(BNM_rec.hectad.years.inter2, .(HECTAD), numcolwise(mean))


# plot these on a map

blighty()
points(BNM_rec.hectads.inter2$east.km, BNM_rec.hectads.inter2$north.km, pch=16, col="red")
points(BNM_rec.hectads.inter1$east.km, BNM_rec.hectads.inter1$north.km, pch=16)


# we can see that almost every hectad was recorded in interval 2, and only about 15 fewer in interval 1

## now we want to recombine these datasets in a way that ditches anything not present in both

BNM_rec.hectads.inter1.vec <- as.character(BNM_rec.hectads.inter1$HECTAD)
BNM_rec.hectads.inter2.vec <- as.character(BNM_rec.hectads.inter2$HECTAD)

BNM_rec.hectads.full <- append(BNM_rec.hectads.inter1.vec,BNM_rec.hectads.inter2.vec)

BNM_rec.hectads.full <- data.frame(cbind(BNM_rec.hectads.full,1))
colnames(BNM_rec.hectads.full) <- c("HECTAD","INTERVALS")
BNM_rec.hectads.full$HECTAD <- as.factor(as.character(BNM_rec.hectads.full$HECTAD))
BNM_rec.hectads.full$INTERVALS <- as.numeric(as.character(BNM_rec.hectads.full$INTERVALS))


BNM_rec.hectads.good <- ddply(BNM_rec.hectads.full, .(HECTAD), numcolwise(sum))

BNM_rec.hectads.good <- BNM_rec.hectads.good[which(BNM_rec.hectads.good$INTERVALS > 1), ]

# this leaves us with a list of 2748 'good hectads'


## now we can use this list of 'good' hectads to select out data from the main set to use
# first create a vector of 'good' hectads

# drop unused levels
BNM_rec.hectads.good <- droplevels(BNM_rec.hectads.good)

# generate the vector
BNM_good.hectads <- levels(BNM_rec.hectads.good$HECTAD)

# pick out the good records
BNM_GB_good <- BNM_GB[which(BNM_GB$HECTAD %in% BNM_good.hectads), ]

# it's promising that we lose only < 0.1% of records

## we can also further restrict to the hectads that have records of a minimum percentage of the regional species richness in each year
# this is going to take quite a complex loop to construct:
# for every hectad in the BNM_GB_good dataframe, we need to calculate the observed species richness
# then we need to calculate pairwise distances to all other recorded hectads, pick out the nearest 100, 
# and calculate the regional species richness of these 100
# and finally turn this into a percentage of regional SR that is recorded in the hectad

# start with the list of recorded hectads, re-attach the eastings and northings
BNM_good <- ddply(BNM_rec.hectad.years, .(HECTAD,EASTING,NORTHING), summarise,
                  COUNT = sum(PRESENCE))

# select out only the good hectads
BNM_good <- BNM_good[which(BNM_good$HECTAD %in% BNM_good.hectads), ]

# now start to construct the loop:
hec.records <- data.frame(YEAR = numeric(),
                          HECTAD = factor(),
                         EASTING = numeric(),
                         NORTHING = numeric(),
                         hec.SR = numeric(),
                         reg.SR = numeric(),
                         hec.perc = numeric())


for (x in BNM_good.hectads){
  hec <- BNM_good[which(BNM_good$HECTAD == x), ]    # pick out details of focal hectad
  candidates <- BNM_good[which(BNM_good$HECTAD != x), ]   # pick out details of all others
  
  hec_east <- hec[[1,2]]     # easting of target
  hec_north <- hec[[1,3]]    # northing of target
  
  candidates$east_diff <- candidates$EASTING - hec_east             # longitudinal difference
  candidates$north_diff <- candidates$NORTHING - hec_north          # latitudinal difference
  
  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
  
  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
  
  closest <- candidates[1:100,] # select out closest 100

  # now we want to calculate species richness from the hectad and from the region in each year
  for (n in years){
    year.recs <- BNM_GB_good[which(BNM_GB_good$YEAR == n), ]
    
  hec.recs <- year.recs[which(year.recs$HECTAD == x), ] # pull out hectad records
  hec.SR <- nlevels(droplevels(hec.recs$COMMON_NAME))       # calculate species richness
  
  reg.recs <- year.recs[which(year.recs$HECTAD %in% droplevels(closest$HECTAD)), ]  # pull out region records
  reg.recs <- rbind(reg.recs,hec.recs)    # add in hectad records (they're part of the regional richness too!)
  reg.SR <- nlevels(droplevels(reg.recs$COMMON_NAME))
  
  hec.perc <- hec.SR*100/reg.SR

  out <- cbind(x,n,hec_east,hec_north,hec.SR,reg.SR,hec.perc)
  hec.records <- rbind(hec.records,out)
  }    
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

# it's pretty promising (but not necessarily surprising for butterflies) 
# that the majority (> 1/2) of the hectad*year combinations are "heavily recorded" (i.e. > 25% of regional species richness) 
# this may prove not to be the case for moths, later...

# this data takes a while to generate so let's back it up
write.table(hec.records, "Data/BNM/HectadRecording/RecordingLevels.txt", row.names = F)

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
legend("topright", col=c("Grey90","Grey70","Red","Black"),pch=16,
       legend=c("Not recorded","Recorded","Well recorded","Heavily recorded"))

# almost all of the squares that are "well" rather than "heavily" recorded are in Highlands & Islands
# especially Shetland which is unsurprising, as it's species-poor but its nearest 100 hectads includes a large chunk of the mainland

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
png(paste0("Data/BNM/HectadRecording/RecordingLevels.png"), width = 1600, height = 800, units = "px", bg = "white")
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
# and also towards some northern rare-butterfly hotspots (e.g. Morecambe Bay, North York Moors, Oban)



# we can also plot out the hectads so we can visualise where their region lies...
# but this takes a LONG time so once I've done it I'll blank it out!
#for (x in BNM_good.hectads){
#  hec <- BNM_good[which(BNM_good$HECTAD == x), ]    # pick out details of focal hectad
#  candidates <- BNM_good[which(BNM_good$HECTAD != x), ]   # pick out details of all others
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
#  png(paste0("Data/BNM/HectadRecording/HectadRegions/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
#  
#  blighty()
#  points(5+(candidates$EASTING/1000), 5+(candidates$NORTHING/1000),pch=16)
#  points(5+(closest$EASTING/1000), 5+(closest$NORTHING/1000),pch=16,col="red")
#  points(5+(hec_east/1000),5+(hec_north/1000),pch=16,col="royalblue")
#  title(main = x)
#  
#  dev.off()
#}

# so the main dataframe (BNM_GB_good) is effectively our dataframe for 'recorded'
# let's generate an additional one for each of 'well recorded' and 'heavily recorded'
hec.records.overall$HECTAD <- as.factor(hec.records.overall$HECTAD)

rec.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC > 0), ]
well.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC >= 10), ]
heavy.hecs <- hec.records.overall[which(hec.records.overall$MED.PERC >= 25), ]

BNM_GB_rec <- BNM_GB_good[which(BNM_GB_good$HECTAD %in% droplevels(rec.hecs$HECTAD)), ]
BNM_GB_well <- BNM_GB_good[which(BNM_GB_good$HECTAD %in% droplevels(well.hecs$HECTAD)), ]
BNM_GB_heavy <- BNM_GB_good[which(BNM_GB_good$HECTAD %in% droplevels(heavy.hecs$HECTAD)), ]


### selecting species ####
## so, we have our good hectads: now we need to select which species to use from those that remain in the data

## we'll need to repeat this for all three recording levels. Let's first do recorded.

### recorded hectads ####
# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
BNM_GB_good.inters <- BNM_GB_rec[which(BNM_GB_rec$YEAR %in% years), ]

# get rid of the year information, keeping only which hectads each species was recorded in:
BNM_GB_good.hecs <- ddply(BNM_GB_good.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                          PRESENCE = 1)

# and now count the number of hectads per species
species.hecs <- ddply(BNM_GB_good.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                      HECTADS = sum(PRESENCE))

print(species.hecs)

# here, everything gets retained except for Large Chequered Skipper (which bizarrely has a single record near Derby in 2006?),
# Lulworth Skipper, which is only recorded from 15 hectads,
# and Heath Fritillary, which is only recorded from 19
# let's add these to a list of species to remove; whilst we're here, we'll also manually stick migrants on this list

under.recorded.species <- species.hecs[which(species.hecs$HECTADS < 20), ]

print(under.recorded.species)

under.recorded <- levels(droplevels(under.recorded.species$COMMON_NAME))

migrants <- c("Red Admiral","Clouded Yellow","Painted Lady","Swallowtail","Large Tortoiseshell")
spec.to.cut <- append(under.recorded,migrants)

# n.b. Swallowtail P. machaon ssp. britannicus is not a migrant, but this data does not distinguish between records of this
# and records of P. machaon ssp. gorganus, which is, leading to some weird results - so it's safer just to exclude it

# next we want to focus on southerly species: i.e. those that reach a northern range margin in the UK, are likely climate-limited here,
# and have the capacity to expand northwards

# the rule for this is that the species' northern range margin across the entire period can't be within 100 km of John O'Groats
# so obviously for this, we need to estimate the northern range margins!

# these are the mean northing of the 10 most northerly occupied hectads in the interval

# calculate the margins for each species
margins <- data.frame(SCI_NAME = factor(),
                      NORTHING = numeric(),
                      COMMON_NAME = factor())

species.list <- levels(droplevels(species.hecs$COMMON_NAME))
species.list <- species.list[which(!(species.list %in% spec.to.cut))]

for (x in species.list){
  sp.dat <- BNM_GB_good.hecs[which(BNM_GB_good.hecs$COMMON_NAME==x), ]
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

print(margins)

# here we lose quite a range of species with distributions that extend into Scotland
# some are truly ubiquitous -
# e.g. Common Blue, Green-veined White and Meadow Brown
# some have patchy distributions but happen to be present in northern Scotland -
# e.g. Grayling, Small Blue
# and some are northerly-distributed
# e.g. Large Heath, Scotch Argus

# none of these are climate-limited (southerly-distributed) so we're fine to get rid of all of them 
# however we *will* want the data on them, so don't append them to the spec.to.cut list

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
  spec <- BNM_GB[which(BNM_GB$COMMON_NAME == x),]
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

# finally, we only want species that have been recorded in a minimum proportion of the 20 years
# count how many years each species has a record (in any hectad) for 

BNM_GB_good.years <- ddply(BNM_GB_good.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

BNM_GB_good.years.count <- ddply(BNM_GB_good.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(BNM_GB_good.years.count)

# now retain only species recorded in at least 15/20 years
insufficient.years.spec <- BNM_GB_good.years.count[which(BNM_GB_good.years.count$YEARS < 20), ]

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
# and the abundance change and phenology advance from the UKBMS data

# since we have the data read in already, let's start with the distribution change
# we need to calculate the percentage of recorded hectads (from those chosen for retention based on recording level) that were occupied in each year
# and plot a slope through time to establish whether this has (significantly) changed

## select out data from years of interest
hectads.years <- BNM_GB_good[which(BNM_GB_good$YEAR %in% years), ]

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
  
  png(paste0("Data/Derived/Distribution/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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
  distrib.slopes <- rbind(distrib.slopes, out)
  
}

distrib.slopes$DISTRIB.SLOPE <- as.numeric(as.character(distrib.slopes$DISTRIB.SLOPE))
distrib.slopes$DISTRIB.SE <- as.numeric(as.character(distrib.slopes$DISTRIB.SE))
distrib.slopes$DISTRIB.CHI <- as.numeric(as.character(distrib.slopes$DISTRIB.CHI))
distrib.slopes$DISTRIB.P <- as.numeric(as.character(distrib.slopes$DISTRIB.P))

summary(distrib.slopes)


# this puts the biggest losers as being High Brown Fritillary, Wood White and Wall (which could match our expectations)
# and the biggest winners as Purple Emperor - perhaps surprising - 
# but then Silver-washed Fritillary, Ringlet, Silver-spotted Skipper, Speckled Wood and Adonis Blue (which matches our expectations)



## and now northern range margins
# we want to select out the ten northernmost hectads for each species in each year
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
  png(paste0("Data/Derived/RangeMargins/Margins/Butterflies/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
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

# this looks pretty decent!
# strong advances in the species we'd expect - e.g. Comma, Brown Argus
# and evidence of retractions in some examples that might also be expected - e.g. Grizzled Skipper (disappearing from Lincs sites)



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



# we have a couple of retractions e.g. Grizzled Skipper where northerly populations have gone extinct
# and some big expansions of which Comma is unsurprisingly the biggest

# merge these variables into a single data frame to output

derived.BNM <- merge(margin.slopes,distrib.slopes)

# write it out
write.csv(derived.BNM, file = "Data/BNM/Derived/derived.csv", row.names = F)



### well recorded hectads ####
# we now want to repeat all of the above using only data from the well-recorded, and later the heavily-recorded, hectads
# although this only trims out a couple of hectads it might make a small difference somewhere
# and might be more important for the moths (so important for reproducibility)

# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
BNM_GB_well.inters <- BNM_GB_well[which(BNM_GB_well$YEAR %in% years), ]

# get rid of the year information, keeping only which hectads each species was recorded in:
BNM_GB_well.hecs <- ddply(BNM_GB_well.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                          PRESENCE = 1)

# and now count the number of hectads per species
species.hecs.well <- ddply(BNM_GB_well.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                      HECTADS = sum(PRESENCE))

print(species.hecs.well)

# here, everything gets retained except for Large Chequered Skipper (which bizarrely has a single record near Derby in 2006?),
# Lulworth Skipper, which is only recorded from 15 hectads,
# and Heath Fritillary, which is only recorded from 19
# let's add these to a list of species to remove; whilst we're here, we'll also manually stick migrants on this list

under.recorded.species.well <- species.hecs.well[which(species.hecs.well$HECTADS < 20), ]

print(under.recorded.species.well)

under.recorded.well <- levels(droplevels(under.recorded.species.well$COMMON_NAME))
spec.to.cut.well <- append(under.recorded.well,migrants)

# n.b. Swallowtail P. machaon ssp. britannicus is not a migrant, but this data does not distinguish between records of this
# and records of P. machaon ssp. gorganus, which is, leading to some weird results - so it's safer just to exclude it

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
  sp.dat <- BNM_GB_well.hecs[which(BNM_GB_well.hecs$COMMON_NAME==x), ]
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

BNM_GB_well.years <- ddply(BNM_GB_well.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

BNM_GB_well.years.count <- ddply(BNM_GB_well.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(BNM_GB_well.years.count)

# now retain only species recorded in at least 15/20 years
insufficient.years.spec.well <- BNM_GB_well.years.count[which(BNM_GB_well.years.count$YEARS < 20), ]

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
hectads.years.well <- BNM_GB_well[which(BNM_GB_well$YEAR %in% years), ]

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
  
  png(paste0("Data/Derived/Distribution/Butterflies/WellRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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
  png(paste0("Data/Derived/RangeMargins/Margins/Butterflies/WellRecorded/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
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
  
  png(paste0("Data/Derived/RangeMargins/Change/Butterflies/WellRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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

derived.BNM.well <- merge(margin.slopes.well,distrib.slopes.well)

# write it out
write.csv(derived.BNM.well, file = "Data/BNM/Derived/derived_wellrecorded.csv", row.names = F)




### heavily recorded hectads ####

# the first step is to remove species which were recorded in fewer than 20 hectads across both time intervals

# select out the data from years of interest
BNM_GB_heavy.inters <- BNM_GB_heavy[which(BNM_GB_heavy$YEAR %in% years), ]

# get rid of the year information, keeping only which hectads each species was recorded in:
BNM_GB_heavy.hecs <- ddply(BNM_GB_heavy.inters, .(HECTAD,SCI_NAME,COMMON_NAME,SQUARE,EASTING,NORTHING), summarise,
                          PRESENCE = 1)

# and now count the number of hectads per species
species.hecs.heavy <- ddply(BNM_GB_heavy.hecs, .(SCI_NAME,COMMON_NAME), summarise,
                           HECTADS = sum(PRESENCE))

print(species.hecs.heavy)

# here, everything gets retained except for Large Chequered Skipper (which bizarrely has a single record near Derby in 2006?),
# Lulworth Skipper, which is only recorded from 15 hectads,
# and Heath Fritillary, which is only recorded from 19
# let's add these to a list of species to remove; whilst we're here, we'll also manually stick migrants on this list

under.recorded.species.heavy <- species.hecs.heavy[which(species.hecs.heavy$HECTADS < 20), ]

print(under.recorded.species.heavy)

under.recorded.heavy <- levels(droplevels(under.recorded.species.heavy$COMMON_NAME))
spec.to.cut.heavy <- append(under.recorded.heavy,migrants)

# n.b. Swallowtail P. machaon ssp. britannicus is not a migrant, but this data does not distinguish between records of this
# and records of P. machaon ssp. gorganus, which is, leading to some weird results - so it's safer just to exclude it

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
  sp.dat <- BNM_GB_heavy.hecs[which(BNM_GB_heavy.hecs$COMMON_NAME==x), ]
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
# some are truly ubiquitous -
# e.g. Common Blue, Green-veined White and Meadow Brown
# some have patchy distributions but happen to be present in northern Scotland -
# e.g. Grayling, Small Blue
# and some are northerly-distributed
# e.g. Large Heath, Scotch Argus

# none of these are climate-limited (southerly-distributed) so we're fine to get rid of all of them 
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

BNM_GB_heavy.years <- ddply(BNM_GB_heavy.inters, .(SCI_NAME,COMMON_NAME,YEAR), summarise, PRESENCE = 1)

BNM_GB_heavy.years.count <- ddply(BNM_GB_heavy.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = sum(PRESENCE))

summary(BNM_GB_heavy.years.count)

# now retain only species recorded in at least 15/20 years
insufficient.years.spec.heavy <- BNM_GB_heavy.years.count[which(BNM_GB_heavy.years.count$YEARS < 20), ]

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
hectads.years.heavy <- BNM_GB_heavy[which(BNM_GB_heavy$YEAR %in% years), ]

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

## we also want a scaled percentage (i.e. PERC compared to PERC in 1995)
# this gives us an estimate of distribution size change relative to the size of the starting distribution

# extract out the 1995 estimates
hectads.1995.heavy <- hectads.annualtotals.heavy[which(hectads.annualtotals.heavy$YEAR == 1995), c(2:3,6)]
colnames(hectads.1995.heavy)[3] <- "START.PERC.REC"

# merge them back in 
hectads.annualtotals.heavy <- merge(hectads.annualtotals.heavy,hectads.1995.heavy)

# and calculate the scaled distribution change
hectads.annualtotals.heavy$PERC.CHANGE <- hectads.annualtotals.heavy$PERC.REC*100/hectads.annualtotals.heavy$START.PERC.REC


### at this point we want to quickly calculate the total number of records (i.e. hectads x years x species)
sum(hectads.annualtotals.heavy$HECTADS)


## and generate species-level summaries
BNM.summary.heavy <- data.frame(SCI_NAME = factor(),
                                COMMON_NAME = factor(),
                                HECTAD.RECORDS = numeric(),    # i.e. total number of records summed across all years
                                RECORDED.HECTADS = numeric())  # i.e. total number of hectads recorded in any year

for (x in spec.to.use.heavy){
  spec.anntot <- hectads.annualtotals.heavy[which(hectads.annualtotals.heavy$COMMON_NAME==x), ]
  
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec.anntot[[1,1]])
  
  HECTAD.RECORDS <- sum(spec.anntot$HECTADS)
  
  spec <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME==x), ]
  spec.coll <- ddply(spec, .(HECTAD,SCI_NAME,COMMON_NAME), summarise,
                     YEARS = sum(PRESENCE))
  
  RECORDED.HECTADS <- nrow(spec.coll)

  out <- data.frame(cbind(SCI_NAME,COMMON_NAME,HECTAD.RECORDS,RECORDED.HECTADS))
  BNM.summary.heavy <- rbind(BNM.summary.heavy,out)
}

summary(BNM.summary.heavy)


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
  
  png(paste0("Data/Derived/Distribution/Butterflies/HeavyRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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
  png(paste0("Data/Derived/RangeMargins/Margins/Butterflies/HeavyRecorded/",x,".png"),width = 612, height = 627, units = "px", bg = "white")
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
  
  png(paste0("Data/Derived/RangeMargins/Change/Butterflies/HeavyRecorded/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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



# we have a couple of retractions e.g. Grizzled Skipper where northerly populations have gone extinct
# and some big expansions of which Comma is unsurprisingly the biggest

# merge these variables into a single data frame to output

derived.BNM.heavy.pre <- merge(BNM.summary.heavy,margin.slopes.heavy)
derived.BNM.heavy <- merge(derived.BNM.heavy.pre,distrib.slopes.heavy)

# write it out
write.csv(derived.BNM.heavy, file = "Data/BNM/Derived/derived_heavyrecorded.csv", row.names = F)






### UKBMS data ####
# now let's pull in the UKBMS data and see how much of that passes its quality filters - and what species remain when we combine both datasets

# read in the data
BMS_raw <- read.csv("../../Data/UKBMS/GAMphenology_refined.csv", header = T)

summary(BMS_raw)

# we have already applied a couple of the QC filters in generating this file -
# we only have rows of data here where a transect was recorded 10 times in a year,
# and the species of interest was recorded on 3 of those
# we also have "Fail" entered under PEAKDAY wherever the GAM either doesn't start the year at zero, 
# or fails to peak and decrease before the end of the year
# and we have quite a few instances of "Fail" under the voltinism ratio, which we don't want to chuck

# we can dispose of the 'fails' now, and make the response variable numeric again
BMS_raw <- BMS_raw[which(BMS_raw$PEAKDAY != "Fail"), ]
BMS_raw$PEAKDAY <- as.numeric(as.character(BMS_raw$PEAKDAY)) 

summary(BMS_raw)

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
BMS_prequel <- BMS_raw[which(BMS_raw$YEAR %in% prequel), ]

# collapse down to a list of populations, counting how many years each was recorded in
BMS_prequel_popns <- ddply(BMS_prequel, .(siteXspecies,SCI_NAME,SPECIES,COMMON_NAME,SITE), summarise, YEARS.PREQUEL = length(PEAKDAY))


## interval: present in all years 1995-2014

# select out the data
BMS_years <- BMS_raw[which(BMS_raw$YEAR %in% years), ]

# collapse down, counting years
BMS_years_popns <- ddply(BMS_years, .(siteXspecies,SCI_NAME,SPECIES,COMMON_NAME,SITE), summarise, YEARS = length(PEAKDAY))


# try a few different threshold levels -
# recorded in *every* year
BMS_years_t20 <- BMS_years_popns[which(BMS_years_popns$YEARS>=20), ]

# recorded in at least 18/20 years
BMS_years_t18 <- BMS_years_popns[which(BMS_years_popns$YEARS>=18), ]

# recorded in at least 15/20 years
BMS_years_t15 <- BMS_years_popns[which(BMS_years_popns$YEARS>=15), ]

# recorded in at least 10/20 years
BMS_years_t10 <- BMS_years_popns[which(BMS_years_popns$YEARS>=10), ]




## try merging these together to produce candidate datasets at each threshold

BMS_t20 <- merge(BMS_years_t20, BMS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
BMS_t20_species <- ddply(BMS_t20, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

BMS_t20_good <- BMS_t20_species[which(BMS_t20_species$POPULATIONS >= 3), ]

BMS_t20_use <- BMS_t20_good[which(BMS_t20_good$COMMON_NAME %in% spec.to.use), ]

# setting the threshold at 20 (i.e. must have been recorded in *every* year) leaves very few populations behind,
# and most of those are populations of ubiquitous species that won't be included in the NRM analysis
# this leaves us with only 11 species using this approach

# ideally we want to retain as many of the good BNM species as possible
# therefore we'll need to use a less strict threshold - let's try some of the other options

# t18
BMS_t18 <- merge(BMS_years_t18,BMS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
BMS_t18_species <- ddply(BMS_t18, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

BMS_t18_good <- BMS_t18_species[which(BMS_t18_species$POPULATIONS >= 3), ]

BMS_t18_use <- BMS_t18_good[which(BMS_t18_good$COMMON_NAME %in% spec.to.use), ]

# this immediately doubles the number of species to 20!
# let's delve even further...

# t15
BMS_t15 <- merge(BMS_years_t15,BMS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
BMS_t15_species <- ddply(BMS_t15, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

BMS_t15_good <- BMS_t15_species[which(BMS_t15_species$POPULATIONS >= 3), ]

BMS_t15_use <- BMS_t15_good[which(BMS_t15_good$COMMON_NAME %in% spec.to.use), ]

# now up to 29 species...

# t10
BMS_t10 <- merge(BMS_years_t10,BMS_prequel_popns)

# collapse down further to count populations per species - we need 3 'good' populations to use a species
BMS_t10_species <- ddply(BMS_t10, .(SCI_NAME,SPECIES,COMMON_NAME), summarise, POPULATIONS = length(YEARS.PREQUEL))

BMS_t10_good <- BMS_t10_species[which(BMS_t10_species$POPULATIONS >= 3), ]

BMS_t10_use <- BMS_t10_good[which(BMS_t10_good$COMMON_NAME %in% spec.to.use), ]

# even dropping it this far adds even more species to the dataset - up to 37



# summarise which species to use/exclude at all possible levels

butterflies <- data.frame(COMMON_NAME = levels(droplevels(BNM_raw$COMMON_NAME)))

butterflies$BNM.QC <- ifelse(butterflies$COMMON_NAME %in% spec.to.use, T,F)
butterflies$BMS.QC.t20 <- ifelse(butterflies$COMMON_NAME %in% BMS_t20_good$COMMON_NAME, T,F)
butterflies$BMS.QC.t18 <- ifelse(butterflies$COMMON_NAME %in% BMS_t18_good$COMMON_NAME, T,F)
butterflies$BMS.QC.t15 <- ifelse(butterflies$COMMON_NAME %in% BMS_t15_good$COMMON_NAME, T,F)
butterflies$BMS.QC.t10 <- ifelse(butterflies$COMMON_NAME %in% BMS_t10_good$COMMON_NAME, T,F)


butterflies$USE <- as.factor(ifelse(butterflies$BNM.QC == T & butterflies$BMS.QC.t20 == T, "Use (threshold 20)",
                             ifelse(butterflies$BNM.QC == T & butterflies$BMS.QC.t20 == F & butterflies$BMS.QC.t18 == T, "Use? (threshold 18)",
                             ifelse(butterflies$BNM.QC == T & butterflies$BMS.QC.t18 == F & butterflies$BMS.QC.t15 == T, "Use? (threshold 15)",
                             ifelse(butterflies$BNM.QC == T & butterflies$BMS.QC.t15 == F & butterflies$BMS.QC.t10 == T, "Exclude? (threshold 10)",
                             ifelse(butterflies$BNM.QC == F & butterflies$COMMON_NAME %in% migrants, "Exclude (migrant)",
                             ifelse(butterflies$BNM.QC == F & butterflies$COMMON_NAME %in% under.recorded, "Exclude (small distribution)",
                             ifelse(butterflies$BNM.QC == F & butterflies$COMMON_NAME %in% elevs.drop, "Exclude (upland species)",
                             "Exclude (no BMS data)"))))))))

butterflies$NORTHERLY <- ifelse(butterflies$COMMON_NAME %in% margins.spec.drop, T, F)


# write out this list
write.csv(butterflies, "Data/butterflyspeciestouse.csv", row.names = F)



### extract data for good BMS populations ####

# we want to use populations for species labelled "Use (threshold 20)", "Use (threshold 18)" and "Use (threshold 15)"

butterflies.to.use.spec <- butterflies[which(butterflies$USE %in% c("Use (threshold 20)","Use? (threshold 18)","Use? (threshold 15)")), ]
butterflies.to.use <- levels(droplevels(butterflies.to.use.spec$COMMON_NAME))

# now pull out the good sites for these species
# these are held in BMS_t15 - the lowest acceptable threshold for inclusion

popns.to.use <- BMS_t15[which(BMS_t15$COMMON_NAME %in% butterflies.to.use), ]

# for each of these populations we want to extract the raw data
# we need a single column for population to do this, which we can create by combining species & site names

popns.to.use$POPULATION <- paste(popns.to.use$COMMON_NAME, popns.to.use$SITE, sep=".")

# now repeat this for the raw data, which is in BMS_raw

BMS_years$POPULATION <- paste(BMS_years$COMMON_NAME, BMS_years$SITE, sep=".")

# and use this to extract only the data from good populations

BMS_popns <- BMS_years[which(BMS_years$POPULATION %in% popns.to.use$POPULATION), ]

## now we want a final cleaning step - we have taken forward species with at least 3 populations, where each population is recorded in at least 15 years
# but we additionally want to take forward only species that were recorded *somewhere* in all 20 years

BMS_spec.years <- ddply(BMS_popns, .(YEAR,SCI_NAME,COMMON_NAME), summarise, POPULATIONS = length(PEAKDAY))

BMS_specs <- ddply(BMS_spec.years, .(SCI_NAME,COMMON_NAME), summarise, YEARS = length(POPULATIONS))

BMS_specs.good <- BMS_specs[which(BMS_specs$YEARS==20),]
BMS.butterflies <- levels(droplevels(BMS_specs.good$COMMON_NAME))

# happily this means that we retain *all* species for the butterflies

BMS_popns <- BMS_popns[which(BMS_popns$COMMON_NAME %in% BMS.butterflies), ]

# now we have our final dataset, we can start producing outputs
# first, though, we want to extract the average abundance per count (rather than raw annual abundance)
# this is to partly address the possible issue of heavily-monitored transects appearing to have higher abundance than occasionally-monitored ones

BMS_popns$MEAN.ABUND <- BMS_popns$COUNT/BMS_popns$RECS


# let's generate some summary stats about how many BMS populations and how many population*year records go into the final dataset
# first, the overall stat of how many sites in total are included
length(levels(droplevels(BMS_popns$SITE)))

# and now some species-level stats

summary.stats <- data.frame(SCI_NAME = factor(),
                            COMMON_NAME = factor(),
                            SITES = numeric(),
                            SITE.YEARS = numeric(),
                            RECORDS = numeric(),
                            INDIVIDUALS = numeric())

for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
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


for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Phenology/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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

## once again, we have a pattern where most species are getting earlier - half a day per year for White Admiral and Small Heath
# and a couple are doing the opposite (Wall, Small Copper)


### now we want to replicate this process for abundance data from the same set of BMS populations
# (we will have good abundance )
# this data is already read in - 
summary(BMS_popns$MEAN.ABUND)

abund.slopes <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           ABUND.EXP.SLOPE = numeric(),
                           ABUND.SLOPE = numeric(),
                           ABUND.SE = numeric(),
                           ABUND.CHI = numeric(),
                           ABUND.P = numeric())


for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
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
  
  png(paste0("Data/Derived/Abundance/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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
  abund.slopes <- rbind(abund.slopes, out)
  
}

abund.slopes$ABUND.EXP.SLOPE <- as.numeric(as.character(abund.slopes$ABUND.EXP.SLOPE))
abund.slopes$ABUND.SLOPE <- as.numeric(as.character(abund.slopes$ABUND.SLOPE))
abund.slopes$ABUND.SE <- as.numeric(as.character(abund.slopes$ABUND.SE))
abund.slopes$ABUND.CHI <- as.numeric(as.character(abund.slopes$ABUND.CHI))
abund.slopes$ABUND.P <- as.numeric(as.character(abund.slopes$ABUND.P))

summary(abund.slopes)

## we have a range from species in sharp decline (e.g. PB Fritillary) to strongly-increasing species (e.g. Silver-washed Fritillary)

## finally, do the same for the voltinism ratio
# here we have a much-reduced set of data because (a) there were more failures (where the GAM failed to return to 0 by the end of the year),
# and (b) there are loads of zeroes (which are univoltine emergences)
# we clearly want to chop out (a), and we also want rid of (b) because we are specifically interested in the relationship between gen 1 and gen 2+
# for some univoltine species this will leave little or no data, but that's fine as we're only going to analyse multivoltines anyway for this one
# however, we want to trim out things with much too little data as they just won't fit a line properly - whilst also retaining all true multivoltines

summary(BMS_popns$VOLTINISM.RATIO) # at the moment the failures are in and it's being treated as a factor!

BMS_volti_pre <- BMS_popns[which(BMS_popns$VOLTINISM.RATIO != "Fail"), ]
BMS_volti_pre$VOLTINISM.RATIO <- as.numeric(as.character(BMS_volti_pre$VOLTINISM.RATIO))
BMS_volti <- BMS_volti_pre[which(BMS_volti_pre$VOLTINISM.RATIO != 0), ]

summary(BMS_volti)

# read in a list of what voltinism things have

voltinism <- read.csv("Data/voltinism.csv", header = T)


# check amount of data going in to this analysis per species, compared to the total amount
volti.spec <- data.frame(cbind(summary(BMS_volti$COMMON_NAME),summary(BMS_popns$COMMON_NAME)))
colnames(volti.spec) <- c("VOLTI_RECS","ALL_RECS")
volti.spec$ratio <- volti.spec$VOLTI_RECS/volti.spec$ALL_RECS
volti.spec$COMMON_NAME <- rownames(volti.spec)

volti.spec <- merge(voltinism[,c(3,5)],volti.spec)


# we can see there's not much attrition for true multivoltines but very heavy attrition for true univoltines
# we can actually split the two almost perfectly by asking for <50% attrition - 
# the only incorrect assignment being Brimstone, which is univoltine but has a life-history that appears bivoltine (with midsummer aestivation)

volti.spec <- volti.spec[which(volti.spec$ratio > 0.5), ]

BMS.butterflies.volti <- as.character(volti.spec$COMMON_NAME)

# prep the analysis
# the data are very clearly overdispersed (log-normal)...
hist(BMS_volti$VOLTINISM.RATIO)
hist(log(BMS_volti$VOLTINISM.RATIO))


volti.slopes <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           VOLTI.EXP.SLOPE = numeric(),
                           VOLTI.SLOPE = numeric(),
                           VOLTI.SE = numeric(),
                           VOLTI.CHI = numeric(),
                           VOLTI.P = numeric())


for (x in BMS.butterflies.volti){
  spec <- BMS_volti[which(BMS_volti$COMMON_NAME==x), ]
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
  
  png(paste0("Data/Derived/Voltinism/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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

BMS_out_pre <- merge(pheno.slopes,abund.slopes)
BMS_out_sum <- merge(summary.stats,BMS_out_pre)
BMS_out <- merge(BMS_out_sum,volti.slopes, all=T)



# and merge this into the BNM derived datasets
# some species might have different sci names between the two datasets
# the BNM one is more up-to-date so ditch the sci name column from the BMS dataset before merging
BMS_BNM <- merge(BMS_out[,-1], derived.BNM)
BMS_BNM_well <- merge(BMS_out[,-1], derived.BNM.well)
BMS_BNM_heavy <- merge(BMS_out[,-1], derived.BNM.heavy)


# and write out all versions
write.csv(BMS_BNM, "Data/Derived/BMS_BNM.csv", row.names = F)
write.csv(BMS_BNM_well, "Data/Derived/BMS_BNM_well.csv", row.names = F)
write.csv(BMS_BNM_heavy, "Data/Derived/BMS_BNM_heavy.csv", row.names = F)




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

for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
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
    
    ifelse(!dir.exists(paste0("Data/Derived/Intraspecific/Phenology/Butterflies/",x,"/")), 
           dir.create(paste0("Data/Derived/Intraspecific/Phenology/Butterflies/",x,"/"), recursive = T),
           FALSE)
    
    png(paste0("Data/Derived/Intraspecific/Phenology/Butterflies/",x,"/",y,".png"), width = 800, height = 800, units = "px", bg = "white")
    
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

for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
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
    
    ifelse(!dir.exists(paste0("Data/Derived/Intraspecific/Abundance/Butterflies/",x,"/")), 
           dir.create(paste0("Data/Derived/Intraspecific/Abundance/Butterflies/",x,"/"), recursive = T),
           FALSE)
    
    png(paste0("Data/Derived/Intraspecific/Abundance/Butterflies/",x,"/",y,".png"), width = 800, height = 800, units = "px", bg = "white")
    
    plot(log(site$MEAN.ABUND) ~ site$YEAR,
         ylim = c(0,10),
         xlab = "YEAR", ylab = "Log (mean no. recorded individuals per transect)")
    
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

write.csv(intraspecific, "Data/Derived/BMS_intraspec.csv", row.names = F)




#### temporal ####

# we now want to divide up the dataset into two 12-year periods, with a slight overlap (1995-2006 & 2003-2014)
# and see whether the pattern changes between them

# first divide up the dataset
inter1 <- 1995:2006
inter2 <- 2003:2014

BMS_inter1 <- BMS_popns[which(BMS_popns$YEAR %in% inter1), ]
BMS_inter2 <- BMS_popns[which(BMS_popns$YEAR %in% inter2), ]


# let's generate some summary stats about how many BMS populations and how many population*year records go into the final dataset
# interval 1
summary.stats.inter1 <- data.frame(SCI_NAME = factor(),
                            COMMON_NAME = factor(),
                            SITES = numeric(),
                            SITE.YEARS = numeric(),
                            RECORDS = numeric(),
                            INDIVIDUALS = numeric())

for (x in BMS.butterflies){
  spec <- BMS_inter1[which(BMS_inter1$COMMON_NAME==x), ]
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

for (x in BMS.butterflies){
  spec <- BMS_inter2[which(BMS_inter2$COMMON_NAME==x), ]
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


for (x in BMS.butterflies){
  spec <- BMS_inter1[which(BMS_inter1$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval1/Phenology/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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


for (x in BMS.butterflies){
  spec <- BMS_inter2[which(BMS_inter2$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Temporal/Interval2/Phenology/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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




### now we want to replicate this process for abundance data from the same set of BMS populations
# interval 1
abund.slopes.inter1 <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           ABUND.EXP.SLOPE = numeric(),
                           ABUND.SLOPE = numeric(),
                           ABUND.SE = numeric(),
                           ABUND.CHI = numeric(),
                           ABUND.P = numeric())


for (x in BMS.butterflies){
  spec <- BMS_inter1[which(BMS_inter1$COMMON_NAME==x), ]
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
  
  png(paste0("Data/Derived/Temporal/Interval1/Abundance/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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


for (x in BMS.butterflies){
  spec <- BMS_inter2[which(BMS_inter2$COMMON_NAME==x), ]
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
  
  png(paste0("Data/Derived/Temporal/Interval2/Abundance/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
  plot(log(spec$COUNT) ~ spec$YEAR,
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
write.csv(inters.derived, "Data/Derived/BMS_intervals.csv", row.names = F)


### site details ####
# I want to assess the annual climatic conditions at each of our sites, using GDD (growing degree days)
# so I need to put together a list of the sites that we're actually using, and attach grid references and target dates (for the GDD)

# first gather a list of sites
summary(BMS_popns)

BMS_sites <- ddply(BMS_popns, .(SITE), summarise,
                   SPECIES = nlevels(droplevels(SCI_NAME)))


# and generate a dataframe containing all combinations of site, year and three target dates: meteorological spring, summer and autumn
# these are Julian day 60, 152, 244

GDD_targets <- expand.grid(SITENAME = levels(droplevels(BMS_sites$SITE)),
                           YEAR = c(1995:2014),
                           JULIAN_DAY = c(60,152,244))

# we want to quickly adjust these for leap years

GDD_targets$JULIAN_DAY <- ifelse(GDD_targets$YEAR/4 == round(GDD_targets$YEAR/4), 1+GDD_targets$JULIAN_DAY, GDD_targets$JULIAN_DAY)


# now read in a dataframe containing details of all the sites
site_details <- read.csv("../../Data/UKBMS/sites_full.csv", header = T)


GDD_target_sites <- merge(GDD_targets,site_details, all.x = T)

# we can now write this out

write.csv(GDD_target_sites, "Data/UKBMS_sites.csv", row.names = F)


# for the purpose of plotting a figure, we want a version of the BMS_sites dataset with voltinism included
# merge the voltinism into the BMS_popns

BMS_popns_volti <- merge(BMS_popns, voltinism[,c(1,3,5)])

BMS_sites_volti <- ddply(BMS_popns_volti, .(SITE,VOLTINISM), summarise,
                         SPECIES = nlevels(droplevels(SCI_NAME)), .drop=FALSE)

BMS_sites_volti <- BMS_sites_volti[which(BMS_sites_volti$SITE %in% levels(droplevels(BMS_popns_volti$SITE))), ]

BMS_sites_volti$SITENAME <- BMS_sites_volti$SITE

BMS_sites_uv <- ddply(BMS_sites_volti[which(BMS_sites_volti$VOLTINISM == 1), ], .(SITENAME), summarise,
                      UNIVOLTINES = sum(SPECIES))

BMS_sites_mv <- ddply(BMS_sites_volti[which(BMS_sites_volti$VOLTINISM == 2), ], .(SITENAME), summarise,
                      MULTIVOLTINES = sum(SPECIES))

BMS_sites_count <- merge(BMS_sites_uv, BMS_sites_mv)

BMS_sites_count$SPECIES <- BMS_sites_count$UNIVOLTINES + BMS_sites_count$MULTIVOLTINES


# now merge in the site details
BMS_site_richness <- merge(BMS_sites_count, site_details)

# and write it out
write.csv(BMS_site_richness, "Data/UKBMS_sitecounts.csv", row.names = F)


### yearly analysis ####

# as well as these trends over time, I also want to extract a range of annual site-level variables for analysis;
# we want this to include phenology, overall abundance, and voltinism, all of which are in the raw data, so
# essentially this just means exporting the trimmed version of the dataset containing only good sites (RIS_popns)

summary(BMS_popns)

write.csv(BMS_popns, file = "Data/UKBMS_annual.csv", row.names = F)






### immediate response ####
# I want to assess whether an early emergence in year t leads to an increase in abundance in t+1
# this will require generation of another derived dataset - for each site * species * year combination,
# we need the emergence date and the abundance change to the next year

combos <- levels(droplevels(BMS_popns$siteXspecies))

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
  combo <- BMS_popns[which(BMS_popns$siteXspecies==x), ]
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
write.csv(indices, "Data/Derived/BMS_indices.csv", row.names = F)


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

for (x in BMS.butterflies){
  spec <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
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
write.csv(indices, "Data/Derived/BNM_NRMindices.csv", row.names = F)








#### exemplars ####
# now we want to plot a nicer job of each of the figures showing trends, for a pair of exemplar species:
# Adonis Blue and Large Skipper
# we'll use the .heavy BNM data

# let's do the whole lot for one species, then the whole lot for the other
### first, Adonis Blue

## let's start with a map plot showing change in distribution and northern range margin

# first, extract out the Adonis Blue data from the dataframes containing observed hectads and annual range margins

hectads.heavy.ab <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME=="Adonis Blue"), ]
margins.heavy.ab <- margins.heavy[which(margins.heavy$COMMON_NAME=="Adonis Blue"), ]

# extract out the years of interest
ab.1995.hecs <- hectads.heavy.ab[which(hectads.heavy.ab$YEAR %in% 1995), ]
ab.2014.hecs <- hectads.heavy.ab[which(hectads.heavy.ab$YEAR %in% 2014), ]

# merge them with all=F (default) to get hectads recorded in both
ab.both.hecs <- merge(ab.1995.hecs[,c(1:3,8:9)],ab.2014.hecs[,c(1:3,8:9)])


# 
set.mine <- set.British.Isles$Object[c(1:84)]

svg("Results/Plots/Exemplars/AdonisBlue/Distribution.svg",width = 12, height = 12,pointsize=12, bg = "white")

blighty(place="set.mine", set=FALSE)
points((5+ab.1995.hecs$EASTING/1000),(5+ab.1995.hecs$NORTHING/1000),pch=16,col="red")
points((5+ab.2014.hecs$EASTING/1000),(5+ab.2014.hecs$NORTHING/1000),pch=16,col="grey50")
points((5+ab.both.hecs$EASTING/1000),(5+ab.both.hecs$NORTHING/1000),pch=16,col="black")
abline(h=(margins.heavy.ab[which(margins.heavy.ab$YEAR==1995),4]/1000), col="red", lwd = 8)
abline(h=(margins.heavy.ab[which(margins.heavy.ab$YEAR==2014),4]/1000), col="grey50", lwd = 8)

dev.off()
 

## next, we want the trends in abundance
ab.abund <- BMS_popns[which(BMS_popns$COMMON_NAME=="Adonis Blue"), ]

test.ab.ab <- glmer(COUNT ~ YEAR + (1|SITE), data = ab.abund,
                 family = poisson (link = "log"))
  
coef.ab.ab <- fixef(test.ab.ab)
  

g.ab.ab <- ggplot(ab.abund)+
  geom_point(aes(x = YEAR, y = log(COUNT)), colour = "grey60")+
  geom_line(aes(x = YEAR, y = log(COUNT), group = SITE), colour = "grey80")+
  geom_abline(intercept = coef.ab.ab[1], slope = coef.ab.ab[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(0,10)+
  xlab("Year") + ylab("Log (abundance)")

g.ab.ab

ggsave("Results/Plots/Exemplars/AdonisBlue/Abundance.svg", plot = g.ab.ab, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## ...and phenology
test.phen.ab <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = ab.abund)

coef.phen.ab <- fixef(test.phen.ab)


g.phen.ab <- ggplot(ab.abund)+
  geom_point(aes(x = YEAR, y = PEAKDAY), colour = "grey60")+
  geom_line(aes(x = YEAR, y = PEAKDAY, group = SITE), colour = "grey80")+
  geom_abline(intercept = coef.phen.ab[1], slope = coef.phen.ab[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(100,200)+
  xlab("Year") + ylab("Julian day of emergence")

g.phen.ab

ggsave("Results/Plots/Exemplars/AdonisBlue/Emergence.svg", plot = g.phen.ab, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## finally, we want an overall depiction of voltinism across all (included) sites and years
## read in the daily data

dailies <- read.csv("../../Data/UKBMS/DailyCounts.csv", header=T)
summary(dailies)

# this gives the abundance index for each combination of species*site*year - we have information about those sites and species in separate files

species <- read.csv("../../Data/UKBMS/species.txt", header=TRUE)
summary(species)

# sites has to be treated a little differently - some site names have commas in them, so it can't simply be read in as a .csv file
# to fix this I've previously opened the sites .txt file in Notepad++, and changed all instances of ", " to "; "
sites <- read.csv("../../Data/UKBMS/sites.txt")
summary(sites)

# we want to take some means amongst all these, but some label variables are currently numeric (e.g. sites$SITENO)
# some columns containing the same information also have non-matching names
# let's make everything that needs to be a factor into a factor and rename where necessary
dailies$SPECIES <- as.factor(dailies$SPECIES)
dailies$SITE <- as.factor(dailies$SITE)
dailies$YEAR <- as.factor(dailies$YEAR)

species$SPECIES <- as.factor(species$BMSCODE)
species <- species[,-1]

sites$SITE <- as.factor(sites$SITENO)
sites <- sites[,2:4]

# some site names have a / in them which really screws up naming of files further down, so let's swap it out
sites$SITENAME <- gsub('/', '-', sites$SITENAME)

# we need to do the same for Thymelicus lineola/sylvestris in the species frame
species$SCI_NAME <- gsub('/', '-', species$SCI_NAME)

summary(dailies)
summary(species)
summary(sites)

# for inspection, let's merge species names into the main indices
dailies_merged <- merge(species,dailies)

## trim down to years and sites of interest
dailies_merged <- dailies_merged[which(dailies_merged$YEAR %in% 1995:2014), ]

sites.ab <- sites[which(sites$SITENAME %in% levels(droplevels(ab.abund$SITE))), ]


## then merge together, which will restrict it to *just* those sites and years, dropping everything else
dailies_merged.ab <- merge(dailies_merged,sites.ab)

# now we need to do some wrangling of the date
dailies_merged.ab$CalDate <- paste(dailies_merged.ab$DAY, dailies_merged.ab$MONTH, dailies_merged.ab$YEAR, sep="/")
dailies_merged.ab$CalDate <- as.Date(dailies_merged.ab$CalDate, format = "%d/%m/%Y")

# set up a reference date
dailies_merged.ab$YEAR <- as.numeric(as.character(dailies_merged.ab$YEAR))
dailies_merged.ab$refdate <- dmy(paste0("31-12-",(dailies_merged.ab$YEAR-1)))

dailies_merged.ab$julianday <- as.numeric(difftime(dailies_merged.ab$CalDate, dailies_merged.ab$refdate))

# add in whether a record is presence or absence
dailies_merged.ab$OBS <- ifelse(dailies_merged.ab$COUNT == 0, 0, 1)


# restrict only to species of interest
dailies_merged.ab <- dailies_merged.ab[which(dailies_merged.ab$COMMON_NAME=="Adonis Blue"), ]

# fit the gam
pheno.ab <- gam(COUNT ~ s(julianday), data = dailies_merged.ab,
                family = "poisson", method = "REML", gamma = 1)

# plot the gam
pred.ab <- data.frame(day = c(1:365))
pred.ab$pred <- predict(pheno.ab, type="response",newdata=data.frame(julianday=1:365))


plot(pred.ab$pred, type = "l",
     xlab = "Julian day",
     ylab = "Predicted abundance")

g.v.ab <- ggplot(pred.ab)+
  geom_point(data = dailies_merged.ab[which(dailies_merged.ab$COUNT > 0), ],
    aes(x = julianday, y = COUNT),
    colour = "grey60")+
  scale_y_continuous(limits = c(0,150))+
  geom_line(aes(x = day, y = pred), colour = "royalblue", size= 1.25)+
  theme_bw()+
  xlab("Julian day") + ylab("Abundance")

g.v.ab

ggsave("Results/Plots/Exemplars/AdonisBlue/Voltinism.svg", plot = g.v.ab, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


# these are all the necessary plots for Adonis Blue, so now let's repeat all of this for...

### Large Skipper
## let's start with a map plot showing change in distribution and northern range margin

# first, extract out the Large Skipper data from the dataframes containing observed hectads and annual range margins

hectads.heavy.ls <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME=="Large Skipper"), ]
margins.heavy.ls <- margins.heavy[which(margins.heavy$COMMON_NAME=="Large Skipper"), ]

# extract out the years of interest
ls.1995.hecs <- hectads.heavy.ls[which(hectads.heavy.ls$YEAR %in% 1995), ]
ls.2014.hecs <- hectads.heavy.ls[which(hectads.heavy.ls$YEAR %in% 2014), ]

# merge them with all=F (default) to get hectads recorded in both
ls.both.hecs <- merge(ls.1995.hecs[,c(1:3,8:9)],ls.2014.hecs[,c(1:3,8:9)])


# 
svg("Results/Plots/Exemplars/LargeSkipper/Distribution.svg",width = 12, height = 12,pointsize=12, bg = "white")

blighty(place="set.mine", set=FALSE)
points((5+ls.1995.hecs$EASTING/1000),(5+ls.1995.hecs$NORTHING/1000),pch=16,col="red")
points((5+ls.2014.hecs$EASTING/1000),(5+ls.2014.hecs$NORTHING/1000),pch=16,col="grey50")
points((5+ls.both.hecs$EASTING/1000),(5+ls.both.hecs$NORTHING/1000),pch=16,col="black")
abline(h=(margins.heavy.ls[which(margins.heavy.ls$YEAR==1995),4]/1000), col="red", lwd = 8)
abline(h=(margins.heavy.ls[which(margins.heavy.ls$YEAR==2014),4]/1000), col="grey50", lwd = 8)

dev.off()


## next, we want the trends in abundance
ls.abund <- BMS_popns[which(BMS_popns$COMMON_NAME=="Large Skipper"), ]

test.ab.ls <- glmer(COUNT ~ YEAR + (1|SITE), data = ls.abund,
                    family = poisson (link = "log"))

coef.ab.ls <- fixef(test.ab.ls)


g.ab.ls <- ggplot(ls.abund)+
  geom_point(aes(x = YEAR, y = log(COUNT)), colour = "grey60")+
  geom_line(aes(x = YEAR, y = log(COUNT), group = SITE), colour = "grey80")+
  geom_abline(intercept = coef.ab.ls[1], slope = coef.ab.ls[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(0,10)+
  xlab("Year") + ylab("Log (abundance)")

g.ab.ls

ggsave("Results/Plots/Exemplars/LargeSkipper/Abundance.svg", plot = g.ab.ls, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## ...and phenology
test.phen.ls <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = ls.abund)

coef.phen.ls <- fixef(test.phen.ls)


g.phen.ls <- ggplot(ls.abund)+
  geom_point(aes(x = YEAR, y = PEAKDAY), colour = "grey60")+
  geom_line(aes(x = YEAR, y = PEAKDAY, group = SITE), colour = "grey80")+
  geom_abline(intercept = coef.phen.ls[1], slope = coef.phen.ls[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(125,225)+
  xlab("Year") + ylab("Julian day of emergence")

g.phen.ls

ggsave("Results/Plots/Exemplars/LargeSkipper/Emergence.svg", plot = g.phen.ls, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## finally, we want an overall depiction of voltinism across all (included) sites and years
## we already have the right data for the years of interest
## trim down to sites of interest
sites.ls <- sites[which(sites$SITENAME %in% levels(droplevels(ls.abund$SITE))), ]

## then merge together, which will restrict it to *just* those sites and years, dropping everything else
dailies_merged.ls <- merge(dailies_merged,sites.ls)

# now we need to do some wrangling of the date
dailies_merged.ls$CalDate <- paste(dailies_merged.ls$DAY, dailies_merged.ls$MONTH, dailies_merged.ls$YEAR, sep="/")
dailies_merged.ls$CalDate <- as.Date(dailies_merged.ls$CalDate, format = "%d/%m/%Y")

# set up a reference date
dailies_merged.ls$YEAR <- as.numeric(as.character(dailies_merged.ls$YEAR))
dailies_merged.ls$refdate <- dmy(paste0("31-12-",(dailies_merged.ls$YEAR-1)))

dailies_merged.ls$julianday <- as.numeric(difftime(dailies_merged.ls$CalDate, dailies_merged.ls$refdate))

# add in whether a record is presence or absence
dailies_merged.ls$OBS <- ifelse(dailies_merged.ls$COUNT == 0, 0, 1)


# restrict only to species of interest
dailies_merged.ls <- dailies_merged.ls[which(dailies_merged.ls$COMMON_NAME=="Large Skipper"), ]

# fit the gam
pheno.ls <- gam(COUNT ~ s(julianday), data = dailies_merged.ls,
                family = "poisson", method = "REML", gamma = 1)

# plot the gam
pred.ls <- data.frame(day = c(1:365))
pred.ls$pred <- predict(pheno.ls, type="response",newdata=data.frame(julianday=1:365))


plot(pred.ls$pred, type = "l",
     xlab = "Julian day",
     ylab = "Predicted abundance")

g.v.ls <- ggplot(pred.ls)+
  geom_point(data = dailies_merged.ls[which(dailies_merged.ls$COUNT > 0), ],
             aes(x = julianday, y = COUNT),
             colour = "grey60")+
  scale_y_continuous(limits = c(0,150))+
  geom_line(aes(x = day, y = pred), colour = "royalblue", size= 1.25)+
  theme_bw()+
  xlab("Julian day") + ylab("Abundance")

g.v.ls

ggsave("Results/Plots/Exemplars/LargeSkipper/Voltinism.svg", plot = g.v.ls, device = svg, width = 100, height = 70, units = "mm", limitsize = T)



## try stitching these together
m.exemplars <- grid.arrange(g.v.ls,g.v.ab + ylab(""),
                            g.phen.ls + xlab(""),g.phen.ab + xlab("") + ylab(""),
                            g.ab.ls,g.ab.ab + ylab(""),
                            ncol = 2)


ggsave("Results/Plots/Exemplars/StitchedPlot.svg", plot = m.exemplars, device = svg, width = 200, height = 210, units = "mm", limitsize = T)






# having done all of that, we now want to recreate this plot for a different pair of species:
# Small Blue and Silver-studded Blue
# we'll use the .heavy BNM data

# let's do the whole lot for one species, then the whole lot for the other
### first, Small Blue

## let's start with a map plot showing change in distribution and northern range margin

# first, extract out the Adonis Blue data from the dataframes containing observed hectads and annual range margins

hectads.heavy.sb <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME=="Small Blue"), ]
margins.heavy.sb <- margins.heavy[which(margins.heavy$COMMON_NAME=="Small Blue"), ]

# extract out the years of interest
sb.1995.hecs <- hectads.heavy.sb[which(hectads.heavy.sb$YEAR %in% 1995), ]
sb.2014.hecs <- hectads.heavy.sb[which(hectads.heavy.sb$YEAR %in% 2014), ]

# merge them with all=F (default) to get hectads recorded in both
sb.both.hecs <- merge(sb.1995.hecs[,c(1:3,8:9)],sb.2014.hecs[,c(1:3,8:9)])


# 
set.mine <- set.British.Isles$Object[c(1:84)]
set.mine.res <- set.British.Isles$Object[c(1:3,16:58,77:81)]

svg("Results/Plots/Exemplars/SmallBlue/Distribution.svg",width = 12, height = 12,pointsize=12, bg = "white")

blighty(place="set.mine.res", set=FALSE)
points((5+sb.1995.hecs$EASTING/1000),(5+sb.1995.hecs$NORTHING/1000),pch=16,col="goldenrod")
points((5+sb.2014.hecs$EASTING/1000),(5+sb.2014.hecs$NORTHING/1000),pch=16, col="royalblue")
points((5+sb.both.hecs$EASTING/1000),(5+sb.both.hecs$NORTHING/1000),pch=16,col="black")
abline(h=(margins.heavy.sb[which(margins.heavy.sb$YEAR==1995),4]/1000), col="goldenrod", lwd = 8)
abline(h=(margins.heavy.sb[which(margins.heavy.sb$YEAR==2014),4]/1000), col="black", lty = 2, lwd = 8)

dev.off()


## next, we want the trends in abundance
sb.abund <- BMS_popns[which(BMS_popns$COMMON_NAME=="Small Blue"), ]

test.ab.sb <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = sb.abund)

coef.ab.sb <- fixef(test.ab.sb)


# we want confidence intervals around the line, so we need to predict data for these

errors.ab.sb <- data.frame(effect(c("YEAR"), test.ab.sb, xlevels=20))



g.ab.sb <- ggplot(sb.abund)+
  geom_ribbon(data = errors.ab.sb,
              aes(x = YEAR, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey80")+
  geom_point(aes(x = YEAR, y = log(MEAN.ABUND)), colour = "grey50")+
  geom_line(aes(x = YEAR, y = log(MEAN.ABUND), group = SITE), colour = "grey70")+
  geom_abline(intercept = coef.ab.sb[1], slope = coef.ab.sb[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(-2,4.5)+
  xlab("Year") + ylab("Log (mean abundance)")

g.ab.sb

ggsave("Results/Plots/Exemplars/SmallBlue/Abundance.svg", plot = g.ab.sb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## ...and phenology
test.phen.sb <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = sb.abund)

coef.phen.sb <- fixef(test.phen.sb)


errors.phen.sb <- data.frame(effect(c("YEAR"), test.phen.sb, xlevels=20))


g.phen.sb <- ggplot(sb.abund)+
  geom_ribbon(data = errors.phen.sb,
              aes(x = YEAR, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey80")+
  geom_point(aes(x = YEAR, y = PEAKDAY), colour = "grey50")+
  geom_line(aes(x = YEAR, y = PEAKDAY, group = SITE), colour = "grey70")+
  geom_abline(intercept = coef.phen.sb[1], slope = coef.phen.sb[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(125,250)+
  xlab("Year") + ylab("Julian day of emergence")

g.phen.sb

ggsave("Results/Plots/Exemplars/SmallBlue/Emergence.svg", plot = g.phen.sb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## finally, we want an overall depiction of voltinism across all (included) sites and years

sites.sb <- sites[which(sites$SITENAME %in% levels(droplevels(sb.abund$SITE))), ]


## then merge together, which will restrict it to *just* those sites and years, dropping everything else
dailies_merged.sb <- merge(dailies_merged,sites.sb)

# now we need to do some wrangling of the date
dailies_merged.sb$CalDate <- paste(dailies_merged.sb$DAY, dailies_merged.sb$MONTH, dailies_merged.sb$YEAR, sep="/")
dailies_merged.sb$CalDate <- as.Date(dailies_merged.sb$CalDate, format = "%d/%m/%Y")

# set up a reference date
dailies_merged.sb$YEAR <- as.numeric(as.character(dailies_merged.sb$YEAR))
dailies_merged.sb$refdate <- dmy(paste0("31-12-",(dailies_merged.sb$YEAR-1)))

dailies_merged.sb$julianday <- as.numeric(difftime(dailies_merged.sb$CalDate, dailies_merged.sb$refdate))

# add in whether a record is presence or absence
dailies_merged.sb$OBS <- ifelse(dailies_merged.sb$COUNT == 0, 0, 1)


# restrict only to species of interest
dailies_merged.sb <- dailies_merged.sb[which(dailies_merged.sb$COMMON_NAME=="Small Blue"), ]

# fit the gam
pheno.sb <- gam(COUNT ~ s(julianday), data = dailies_merged.sb,
                family = "poisson", method = "REML", gamma = 1)

# plot the gam
pred.sb <- data.frame(day = c(1:365))
pred.sb$pred <- predict(pheno.sb, type="response",newdata=data.frame(julianday=1:365))


plot(pred.sb$pred, type = "l",
     xlab = "Julian day",
     ylab = "Predicted abundance")

g.v.sb <- ggplot(pred.sb)+
  geom_point(data = dailies_merged.sb[which(dailies_merged.sb$COUNT > 0), ],
             aes(x = julianday, y = COUNT),
             colour = "grey60")+
  scale_y_continuous(limits = c(0,100))+
  geom_line(aes(x = day, y = pred), colour = "royalblue", size= 1.25)+
  theme_bw()+
  xlab("Julian day") + ylab("Abundance")

g.v.sb

ggsave("Results/Plots/Exemplars/SmallBlue/Voltinism.svg", plot = g.v.sb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


# these are all the necessary plots for Small Blue, so now let's repeat all of this for...

### Silver-studded Blue
## let's start with a map plot showing change in distribution and northern range margin

# first, extract out the SSB data from the dataframes containing observed hectads and annual range margins

hectads.heavy.ssb <- hectads.years.heavy[which(hectads.years.heavy$COMMON_NAME=="Silver-studded Blue"), ]
margins.heavy.ssb <- margins.heavy[which(margins.heavy$COMMON_NAME=="Silver-studded Blue"), ]

# extract out the years of interest
ssb.1995.hecs <- hectads.heavy.ssb[which(hectads.heavy.ssb$YEAR %in% 1995), ]
ssb.2014.hecs <- hectads.heavy.ssb[which(hectads.heavy.ssb$YEAR %in% 2014), ]

# merge them with all=F (default) to get hectads recorded in both
ssb.both.hecs <- merge(ssb.1995.hecs[,c(1:3,8:9)],ssb.2014.hecs[,c(1:3,8:9)])


# 
svg("Results/Plots/Exemplars/SilverStuddedBlue/Distribution.svg",width = 12, height = 12,pointsize=12, bg = "white")

blighty(place="set.mine.res", set=FALSE)
points((5+ssb.1995.hecs$EASTING/1000),(5+ssb.1995.hecs$NORTHING/1000),pch=16,col="goldenrod")
points((5+ssb.2014.hecs$EASTING/1000),(5+ssb.2014.hecs$NORTHING/1000),pch=16,col="royalblue")
points((5+ssb.both.hecs$EASTING/1000),(5+ssb.both.hecs$NORTHING/1000),pch=16,col="black")
abline(h=(margins.heavy.ssb[which(margins.heavy.ssb$YEAR==1995),4]/1000), col="goldenrod", lwd = 8)
abline(h=(margins.heavy.ssb[which(margins.heavy.ssb$YEAR==2014),4]/1000), col="black", lty = 2, lwd = 8)

dev.off()


## next, we want the trends in abundance
ssb.abund <- BMS_popns[which(BMS_popns$COMMON_NAME=="Silver-studded Blue"), ]

test.ab.ssb <- lmer(log(MEAN.ABUND) ~ YEAR + (1|SITE), data = ssb.abund)

coef.ab.ssb <- fixef(test.ab.ssb)

errors.ab.ssb <- data.frame(effect(c("YEAR"), test.ab.ssb, xlevels=20))


g.ab.ssb <- ggplot(ssb.abund)+
  geom_ribbon(data = errors.ab.ssb,
              aes(x = YEAR, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey80")+
  geom_point(aes(x = YEAR, y = log(MEAN.ABUND)), colour = "grey50")+
  geom_line(aes(x = YEAR, y = log(MEAN.ABUND), group = SITE), colour = "grey70")+
  geom_abline(intercept = coef.ab.ssb[1], slope = coef.ab.ssb[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(-2,4.5)+
  xlab("Year") + ylab("Log (mean abundance)")

g.ab.ssb

ggsave("Results/Plots/Exemplars/SilverStuddedBlue/Abundance.svg", plot = g.ab.ssb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## ...and phenology
test.phen.ssb <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = ssb.abund)

coef.phen.ssb <- fixef(test.phen.ssb)

errors.phen.ssb <- data.frame(effect(c("YEAR"), test.phen.ssb, xlevels=20))


g.phen.ssb <- ggplot(ssb.abund)+
  geom_ribbon(data = errors.phen.ssb,
              aes(x = YEAR, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey80")+
  geom_point(aes(x = YEAR, y = PEAKDAY), colour = "grey50")+
  geom_line(aes(x = YEAR, y = PEAKDAY, group = SITE), colour = "grey70")+
  geom_abline(intercept = coef.phen.ssb[1], slope = coef.phen.ssb[2], size = 1, linetype = "dashed")+
  theme_bw()+
  ylim(125,250)+
  xlab("Year") + ylab("Julian day of emergence")

g.phen.ssb

ggsave("Results/Plots/Exemplars/SilverStuddedBlue/Emergence.svg", plot = g.phen.ssb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)


## finally, we want an overall depiction of voltinism across all (included) sites and years
## we already have the right data for the years of interest
## trim down to sites of interest
sites.ssb <- sites[which(sites$SITENAME %in% levels(droplevels(ssb.abund$SITE))), ]

## then merge together, which will restrict it to *just* those sites and years, dropping everything else
dailies_merged.ssb <- merge(dailies_merged,sites.ssb)

# now we need to do some wrangling of the date
dailies_merged.ssb$CalDate <- paste(dailies_merged.ssb$DAY, dailies_merged.ssb$MONTH, dailies_merged.ssb$YEAR, sep="/")
dailies_merged.ssb$CalDate <- as.Date(dailies_merged.ssb$CalDate, format = "%d/%m/%Y")

# set up a reference date
dailies_merged.ssb$YEAR <- as.numeric(as.character(dailies_merged.ssb$YEAR))
dailies_merged.ssb$refdate <- dmy(paste0("31-12-",(dailies_merged.ssb$YEAR-1)))

dailies_merged.ssb$julianday <- as.numeric(difftime(dailies_merged.ssb$CalDate, dailies_merged.ssb$refdate))

# add in whether a record is presence or absence
dailies_merged.ssb$OBS <- ifelse(dailies_merged.ssb$COUNT == 0, 0, 1)


# restrict only to species of interest
dailies_merged.ssb <- dailies_merged.ssb[which(dailies_merged.ssb$COMMON_NAME=="Silver-studded Blue"), ]

# fit the gam
pheno.ssb <- gam(COUNT ~ s(julianday), data = dailies_merged.ssb,
                family = "poisson", method = "REML", gamma = 1)

# plot the gam
pred.ssb <- data.frame(day = c(1:365))
pred.ssb$pred <- predict(pheno.ssb, type="response",newdata=data.frame(julianday=1:365))


plot(pred.ssb$pred, type = "l",
     xlab = "Julian day",
     ylab = "Predicted abundance")

g.v.ssb <- ggplot(pred.ssb)+
  geom_point(data = dailies_merged.ssb[which(dailies_merged.ssb$COUNT > 0), ],
             aes(x = julianday, y = COUNT),
             colour = "grey60")+
  scale_y_continuous(limits = c(0,100))+
  geom_line(aes(x = day, y = pred), colour = "royalblue", size= 1.25)+
  theme_bw()+
  xlab("Julian day") + ylab("Abundance")

g.v.ssb

ggsave("Results/Plots/Exemplars/SilverStuddedBlue/Voltinism.svg", plot = g.v.ssb, device = svg, width = 100, height = 70, units = "mm", limitsize = T)



## try stitching these together
m.exemplars.2 <- grid.arrange(g.v.ssb,g.v.sb + ylab(""),
                            g.phen.ssb + xlab(""),g.phen.sb + xlab("") + ylab(""),
                            g.ab.ssb,g.ab.sb + ylab(""),
                            ncol = 2)


ggsave("Results/Plots/Exemplars/StitchedPlot2.svg", plot = m.exemplars.2, device = svg, width = 200, height = 210, units = "mm", limitsize = T)






### loop the voltinism plots to create one for every species
## finally, we want an overall depiction of voltinism across all (included) sites and years
## we already have the right data for the years of interest


# wrangle the dates
dailies_merged$CalDate <- paste(dailies_merged$DAY, dailies_merged$MONTH, dailies_merged$YEAR, sep="/")
dailies_merged$CalDate <- as.Date(dailies_merged$CalDate, format = "%d/%m/%Y")

# set up a reference date
dailies_merged$YEAR <- as.numeric(as.character(dailies_merged$YEAR))
dailies_merged$refdate <- dmy(paste0("31-12-",(dailies_merged$YEAR-1)))

dailies_merged$julianday <- as.numeric(difftime(dailies_merged$CalDate, dailies_merged$refdate))

# add in whether a record is presence or absence
dailies_merged$OBS <- ifelse(dailies_merged$COUNT == 0, 0, 1)


# now construct a loop to fit and plot out the GAM

for (x in levels(droplevels(BMS_popns$COMMON_NAME))){
  spec <- dailies_merged[which(dailies_merged$COMMON_NAME==x), ]
  
  abund <- BMS_popns[which(BMS_popns$COMMON_NAME==x), ]
  sites.spec <- sites[which(sites$SITENAME %in% levels(droplevels(abund$SITE))), ]
  
  ## then merge together, which will restrict it to *just* those sites and years, dropping everything else
  dailies_merged.spec <- merge(spec,sites.spec)
  
  
  # fit the gam
  pheno <- gam(COUNT ~ s(julianday), data = dailies_merged.spec,
               family = "poisson", method = "REML", gamma = 1)
  
  # plot the gam
  pred <- data.frame(day = c(1:365))
  pred$pred <- predict(pheno, type="response",newdata=data.frame(julianday=1:365))
  

  g.spec <- ggplot(pred)+
    geom_point(data = dailies_merged.spec[which(dailies_merged.spec$COUNT > 0), ], 
               aes(x = julianday, y = COUNT), 
               colour = "grey60")+
    scale_y_continuous(limits = c(0,100))+
    geom_line(aes(x = day, y = pred), colour = "royalblue", size = 1.25)+
    theme_bw()+
    xlab("Julian day") + ylab("Abundance")
  
  ggsave(paste0("Results/Plots/Exemplars/AllSpecies/Butterflies/",x,".png"), plot = g.spec, width = 100, height = 70, units = "mm", limitsize = F)
  
}



### post reviewers' comments ####

# one of the reviewers wants us to do a supplementary analysis using the full length of the available butterfly data
# i.e. for every population, rather than restricting to 1995-2014, calculate the trends from start to finish of recording


# we still have a list of all the populations to use, in 'popns.to.use'

# so now we return to the completely raw data, with all years included, to extract the relevant data


# now repeat this for the raw data, which is in BMS_raw

BMS_raw$POPULATION <- paste(BMS_raw$COMMON_NAME, BMS_raw$SITE, sep=".")

# and use this to extract only the data from good populations

BMS_allyears <- BMS_raw[which(BMS_raw$POPULATION %in% popns.to.use$POPULATION), ]

summary(BMS_allyears)

## now extract abundance and phenology trends for each of these populations, as before

BMS_allyears$MEAN.ABUND <- BMS_allyears$COUNT/BMS_allyears$RECS


# let's generate some summary stats about how many BMS populations and how many population*year records go into the final dataset
# first, the overall stat of how many sites in total are included (should be 110 as before)
length(levels(droplevels(BMS_allyears$SITE)))

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

for (x in BMS.butterflies){
  spec <- BMS_allyears[which(BMS_allyears$COMMON_NAME==x), ]
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


for (x in BMS.butterflies){
  spec <- BMS_allyears[which(BMS_allyears$COMMON_NAME==x), ]
  COMMON_NAME <- x
  SCI_NAME <- as.character(spec$SCI_NAME[1])
  print(x)
  
  test1 <- lmer(PEAKDAY ~ YEAR + (1|SITE), data = spec)
  
  PHENO.SLOPE <- -round(summary(test1)$coefficients[2,1], 2)
  PHENO.SE <- round(summary(test1)$coefficients[2,2], 2)
  PHENO.CHI <- round(drop1(test1, test = "Chi")[2,3], 2) 
  PHENO.P <- round(drop1(test1, test = "Chi")[2,4], 4)
  coef <- fixef(test1)
  
  png(paste0("Data/Derived/Phenology/FullRange/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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



### now we want to replicate this process for abundance data from the same set of BMS populations
# (we will have good abundance )
# this data is already read in - 
summary(BMS_allyears$MEAN.ABUND)

abund.slopes.ay <- data.frame(SCI_NAME = factor(),
                           COMMON_NAME = factor(),
                           ABUND.EXP.SLOPE = numeric(),
                           ABUND.SLOPE = numeric(),
                           ABUND.SE = numeric(),
                           ABUND.CHI = numeric(),
                           ABUND.P = numeric())


for (x in BMS.butterflies){
  spec <- BMS_allyears[which(BMS_allyears$COMMON_NAME==x), ]
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
  
  png(paste0("Data/Derived/Abundance/FullRange/Butterflies/",x,".png"), width = 800, height = 800, units = "px", bg = "white")
  
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

BMS_out_pre.ay <- merge(pheno.slopes.ay,abund.slopes.ay)
BMS_out.ay <- merge(summary.stats.ay,BMS_out_pre.ay)


# and write out all versions
write.csv(BMS_out.ay, "Data/Derived/BMS_ay.csv", row.names = F)



