#######################################################
####   Script for initial processing of RIS data   ####
#######################################################

##### goals #####

# we have the daily data for all Rothamsted traps in GB, for all moth species recorded
# this is held in a separate .csv file for each year

# we want to process this down into a single .csv file containing the FIRST, LAST, MEANDAY and SD phenology estimates 
# for every combination of site*species*year, much like we have for UKBMS


### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","reshape2","lubridate","mgcv")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### I'm going to do everything inside a massive loop, so it imports the data for each year, wrangles it, 
# and writes the site*species summaries to a single dataframe for eventual export


# prepare the loop

final.list <- NULL

# and run the loop!

for (year in 1980:2016){
  print(year)
  
  raw_data <- read.csv(paste0("../Data/RIS/RIS_raw/AllSitesAllMoths ",year,".csv"),header=T)
  
  # we have three things for 'species identity', a thing for 'site name', a date, an abundance, and a couple of excess columns which can be trimmed off
  raw_data <- raw_data[,c(1:5,7)]
  
  # the first thing to do is to convert the date - which is currently Julian - into a Julian day (i.e. 1-365, or 366 in a leap year):
  # this is a format counting outwards from Jan 1 in each year

  # we can do this conversion using the 'difftime' function in the 'lubridate' package

  # we'll need the CalDate column to be in R date format
  
  raw_data$CalDate <- as.Date(raw_data$CalDate, format = "%d/%m/%Y")
  
  # there is a slight problem which means we'll have to do this cleverly: 
  # lubridate uses a 'day 0' to count forwards/backwards from, and we can't overrule this
  # therefore we need to make this calculation from Dec 31st of the previous year for Jan 1st to be day 1
  
  # first we need to set up that reference date (i.e. Dec 31st of the preceding year)

  Dec31 <- paste0("31-12-",(year-1))
  
  ref_date <- dmy(Dec31)
  
  
  # now we can calculate Julian day
  
  raw_data$julianday <- as.numeric(difftime(raw_data$CalDate, ref_date))
  
  
  
  ## now we are ready to calculate some summary stats!
  
  summary <- ddply(raw_data, .(binomial, Code.RISgeneric, Common.Name, TrapName), summarise,
                   SD = sd(julianday), MEANDAY = weighted.mean(julianday, DailyCount),
                   FIRST = min(julianday), LAST = max(julianday),
                   COUNT = sum(DailyCount), OBS = length(julianday))

  
  colnames(summary) <- c("SCI_NAME","SPECIES","COMMON_NAME","SITE","SD","MEANDAY","FIRST","LAST","COUNT","OBS")

  summary$YEAR <- as.factor(year)
  
  summary[is.na(summary)] <- 0
  
  # we also want a list of all trap nights
  # (this is really important for putting zero catches into species-level estimates so that GAMs don't go wild with low data!)
  
  trapnightssummary <- ddply(raw_data, .(TrapName,julianday), summarise,
                             TotalCatch = sum(DailyCount))
  
  trapnightssummary$TrapNight <- as.factor(paste(trapnightssummary$julianday,trapnightssummary$TrapName,sep = "."))
  
  trapnights <- levels(trapnightssummary$TrapNight)
  
  # now, for each combination of site * species with sufficient observations we want to estimate the phenology of the first generation
  # using a GAM-fitting approach and finding the *first peak* of each GAM to identify the peak emergence of the first generation
  
  # first generate a vector of all site*species combinations from the 'summary' frame
  
  summary$siteXspecies <- factor(paste(summary$SCI_NAME,summary$SITE,sep="."))
  
  # chuck out combinations with too few observations
  summary.good <- summary[which(summary$OBS > 3), ]
  
  # create the vector to loop over
  summary.good$siteXspecies <- droplevels(summary.good$siteXspecies)
  
  loop.vec <- as.character(summary.good$siteXspecies)
  
  # add the same combination columns to the raw data that we'll be working from
  raw_data$siteXspecies <- factor(paste(raw_data$binomial,raw_data$TrapName,sep="."))
  raw_data$TrapNight <- factor(paste(raw_data$julianday,raw_data$TrapName,sep="."))
  
  # seed some output vectors
  RECS <- NULL
  peak.day <- NULL
  first.trough <- NULL
  second.peak <- NULL
  voltinism.ratio <- NULL
  first.gen <- NULL
  second.gen <- NULL
  current.sp <- NULL
  previous.sp <- "None"
  
  # now construct the loop
  for (combo in loop.vec){
    sxs.data <- raw_data[which(raw_data$siteXspecies == combo), ]
    current.sp <- as.character(sxs.data[[1,1]])
    current.site <- as.character(sxs.data[[1,4]])
    
    if (!(current.sp == previous.sp)){
      print(current.sp)
    }

    
    ### input days with traps but no observations (at recorded sites only)
    # first pick out recorded sites
    sxs.data$TrapName <- droplevels(sxs.data$TrapName)
    recsites <- levels(sxs.data$TrapName)
    
    # take only these from the trapnight summary
    alltraps <- trapnightssummary[which(trapnightssummary$TrapName %in% recsites), ]
    
    # merge these into the data
    sxs.data <- merge(sxs.data, alltraps, all=T)
    
    # make NAs in the DailyCount column into zeroes (i.e. there was a trap, but it didn't catch the species of interest)
    sxs.data$DailyCount[is.na(sxs.data$DailyCount)] <- 0
    
    # count the number of records
    recs <- nrow(sxs.data)
    RECS <- append(RECS, recs)
    
    # check if there are enough points (including zeroes) to fit a GAM - need at least 10
    if (nrow(sxs.data) < 10){
      peak1 <- "Fail"
      true.trough <- "Fail"
      peak2 <- "Fail"
    } else {
    
    
    ### fit a GAM to the data
    pheno <- gam(DailyCount ~ s(julianday),data=sxs.data,family="poisson",method="REML",gamma=1)
  
    # predict the abundance on each day of the year from the GAM
    pred <- data.frame(day = c(1:365))
    pred$pred <- predict(pheno,type="response",newdata=data.frame(julianday=1:365))
  
    # save the model output as a simple plot - in case we ever need to check any of the models
    ifelse(!dir.exists(paste0("Plots/GAMs/ForIona/RIS/",current.sp,"/",current.site,"/")), 
           dir.create(paste0("Plots/GAMs/ForIona/RIS/",current.sp,"/",current.site,"/"), recursive = T),
           FALSE)
    
    # very occasionally we might encounter infinite values, which we can't plot - so let's just make them very, very big
    
    pred$predplot <- ifelse(pred$pred > 100, 100, pred$pred)
    
    # there are also cases where the estimate wiggles about at extremely low values and there might be a false peak in amongst that
    # so let's cut these out first and keep only true emergences
    
    pred$pred <- ifelse(pred$pred < 0.0001, 0, pred$pred)
    
    ymax <- max(max(pred$predplot),max(sxs.data$DailyCount))
    
    png(paste0("Plots/GAMs/ForIona/RIS/",current.sp,"/",current.site,"/",year,".png"), width = 800, height = 800, units = "px", bg = "white")
    
    
    plot(pred$predplot, type = "l",
         xlab = "Julian day",
         ylab = "Predicted abundance",
         ylim = c(0,ymax))
    points(DailyCount ~ julianday, data=sxs.data)
    
    dev.off()
    
    # now find the first maximum in that GAM - 
    # run down every day sequentially until it finds a day that has a higher predicted abundance than the day before and the day after
    # sometimes the data just won't support a GAM that both goes up and comes back down - we don't want to use non-peaking cases
    # (in such cases the loop will run all the way to the end of the year)
  
    # n.b. although we have plotted them, we don't want to carry forward estimates from any GAM that doesn't start the year at zero
    # these GAMs indicate that the flight period overlapped with the start of annual recording at the site 
    # (i.e. there are not enough zero/small counts before the first record to get a good estimate of when the first emergence peaked) 
    
    # so let's place this part in an if/else that pulls out of any GAM that breaks that rule
  
    # n.b. this doesn't actually give our final estimate of peak day, but it is important to allow the code to run properly further down
    
      
    peak <- FALSE
    x <- 1
    
    if (!pred[1,2] < 1){
      peak1 <- "Fail"
      true.trough <- "Fail"
      peak2 <- "Fail"
    } else {
      # loop until true or broken...
      
      while (!peak) {
        x <- x+1
        if (x == 365){
          peak <- TRUE
          answer <- "Fail"
        } else {
          peak <- pred[x,2] > pred[x-1,2] & pred[x,2] > pred[x+1,2]
          answer <- x
        }
      }
      
      # we want to cut off this all here if there is no peak
      if (answer == "Fail"){
        peak1 <- "Fail"
        true.trough <- "Fail"
        peak2 <- "Fail"
      } else {
        
        
        # and using the same logic, find the first trough after that peak...
        # run down all days sequentially until either it finds a day that has a lower predicted abundance than the day before and the day after
        # or until it returns to zero (in which case the day after will have the same predicted abundance)
        
        # however, it is more complex than that
        # if there are two peaks (i.e. two generations), we want it to pick the first trough
        # but if there are more than two peaks (i.e. either 3 generations or a big mess that the GAM can't quite resolve) - 
        # we might want a tiny bit more sublety
        
        # first, count the number of peaks and troughs
        
        peaks <- NULL
        
        for (x in answer:364){
          if (pred[x,2] > pred[x-1,2] & pred[x,2] > pred[x+1,2]){    # identified a possible peak
            
            prev.trough <- FALSE
            y <- x
            
            while (!prev.trough) {
              y <- y-1
              if (y == 1){
                prev.trough <- TRUE
              } else {
                prev.trough <- pred[y,2] < pred[y-1,2] & pred[y,2] < pred[y+1,2]
              }
            }
            
            next.trough <- FALSE
            z <- x
            
            while (!next.trough) {
              z <- z+1
              if (z == 365){
                next.trough <- TRUE
              } else {
                next.trough <- pred[z,2] < pred[z-1,2] & pred[z,2] < pred[z+1,2]
              }
            }
            
            
            peakiness <- min((pred[x,2]/pred[y,2]),(pred[x,2]/pred[z,2]))
            
            if (peakiness > 1.25){
              peaks <- append(peaks, x)
            }
          }
        }
        
        
        if (pred[365,2] > pred[364,2]){
          peaks <- append(peaks, 365)
        }
        
        n.peaks <- length(peaks)
        
        
        # now identify all the troughs
        
        troughs <- NULL
        
        for (x in answer:364){
          if (pred[x,2] < pred[x-1,2] & pred[x,2] <= pred[x+1,2]){    # identified a possible trough
            
            prev.peak <- FALSE
            y <- x
            
            while (!prev.peak) {
              y <- y-1
              if (y == 1){
                prev.peak <- TRUE
              } else {
                prev.peak <- pred[y,2] > pred[y-1,2] & pred[y,2] > pred[y+1,2]
              }
            }
            
            next.peak <- FALSE
            z <- x
            
            while (!next.peak) {
              z <- z+1
              if (z == 365){
                next.peak <- TRUE
              } else {
                next.peak <- pred[z,2] > pred[z-1,2] & pred[z,2] > pred[z+1,2]
              }
            }
            
            
            if (pred[z,2] == 0){
              troughs <- append(troughs, x)
            } else {
              
              troughiness <- max((pred[x,2]/pred[y,2]),(pred[x,2]/pred[z,2]))
              
              if (troughiness < 0.8){
                troughs <- append(troughs, x)
              }
            }
          }
        }
        
        
        if (pred[365,2] < pred[364,2]){
          troughs <- append(troughs, 365)
        }
        
        
        # stop everything if there are no legitimate troughs (i.e. the GAM just goes up and up)
        if (length(troughs)==0){
          peak1 <- "Fail"
          true.trough <- "Fail"
          peak2 <- "Fail"
        } else {
          
        
          # now, use a series of ifelse to decide:
          # (1) if there are only one or two peaks, use trough 1
          # (2) if 'trough 2' is actually on day 365 (i.e. the GAM is still declining at the end), use trough 1
          # (3) if there are 3 peaks, choose from trough 1 and trough 2 which is lowest (reasoning that this less likely to be the false trough)
          # (4) if there are more than 3 peaks, just choose the lowest trough
          # this will not give a perfect answer in the messiest cases but should get it right in the vast majority
          
          if (length(troughs) == 1){
            true.trough <- troughs[1]
          } else if (n.peaks <= 2){
            true.trough <- troughs[1]
          } else if (troughs[2] == 365){
            true.trough <- troughs[1]
          } else if (n.peaks == 3){
            true.trough <- ifelse(pred[troughs[1],2] > pred[troughs[2],2], troughs[2], troughs[1])
          } else if (n.peaks > 3){
            trough.preds <- NULL
            
            for (i in 1:length(troughs)-1){
              trough.preds <- append(trough.preds,pred[troughs[i],2])
            }
            
            true.trough <- troughs[which.min(trough.preds)]
          }
          
          
          
          
          # finally, find the peaks of the generations on either side of that trough
          
          # first, find all peaks before the first trough and choose the highest
          peak1s <- peaks[which(peaks < true.trough)]
          
          if (length(peak1s) == 1){
            peak1 <- peak1s[1]
          } else if (length(peak1s) == 2){
            peak1 <- ifelse(pred[peak1s[2],2] < pred[peak1s[1],2], peak1s[1], peak1s[2])
          } else if (length(peak1s) > 2){
            peak1.preds <- NULL
            
            for (i in 1:length(peak1s)){
              peak1.preds <- append(peak1.preds, pred[peak1s[i],2])
            }
            
            peak1 <- peak1s[which.max(peak1.preds)]
          }
          
          # then find the peaks after the first trough and choose the higher of the first two
          # i.e. the peak of the second generation
          # make an adjustment to allow for univoltines (only one peak)
          
          peak2s <- peaks[which(peaks > true.trough)]
          
          
          if (n.peaks==1){
            peak2 <- "Fail"
          } else if (length(peak2s) == 1){
            peak2 <- peak2s[1]
          } else if (length(peak2s) == 2){
            peak2 <- ifelse(pred[peak2s[2],2] < pred[peak2s[1],2], peak2s[1], peak2s[2])
          } else if (length(peak2s) > 2){
            peak2.preds <- NULL
            
            for (i in 1:length(peak2s)){
              peak2.preds <- append(peak2.preds, pred[peak2s[i],2])
            }
            
            peak2 <- peak2s[which.max(peak2.preds)]
          }
          
          
        }
      }
    }
    }
    # append the first peak and the true first trough to the output vector
    peak.day <- append(peak.day,peak1)  
    first.trough <- append(first.trough, true.trough)
    second.peak <- append(second.peak, peak2)
    
    # now calculate the ratio of abundance before and after the first trough
    if (true.trough == "Fail"){
      ratio <- "Fail"
      first.gen <- "Fail"
      second.gen <- "Fail"
    } else if (n.peaks == 1){
      ratio <- 0
      gen1 <- sum(pred$pred)
      gen2 <- "Fail"
    } else if (365 %in% peaks) {
      ratio <- "Fail"
      first.gen <- "Fail"
      second.gen <- "Fail"
    } else {
      
      pre.peak <- pred[which(pred$day < true.trough), ]
      post.peak <- pred[which(pred$day > true.trough), ]
      
      gen1 <- sum(pre.peak$pred)
      gen2 <- sum(post.peak$pred)
      
      ratio <- gen2/gen1
    }
    
    voltinism.ratio <- append(voltinism.ratio,ratio)
    first.gen <- append(first.gen, gen1)
    second.gen <- append(second.gen, gen2)
  
  previous.sp <- current.sp
  
  }
  
  # now stitch this year's vectors together into an output frame
  output <- data.frame(cbind(loop.vec,RECS,peak.day,first.trough,second.peak,voltinism.ratio,first.gen,second.gen))
  
  colnames(output) <- c("siteXspecies","RECS","PEAKDAY","FIRST.TROUGH","SECOND.PEAK","VOLTINISM.RATIO","FIRST.GEN","SECOND.GEN")
  
  # label the year
  output$YEAR <- year
  
  
  # merge this info into the other summary data
  out.summary <- merge(summary.good,output)
  
  

  
  
  
  # and attach it into the overall output
  final.list <- append(final.list,list(out.summary))
  
}

### prepare the data for output

final <- do.call("rbind",final.list)



write.table(final, file = "../Data/RIS/GAMphenology_refined.csv", sep=",", row.names = F)
