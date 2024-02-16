cat("\014") #clears console
rm(list=ls()) 
graphics.off()

#### load relevant packages ####
library("readxl")
library("tidyverse")
library("reshape")
library("RColorBrewer")
library("colorspace")

##### data input ####
dat <- read_excel("00_DataInputs/20221107_OD600_Hmari_CDM Exp 03(Responses).xlsx", 
                  col_types = "guess")
expSetup <- read_excel("00_DataInputs/20221107_experimentalSetup.xlsx", 
                       col_types = "guess")
# if timestamp isn't recognized as Date, run the following code
# Data$Date <- as.POSIXct(Data$Timestamp,
#                               origin="1899-12-30",
#                               tz="GMT")

#### inputs
Name <- "Hmari"
Date <- "20221107" # Date of inoculation


# # to test function
Data <- dat
ExperimentDetails <- expSetup
experimentName <- Name
experimentDate <- Date

# #### define function ####
# PlotODs <- function(experimentName, 
#                     experimentDate, 
#                     Data, 
#                     ExperimentDetails){

#### calculate elapsed time ####   
  Data$ElapsedTime <- c("NA")
  
  for(i in 1:nrow(Data)){
    Data$ElapsedTime[i] <- as.numeric(difftime( as.POSIXct( Data$Timestamp[i], format="%Y-%m-%dT%H:%M"), 
                                    as.POSIXct( Data$Timestamp[1], format="%Y-%m-%dT%H:%M"), 
                                    units="hour"))
  }
  
  Data$ElapsedTime <- as.numeric(Data$ElapsedTime)
  
  #### convert data to longform ####
  
  # col names for all experimental flasks
  datacols <- Data %>% 
    dplyr:: select(grep("Rep", names(Data)), grep("Control", names(Data)))
  flasks <- colnames(datacols)
  
  longData = Data %>%
    gather(flasks, key = Flask, value = OD)
  

  #### add "OD2" column with no negative numbers or 0s ####
  longData$OD2 <- longData$OD

  colstochange <- c("OD2")
  for (i in colstochange){
    longData[,i][longData[,i] <= 0] <- 0.001
  }
  
  # remove rows with NA in OD2 column
  longData <- longData[!is.na(longData$OD2),]

  #### remove outliers    #### 
  # only remove from "OD2" column, leave OD column with original data
  # longData[1,"OD2"] <- NA
  # remove erroneous control reasing

  

  # add inoculation date & exp name to exp details sheet
  ExperimentDetails$InoculationDate <- experimentDate
  ExperimentDetails$ExpName <- experimentName
  
  #### merge long Data with experimental set up ####
  longData <- merge(longData, ExperimentDetails)
  
  
  #### export long data ####
  write.csv(longData, file = paste(experimentDate,
                                  "_",
                                  experimentName,
                                    "_long.csv",
                                    sep = ""))
  
  #### create OD summary based on treatment ####
  
  summ <- longData %>%
    group_by(Treatment, Timestamp) %>%
    summarize(ODmean = mean(OD2, na.rm = T),
              ODsd = sd(OD2, na.rm = T),
              ODcount = length(OD2),
              ElapsedTime = mean(ElapsedTime),
              OD2mean = mean(OD2, na.rm = T),
              OD2sd = sd(OD2, na.rm = T),)
  summ$ODse <- summ$ODsd / (sqrt(summ$ODcount))
  summ$OD2se <- summ$OD2sd / (sqrt(summ$ODcount))
  
  #### merge with ExperimentDetails 
  # remove replicate names from metadata
  ExpMetaData <- ExperimentDetails[ -c(1) ] 
  # remove duplicate rows
  ExpMetaData <- distinct(ExpMetaData)
  summ <- merge(summ, ExpMetaData)
  
  
  #### export summary data ####
  write.csv(summ, file = paste(experimentDate,
                                   "_",
                                   experimentName,
                                   "_summary.csv",
                                   sep = ""))
  
  
  #### PLOTS ####
  #### plot all replicates on one figure ####
  cols <- ExperimentDetails$Col
  symbology <- as.data.frame(cbind(flasks, cols))
  colnames(symbology)[1] <- c("Flask")
  symbology$Flask <- as.character(symbology$Flask)
  symbology$colsToUse <- as.character(symbology$cols)

  Controls <- longData[grepl("Control", longData$Flask),]
  listOfControls <- unique(Controls$Flask)

  # use this code to set symbols for replicates vs controls
  # symbology$pch <- 16
  # for (k in 1:nrow(symbology)){
  #   if(symbology$Flask[k] %in% listOfControls){
  #     symbology$pch[k] = 1
  #   }
  # }

  longData <- merge(longData, symbology)
  longData$colsToUse <- as.character(longData$colsToUse)

   
  maxOD <- max(longData$OD2, na.rm = T) * 1.1
  maxTime <- max(longData$ElapsedTime, na.rm = T)
  experimentStart <- as.Date(format(min(longData$Timestamp), "%Y-%m-%d"))
  
  png(paste("01_DataWranglingFigs/",
            experimentDate,
            "_alldata_plot.png",
            sep = ""),
      width = 9, height = 4, units = 'in', res = 300)
  par(mfrow=c(1,2),
      mar=c(4,6,4,0),
      oma = c(1, 1, 1, 10))
  plot(longData$ElapsedTime, longData$OD2,
       las = 1,
       cex = 1.5,
       type = "p",
       pch = longData$pch,
       col = longData$Col,
       bg = "white",
       ylim = c(0, maxOD),
       xlim = c(0, maxTime),
       ylab = "",
       xlab = "",
       cex.axis = 1
       )
  title(main = experimentName, 
        line = -1, 
        cex.main = 1.5,
        outer = TRUE,
        adj = 0)
  title(main = paste("Inoculated:  ", experimentStart, sep = ""), 
        line = -2.5, 
        cex.main = 1.2,
        outer = TRUE,
        adj = 0)
  title(xlab = "Time (Hours)", line = 2.7, cex.lab = 1.2)
  title(ylab = expression("OD"[600]), 
        line = 2.7, 
        cex.lab = 1.2)
  
  # add lines for each experiment
  for (i in 1:length(flasks)){
    temp <- subset(longData, longData$Flask == flasks[i])
    temp <- temp[!is.na(temp$OD2),]
    # sort temp by elapsed time
    temp <- temp[with(temp, order(ElapsedTime)), ]
    lines(temp$ElapsedTime, temp$OD2,
          col = temp$Col,
          lwd = 1.5)
  }
  
  # add points for each experiment (so that points are on top of lines)
  for (i in 1:length(flasks)){
    temp <- subset(longData, longData$Flask == flasks[i])
    temp <- temp[!is.na(temp$OD2),]
    points(temp$ElapsedTime, temp$OD2,
          col = temp$Col,
          bg = "white",
          cex = 1.5,
          pch = temp$pch)
  }
  
  plot(longData$ElapsedTime, longData$OD2,
       las = 1,
       cex = 1.5,
       type = "p",
       pch = longData$pch,
       col = longData$Col,
       bg = "white",
       ylim = c(0.001, maxOD),
       xlim = c(0, maxTime),
       yaxt = "n",
       ylab = "",
       xlab = "",
       cex.axis = 1,
       log = "y")
  
  axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
       las = 1,
       labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
       tck = -0.05,
       cex.axis = 0.9,
       line = 0)
  axis(2, at = c(seq(0.001, 0.009, 0.001),
                 seq(0.01, 0.09, 0.01),
                 seq(0.1, 0.9, 0.1),
                 seq(1, 9, 1)),
       las = 1,
       labels = F,
       tck = -0.025,
       cex.axis = 0.9,
       line = 0)

  title(xlab = "Time (Hours)", line = 2.7, cex.lab = 1.2)
  title(ylab = expression("Log OD"[600]), line = 3, cex.lab = 1.2)
  
  # add lines for each experiment
  for (i in 1:length(flasks)){
    temp <- subset(longData, longData$Flask == flasks[i])
    temp <- temp[!is.na(temp$OD2),]
    # sort temp by elapsed time
    temp <- temp[with(temp, order(ElapsedTime)), ]
    lines(temp$ElapsedTime, temp$OD2,
          col = temp$Col,
          lwd = 1.5)
  }
  
  # add points for each experiment (so that points are on top of lines)
  for (i in 1:length(flasks)){
    temp <- subset(longData, longData$Flask == flasks[i])
    temp <- temp[!is.na(temp$OD2),]
    points(temp$ElapsedTime, temp$OD2,
           col = temp$Col,
           bg = "white",
           cex = 1.5,
           pch = temp$pch)
  }
  
  # make dataframe for easy legend plotting
  legendDat <- subset(longData, longData$ElapsedTime == 0)
  legend("topright",
         xpd = NA, # plots in outer margin area
         cex = 0.6,
         inset = c(-0.8, 0),
         bty = "n",
         legend = legendDat$Flask,
         text.col = "black",
         col = legendDat$Col,
         pch = legendDat$pch,
         pt.cex = 1.2,
         lwd = c(rep(2, length(legendDat$Flask))))
  dev.off()
  
  
  
  
  #### plot mean ODs for each treatment on one plot ####
  maxmeanOD <- max(summ$ODmean, na.rm = T)*1.2
  
  treatments <- unique(summ$Treatment)

  symbology2 <- as.data.frame(cbind(expSetup$Treatment, expSetup$Col))
  colnames(symbology2) <- c("Treatment", "cols2")
  
  
  symbology2 <- distinct(symbology2)
  
  Controls <- longData[grepl("Control", longData$Treatment),]
  listOfControls <- unique(Controls$Treatment)
  
  symbology2$pch <- 16
  for (k in 1:nrow(symbology2)){
    if(symbology2$Treatment[k] %in% listOfControls){
      symbology2$pch[k] = 21
    }
  }
  
  
  summ <- merge(summ, symbology2)
  summ$cols2 <- as.character(summ$cols2)
  
  # change Nan to NA
  colstochange <- c("ODmean") # changes NaN to NA
  for (k in colstochange){
    summ[,k][is.nan(summ$ODmean)] <- NA
  }
  
  summ <- summ[!is.na(summ$ODmean), ]
  
  
  
  png(paste("01_DataWranglingFigs/",
            experimentDate,
            "_meanOD_plot.png",
            sep = ""),
      width = 9, height = 4, units = 'in', res = 300)
  par(mfrow=c(1,2),
      mar=c(4,6,4,0),
      oma = c(1, 1, 1, 10)) 
  plot(summ$ElapsedTime, summ$OD2mean,
       las = 1,
       cex = 1.5,
       type = "p",
       pch = summ$pch,
       col = summ$cols2,
       bg = "white",
       ylim = c(0, maxmeanOD),
       xlim = c(0, maxTime),
       ylab = "",
       xlab = "",
       cex.axis = 1
  )
  title(main = experimentName, 
        line = -1, 
        cex.main = 1.5,
        outer = T,
        adj = 0)
  title(main = paste("Inoculated:  ", experimentStart, sep = ""), 
        line = -2.5, 
        cex.main = 1.2,
        outer = T,
        adj = 0)
  title(xlab = "Time (Hours)", line = 2.7, cex.lab = 1.2)
  title(ylab = expression("OD"[600]), line = 2.7, cex.lab = 1.2)
  
  # add lines for each experiment
  for (i in 1:length(treatments)){
    temp2 <- subset(summ, summ$Treatment == treatments[i])
    temp2 <- temp2[with(temp2, order(ElapsedTime)), ] # sort by elapsed time
    temp2 <- temp2[!is.na(temp2$OD2mean),] # rempve NAs
    lines(temp2$ElapsedTime, temp2$OD2mean,
          col = temp2$cols2,
          lwd = 1.5)
    arrows(x0 = temp2$ElapsedTime, 
           y0= temp2$ODmean - temp2$OD2se, 
           x1 = temp2$ElapsedTime, 
           y1 = temp2$ODmean + temp2$OD2se, 
           code = 1, # no cap of bars
           angle = 90, 
           length = 0,
           lwd = 2, 
           col = temp2$cols2) 
  }
  
  # add points for each experiment (so points are on top of lines)
  for (i in 1:length(treatments)){
    temp2 <- subset(summ, summ$Treatment == treatments[i])
    temp2 <- temp2[!is.na(temp2$OD2mean),]
    points(temp2$ElapsedTime, temp2$OD2mean,
          col = temp2$cols2,
          bg = "white",
          pch = temp2$pch,
          cex = 1.5)
  }
  
  
  symbology2$Treatment <- as.character(symbology2$Treatment)
  symbology2$cols2 <- as.character(symbology2$cols2)
  
  
  #### mean ODs on log scale ####
  plot(summ$ElapsedTime, summ$OD2mean,
       las = 1,
       cex = 1.5,
       type = "p",
       pch = summ$pch,
       col = summ$cols2,
       bg = "white",
       ylim = c(0.001, maxmeanOD),
       xlim = c(0, maxTime),
       ylab = "",
       xlab = "",
       yaxt = "n",
       cex.axis = 1,
       log = "y")
  axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
       las = 1,
       labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
       tck = -0.05,
       cex.axis = 0.9,
       line = 0)
  axis(2, at = c(seq(0.001, 0.009, 0.001),
                 seq(0.01, 0.09, 0.01),
                 seq(0.1, 0.9, 0.1),
                 seq(1, 9, 1)),
       las = 1,
       labels = F,
       tck = -0.025,
       cex.axis = 0.9,
       line = 0)
  
  title(xlab = "Time (Hours)", line = 2.7, cex.lab = 1.2)
  title(ylab = expression("Log OD"[600]), line = 3, cex.lab = 1.2)
  
  # add lines for each experiment
  for (i in 1:length(treatments)){
    temp2 <- subset(summ, summ$Treatment == treatments[i])
    temp2 <- temp2[!is.na(temp2$OD2mean),]
    temp2 <- temp2[with(temp2, order(ElapsedTime)), ] # sort by elapsed time
    lines(temp2$ElapsedTime, temp2$OD2mean,
          col = temp2$cols2,
          lwd = 1.5)
    arrows(x0 = temp2$ElapsedTime, 
           y0= temp2$OD2mean - temp2$OD2se, 
           x1 = temp2$ElapsedTime, 
           y1 = temp2$OD2mean + temp2$OD2se, 
           code = 1, # no cap of bars
           angle = 90, 
           length = 0,
           lwd = 2,
           col = temp2$cols2) 
  }
  
  # add points for each experiment (so points are on top of lines)
  for (i in 1:length(treatments)){
    temp2 <- subset(summ, summ$Treatment == treatments[i])
    temp2 <- temp2[!is.na(temp2$OD2mean),]
    points(temp2$ElapsedTime, temp2$OD2mean,
           col = temp2$cols2,
           bg = "white",
           pch = temp2$pch,
           cex = 1.5)
  }
  
  legend("topright",
         xpd = NA, # plots in outer margin area
         cex = 0.8,
         inset = c(-0.85, 0),
         bty = "n",
         legend = symbology2$Treatment,
         text.col = "black",
         col = symbology2$cols2,
         pt.cex = 1.2,
         pch = symbology2$pch,
         lwd = c(rep(2, length(symbology2$Treatment))))
  dev.off()
  
  
  #### make one plot for each treatment ####

 for (i in 1:length(treatments)){

   # subset longdata to plot indiv replicates
   temp4 <- subset(longData, longData$Treatment ==treatments[i])
   temp4 <- subset(longData, longData$Treatment ==treatments[i]) # remove controls
   
   temp3 <- subset(summ, summ$Treatment == treatments[i])
   # sort by elapsed time
   temp3 <- temp3[with(temp3, order(ElapsedTime)), ]
   
   tempmax <- max(temp3$ODmean) + (max(temp3$ODsd, na.rm = T)*2)
   if(is.na(temp3$ODcount[length(temp3)]) == TRUE){ # set max for n=1 treatments
     tempmax <- max(temp3$ODmean) * 1.2
   }
   if(grepl("Control", treatments[i]) == TRUE){ # set max for controls
     tempmax <- 0.3

   }
   
   
   png(paste("01_DataWranglingFigs/",
             experimentDate,
             "_",
             treatments[i],
             ".png",
             sep = ""),
       width = 8, height = 4.5, units = 'in', res = 300)
   par(mfrow=c(1,2),
       mar=c(4,5,4,1),
       oma = c(1, 1, 1, 2))
   plot(temp3$ElapsedTime, temp3$OD2mean,
        las = 1,
        cex = 1.5,
        type = "p",
        pch = "_",
        col = temp3$cols2,
        bg = "white",
        ylim = c(0, tempmax),
        xlim = c(0, maxTime),
        ylab = "",
        xlab = "",
        cex.axis = 1
   )
   title(main = experimentName,
         line = -0.5,
         cex.main = 1.8,
         adj = 0, # left justify
         outer = TRUE) 
   title(main = paste("Treatment:  ", treatments[i], sep = ""),
         outer = TRUE,
         line = -2,
         cex.main = 1.5,
         adj = 0) # left justify   
   title(main = paste("Inoculated:  ", experimentStart, sep = ""),
         outer = TRUE,
         line = -3.2,
         cex.main = 1.1,
         adj = 0) # left justify

   title(xlab = "Time (Hours)", 
         line = -1, 
         cex.lab = 1.4, 
         outer = TRUE)
   title(ylab = expression("OD"[600]), line = 2.8, cex.lab = 1.4)

   arrows(x0 = temp3$ElapsedTime,
            y0 = temp3$ODmean - temp3$OD2se,
            x1 = temp3$ElapsedTime,
            y1 = temp3$ODmean + temp3$OD2se,
            code = 1, # no cap of bars
            angle = 90,
            length = 0,
          lwd= 2,
            col = temp3$cols) 
   
   points(temp4$ElapsedTime, temp4$OD2, # add points for each replicate
          col = temp4$Col,
          pch = 1,
          bg = "white",
          lwd = 1.5,
          cex = 0.7)
   legend("topleft",
          cex = 0.8,
          bty = "n",
          legend = c("Mean ± SE",
                     "Reps"),
          text.col = "black",
          col = temp4$Col,
          pt.cex = 1.2,
          pch = c("+", "o"))
   
   plot(temp3$ElapsedTime, temp3$OD2mean,
        las = 1,
        cex = 1.5,
        type = "p",
        pch = "__",
        col = temp3$cols2,
        bg = "white",
        ylim = c(0.001, tempmax),
        xlim = c(0, maxTime),
        ylab = "",
        xlab = "",
        yaxt = "n",
        cex.axis = 1.2,
        log = "y")
   axis(2, at = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
        las = 1,
        labels = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
        tck = -0.05,
        cex.axis = 1,
        line = 0)
   axis(2, at = c(seq(0.001, 0.009, 0.001),
                  seq(0.01, 0.09, 0.01),
                  seq(0.1, 0.9, 0.1),
                  seq(1, 9, 1)),
        las = 1,
        labels = F,
        tck = -0.025,
        cex.axis = 0.9,
        line = 0)
   title(ylab = expression("Log OD"[600]), line = 2.8, cex.lab = 1.4)
   
   arrows(x0 = temp3$ElapsedTime,
          y0 = temp3$ODmean - temp3$OD2se,
          x1 = temp3$ElapsedTime,
          y1 = temp3$ODmean + temp3$OD2se,
          code = 1, # no cap of bars
          angle = 90,
          length = 0,
          lwd= 2,
          col = temp3$cols) 

   points(temp4$ElapsedTime, temp4$OD2, # add points for each replicate
          # col = "red",
          col = temp4$Col,
          pch = 1,
          bg = "white",
          lwd = 1.5,
          cex = .7)
   dev.off()
 }
# }

# PlotODs(Name, Date, dat, expSetup)
  
  