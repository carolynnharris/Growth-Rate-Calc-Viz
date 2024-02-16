cat("\014") #clears console
rm(list=ls()) 
graphics.off()

#### load relevant packages ####
library("readxl")
library("tidyverse")
library("reshape")
library("yarrr")


##### data input ####
dat <- read.csv("20221107_Hmari_long.csv", stringsAsFactors = F)
expSetup <- read_excel("00_DataInputs/20221107_experimentalSetup.xlsx", 
                       col_types = "guess")

Name <- "Hmari"
Date <- "20221107" # Date of inoculation


n = 2 # number of points to use
Sim = 0.95 # should be <1
includeNegControls = "No" # yes or no


# #### remove outliers from growth curve data
# dat[X,"OD"] <- NA 

# # for testing function
NumberObs = 2
Similarity = 0.95
experimentName = Name
ExperimentDetails = expSetup
experimentDate <- "2022-07-27" # date of inoculation 
includeNegControls = "No"

#### make dataframe to store results ####
headings <- c("Replicate",
              "GrowthRate",
              "DoublingTime",
              "r",
              "Intercept",
              "StartingObs",
              "NumberObs",
              "SacrificeOD")
heading.count <- length(headings)
GrowthRate_compile <- as.data.frame(matrix(ncol = heading.count, 
                                                   nrow = 0))
colnames(GrowthRate_compile) <- headings 

list.out <- list(GrowthRate_compile = GrowthRate_compile)

#### define function ####
# GrowthRateCalculator = function(data,
#                                 ExperimentDetails,
#                                 Name,
#                                 Date,
#                                 NumberObs,
#                                 Similarity,
#                                 includeNegControls){

  #### growth curve calculator ####
  
  data = dat
  
  if(includeNegControls == "No"){
    data <- data[!grepl("Control", data$Flask),]
  }
  
  
  Flasks <- unique(data$Flask)
  
  for (i in 1:length(Flasks)){
  # i = 15
  
    dat2 <- subset(data, data$Flask == Flasks[i])
    df <- subset(dat2, OD2 > 0) # remove negative OD values & NAs
    df <- df[with(df, order(ElapsedTime)), ] # order by time
    
    #### calculate moving average of n points ####
    experimentName = Name
    experimentDate = Date
    NumberObs = n # number of points to use
    percent_similarity = Sim # should be <1
    listlen <- nrow(df) # number observations
    counter <- 1
    regprodmax = 0
    model <- NULL
    startCount <- 0
    df$lnOD <- log(df$OD2)
    
    ##Go through every set of n points and get slopes and correlations
    while(counter + n <= listlen){
      temp <- df[counter:(counter + n),] # subset for starting and ending points
      model <- lm(lnOD ~ ElapsedTime, data = temp)
      
      slope <- as.numeric(coef(model)[2])
      correlation <- sqrt(summary(model)$r.squared)
      regprod <- slope * correlation
      
      ##only keep the new model if it performs better than the previous model
      if(regprod > regprodmax){
        regprodmax = regprod 
        startCount <- counter # observation to start growth rate calc at
        bestmodel <- model
        bestmodeldat <- temp
      }
      counter = counter + 1
      
    }
    
    #### add points to model to see if model performs better ####
    #This part of the algorithm takes the regression that you calculated before and expands it one point at a time
    #This continues until either the product of correlation and slope remains above a threshold of the max value, or the list ends.
    endCount <- startCount + n 
    rpmax = regprodmax
    ExtendedModel <- NULL
    while(endCount <= listlen){
      temp2 <- df[startCount:endCount,]
      model <- lm(lnOD ~ ElapsedTime, data = temp2)
      slope <- as.numeric(coef(model)[2])
      correlation <- sqrt(summary(model)$r.squared)
      newregprod <- slope * correlation
      numberPoints = length(startCount:endCount)
      ##If the product of the slope and correlation is less than a given percent of the prevoius max then return what was in the previous max 
      if (newregprod > percent_similarity * rpmax){
        ExtendedModel <- model
        ExtendedModeldat <- temp2
        numberPointstoUse = numberPoints
        rpmax = newregprod
      }
      
      # add one more observation to model
      endCount <- endCount +1
    }
    
    
    #### calculate results ####
    Intercept = as.numeric(coef(bestmodel)[1]) # intercept from final model
    if(is.null(ExtendedModel) == F){
      Intercept = as.numeric(coef(ExtendedModel)[1])
    }
    
    Rate = as.numeric(coef(bestmodel)[2]) # growth rate (slope) from final model
    if(is.null(ExtendedModel) == F){
      Rate = as.numeric(coef(ExtendedModel)[2])
    }
    r = as.numeric(summary(bestmodel)$r.squared) # r squared 
    if(is.null(ExtendedModel) == F){
      r = as.numeric(summary(ExtendedModel)$r.squared)
    }
    DoublingTime = log(2)/Rate
    
    SacrificeOD <- df$OD[listlen]
    
    
    #### plot results ####
    
    yMax <- max(df$lnOD, na.rm = T) + 0.2
    yMin <- min(df$lnOD, na.rm = T)
    
    png(paste("02_GrowthRateCalculatorFigs/",
              experimentDate,
              "_",
              Flasks[i],
              "_GrowthRate.png",
              sep = ""),
        width = 6, height = 4, units = 'in', res = 300)
    par(mfrow=c(1,1),
        mar=c(4,4,3,8),
        oma = c(0, 0, 0, 0))  
    plot(lnOD ~ ElapsedTime, data = df, # plot all observations
         las = 1,
         pch = 16,
         col = "grey75",
         cex = 1.6,
         cex.axis = 0.8,
         xlab = "",
         ylab = "",
         ylim = c(yMin, yMax))
    title(main = Flasks[i], 
          line = 1.7, 
          cex.main = 1)
    title(main = experimentName, 
          line = 0.5, 
          cex.main = 0.8)
    title(xlab = "Elapsed Time (Hours)", 
          line = 2.3, 
          cex.lab = 1)
    title(ylab = "Log OD600", 
          line = 2.3, 
          cex.lab = 1)
    points(lnOD ~ ElapsedTime, data = bestmodeldat, # add best model with n+1 points
           col = "cyan3",
           cex = 1.5,
           lwd = 3)
    points(lnOD ~ ElapsedTime, data = ExtendedModeldat, # add final model with more points
           col = "slateblue3",
           cex = 1.9,
           lwd = 2)
    abline(bestmodel, # add best model with n+1 points
           lwd = 2,
           lty = 2,
           col = "cyan3")
    abline(ExtendedModel, # add final model with more points
           lwd = 2,
           col = "slateblue3")
    mtext(paste("Growth Rate = ", 
                round(Rate,3),
                " hr-1",
                sep = ""),
          col = "slateblue3",
          side = 1,
          adj = 1.75,
          line = -6.3,
          cex = 0.8)
    mtext(paste("Doubling = ", 
                round(DoublingTime,2),
                " hr",
                sep = ""),
          col = "slateblue3",
          side = 1,
          adj = 1.51,
          line = -5.3,
          cex = 0.8)
    mtext(paste("Sacrifice OD = ", 
                round(SacrificeOD,2),
                sep = ""),
          col = "slateblue3",
          side = 1,
          adj = 1.5,
          line = -4.3,
          cex = 0.8)
    mtext(paste("r2 = ", 
                round(r,3),
                sep = ""),
          col = "slateblue3",
          side = 1,
          adj = 1.23,
          line = -3.3,
          cex = 0.8)
    legend("bottomright",
           xpd = TRUE,
           cex = 0.8,
           inset = c(-0.36, 0),
           bty = "n",
           legend = c("Final Model", "Starter Model"),
           text.col = c("slateblue3", "cyan3"),
           col = c("slateblue3", "cyan3"),
           lty = c(1,2),
           lwd = c(2,2))
    dev.off()
    
    #### populate summary table with results ####
    
    FlaskSummary <- as.data.frame(cbind(Flasks[i], # replicate
                                        round(Rate, 5), # growth rate
                                        DoublingTime,
                                        round(r,5),
                                        round(Intercept, 5), # intercept
                                        startCount,
                                        numberPointstoUse,
                                        SacrificeOD)) 
    colnames(FlaskSummary) <- c("Replicate",
                                "GrowthRate",
                                "DoublingTime",
                                "r",
                                "Intercept",
                                "StaringObs",
                                "NumberObs",
                                "SacrificeOD")
    
    FlaskSummary$Replicate <- as.character(FlaskSummary$Replicate)
    FlaskSummary$GrowthRate <- as.numeric(as.character(FlaskSummary$GrowthRate))
    FlaskSummary$DoublingTime <- as.numeric(as.character(FlaskSummary$DoublingTime))
    FlaskSummary$r <- as.numeric(as.character(FlaskSummary$r))
    FlaskSummary$Intercept <- as.numeric(as.character(FlaskSummary$Intercept))
    FlaskSummary$StaringObs <- as.numeric(as.character(FlaskSummary$StaringObs))
    FlaskSummary$NumberObs <- as.numeric(as.character(FlaskSummary$NumberObs))
    FlaskSummary$SacrificeOD <- as.numeric(as.character(FlaskSummary$SacrificeOD))
    
    list.out$GrowthRate_compile[nrow(list.out$GrowthRate_compile)+1,] <- FlaskSummary
    
  } # ends loop through flasks

  
  #### summarize Growth Rate Data
  # merge with experimental details file
  names(ExperimentDetails)[names(ExperimentDetails)=="Flask"] <- "Replicate"

  # add inoculation date & exp name to exp details sheet
  ExperimentDetails$InoculationDate <- experimentDate
  ExperimentDetails$ExpName <- experimentName
  
  list.out$GrowthRate_compile <- merge(list.out$GrowthRate_compile, ExperimentDetails)
  #### save summary as .csv ####
  write.csv(list.out$GrowthRate_compile,
            file = paste(Date,
                         "_GrowthRate_Compile.csv",
                         sep = ""))
  
  Mergedat <- merge(list.out$GrowthRate_compile, ExperimentDetails)
  GrowthRate_summary <- Mergedat %>%
    group_by(Treatment) %>%
    summarize(GrowthRatemean = mean(GrowthRate, na.rm = T),
              GrowthRatesd = sd(GrowthRate, na.rm = T),
              DoublingTimemean = mean(DoublingTime, na.rm = T),
              DoublingTimesd = sd(DoublingTime, na.rm = T),
              rmean = mean(r, na.rm = T),
              rsd = sd(r, na.rm = T),
              SacrificeODmean = mean(SacrificeOD, na.rm = T),
              SacrificeODsd = sd(SacrificeOD, na.rm = T),
              count = length(GrowthRate))
  
  GrowthRate_summary <- merge(GrowthRate_summary, ExperimentDetails)
  GrowthRate_summary <- subset(GrowthRate_summary, select=-c(Replicate)) # remove replicates
  GrowthRate_summary <- distinct(GrowthRate_summary) # remove duplicate rows
  
  # write out summary
  write.csv(GrowthRate_summary,
            file = paste(Date,
                         "_GrowthRate_Summary.csv",
                         sep = ""))
  
  
  #### plot summary ####
  
  parameters <- c("GrowthRate", "DoublingTime", "r", "SacrificeOD")
  
  for (j in 1:length(parameters)){
    
    # make color palette for pirate plots
    symbology <- as.data.frame(cbind(expSetup$Treatment, expSetup$Col))
    symbology <- distinct(symbology)
    
    colnames(symbology) <- c("Treatment", "Col")
    # remove controls from symbology table 
    if(includeNegControls == "No"){ 
      symbology <- symbology[!grepl("Control", symbology$Treatment),]
    }
    sort.sym <- with(symbology,  symbology[order(Treatment) , ])
    cols <- as.character(sort.sym$Col)
    
    # cols <- c("black", "deeppink2", "blueviolet", "cyan3")
    
    parameter <- Mergedat[parameters[j]]
    Names <- Mergedat["Treatment"]
    
    plotdata <- as.data.frame(cbind(Names, parameter))
    colnames(plotdata) <- c("Treatment", "Parameter")
    plotdata <- arrange(plotdata, Treatment)
    
    png(paste("02_GrowthRateCalculatorFigs/",
              experimentDate,
              "_",
              parameters[j],
              "_mean_plot.png",
              sep = ""),
        width = 6, height = 4, units = 'in', res = 300)
    par(mfrow=c(1,1),
        mar=c(4,4.5,3,2),
        oma = c(0, 0, 0, 0))  
    pirateplot(formula = Parameter ~ Treatment,
               data = plotdata,
               theme = 0,
               pal = cols,
               main = "",
               ylab = "",
               cex.names = 0.7, # makes x-axis labels smaller
               xlab = "",
               point.cex = 1.3,
               point.o = 1, # point opacity
               avg.line.o = 0.6,
               inf.method = 'ci',
               inf.disp = 'rect',
               inf.b.col = cols, # Inf border col
               inf.b.o = 1, # Inference border
               # bar.f.o = 0.6,
               # bar.f.col = cols) # makes bars colors
               bar.f.col = gray(0.75)) # makes bar light grey
    title(main = parameters[j],
          line = 2,
          cex.main = 1.1)
    title(main = paste(experimentName),
          line = 0.8,
          cex.main = 0.9)
    title(xlab = "Treatment",
          line = 2.2,
          cex.lab = 1.3)
    # add axis label based on parameter
    if(parameters[j] == "GrowthRate"){
      title(ylab = "Growth Rate (hr-1)",
          line = 3,
          cex.lab = 1.3)}
    if(parameters[j] == "DoublingTime"){
      title(ylab = "Doubling Time (hr)",
            line = 3,
            cex.lab = 1.3)}
    if(parameters[j] == "r"){
      title(ylab = "r-squared",
            line = 3,
            cex.lab = 1.3)}
    if(parameters[j] == "SacrificeOD"){
      title(ylab = "Sacrifice OD",
            line = 3,
            cex.lab = 1.3)}
    dev.off()
    
  } # end summary plot loop
  
  
  #### plot growth rate v. sacrifice OD ####
  
  maxOD = max(Mergedat$SacrificeOD)*1.5
  maxRate = max(Mergedat$GrowthRate)*1.2
  legendDat <- Mergedat[,-(1:6)]
  legendDat <- legendDat[,-(4:11)]
  legendDat <- legendDat[!duplicated(legendDat), ]
  
  png(paste("02_GrowthRateCalculatorFigs/",
            experimentDate,
            "_GrowthRateVsSacrificeOD.png",
            sep = ""),
      width = 6, height = 4, units = 'in', res = 300)
  par(mfrow=c(1,1),
      mar=c(4,4.5,3,10),
      oma = c(0, 0, 0, 0))  
  
  plot(Mergedat$SacrificeOD ~ Mergedat$GrowthRate,
       las = 1,
       cex = 1.5,
       xlim = c(0, maxRate),
       ylim = c(0, maxOD),
       type = "p",
       pch = Mergedat$pch,
       col = Mergedat$Col,
       ylab = "Sacrifice OD",
       xlab = "Growth Rate")
  legend("topright",
         xpd = TRUE,
         cex = 0.75,
         inset = c(-0.46, 0),
         bty = "n",
         legend = legendDat$Treatment,
         col = legendDat$Col,
         pch = legendDat$pch)
  dev.off()
  
  #### plot doubling time v. sacrifice OD ####
  
  maxTd = max(Mergedat$DoublingTime)*1.2
  
  png(paste("02_GrowthRateCalculatorFigs/",
            experimentDate,
            "_DoublingTimeVsSacrificeOD.png",
            sep = ""),
      width = 6, height = 4, units = 'in', res = 300)
  par(mfrow=c(1,1),
      mar=c(4,4.5,3,10),
      oma = c(0, 0, 0, 0))  
  
  plot(Mergedat$SacrificeOD ~ Mergedat$DoublingTime,
       las = 1,
       cex = 1.5,
       xlim = c(0, maxTd),
       ylim = c(0, maxOD),
       type = "p",
       pch = Mergedat$pch,
       col = Mergedat$Col,
       ylab = "Sacrifice OD",
       xlab = "Doubling Time (hours)")
  legend("topright",
         xpd = TRUE,
         cex = 0.75,
         inset = c(-0.46, 0),
         bty = "n",
         legend = legendDat$Treatment,
         col = legendDat$Col,
         pch = legendDat$pch)
  dev.off()
  
  
  
  
  return(list.out)
# }  # ends growth rate calculator function


#### run function ####
# list.out <- GrowthRateCalculator(dat, 
#                                  expSetup,
#                                  Name, 
#                                  Date, 
#                                  n, 
#                                  Sim, 
#                                  includeNegControls)






