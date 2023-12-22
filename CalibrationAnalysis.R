# Analysis complete stochastic runs
# load data into R
library(readr)

fnames <- c("CalibrationFullStochastic_1.csv",  "CalibrationFullStochastic_2.csv",  
  "CalibrationFullStochastic_3.csv",  "CalibrationFullStochastic_4.csv",  
  "CalibrationFullStochastic_5.csv",  "CalibrationFullStochastic_6.csv",  
  "CalibrationFullStochastic_7.csv",  "CalibrationFullStochastic_8.csv",  
  "CalibrationFullStochastic_9.csv",  "CalibrationFullStochastic_10.csv")

Simulations <- do.call(rbind, lapply(fnames, read.csv, sep=";", header = FALSE))

cnames <- c("t", "S", "Vhist", "Vnew", "E", "I", "R", "cumulI", "I7", "I30", 
            "I60", "I90", "I150", "sim", "N", "propR", "seed", "VEhist",
            "VEnew1", "upperBeta", "upperGamma1", "upperGamma2", "newvaccine", 
            "StartWeekOutbreak")

colnames(Simulations) <- cnames

# Antibodies
Simulations$Ab <- (Simulations$Vhist + Simulations$R + Simulations$Vnew) /
  Simulations$N
  
# End size is between 1200 and 1800
EpiSize <- Simulations[Simulations$cumulI>1200 & Simulations$cumulI<1800,]

# Endtime is between 16*7 = 112 days after first case (minimum) and 202 days
EndTime <- EpiSize[EpiSize$t > 112 & EpiSize$t <202,]

# infections in time
M1 <- EndTime[EndTime$I30 < 150,]
M2 <- M1[M1$I60 < 850,]
M3 <- M2[M2$I90 < 1550,]

# Antibodies between 35% and 45%
Ab <- M3[M3$Ab>0.35 & M3$Ab<0.55,]

# About 18,000 vaccinations given (first vaccine)
Vaccinated <- Ab[Ab$Vnew>14000 & Ab$Vnew<22000,]

# additional infections in time
# output matlab file by day
TimeAnalysis <- Vaccinated

TimeAnalysis$M1 <- TimeAnalysis$`28`
TimeAnalysis$M2 <- TimeAnalysis$`56` - TimeAnalysis$`28`
TimeAnalysis$M3 <- TimeAnalysis$`84` - TimeAnalysis$`56`
TimeAnalysis$M4 <- TimeAnalysis$`112` - TimeAnalysis$`84`
TimeAnalysis$M5 <- TimeAnalysis$`140` - TimeAnalysis$`112`

# number of diagnosed cases after four weeks is 95
FirstMonth <- TimeAnalysis[TimeAnalysis$M1 < 95 * 1.5,]

# number of diagnosed cases in weeks 4 -8 is 454
SecondMonth <- FirstMonth[FirstMonth$M2>0.5*454 & FirstMonth$M2<1.5*454,]

# number of diagnosed cases in weeks 8-12 is 476
ThirdMonth <- SecondMonth[SecondMonth$M3 > 0.5 * 476 & SecondMonth$M3<1.5*476,]

# number of diagnosed cases in weeks 12-16 is 170
FourthMonth <- ThirdMonth[ThirdMonth$M4 <1.5*170,]

# calibrated set is FourthMonth