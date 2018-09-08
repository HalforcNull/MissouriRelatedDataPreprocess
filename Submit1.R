
#================       Report All Data together      =======================


# 

library(dplyr)

tb184A1 <- read.csv('184A1(Kallisto).csv', check.names = FALSE)
tbC32 <- read.csv('C32(Kallisto).csv', check.names = FALSE)
tbCCD18 <- read.csv('CCD18(Kallisto).csv', check.names = FALSE)
tbHEM <- read.csv('HEM(Kallisto).csv', check.names = FALSE)
tbHT29 <- read.csv('HT29(Kallisto).csv', check.names = FALSE)
tbhTERT <- read.csv('hTERT(Kallisto).csv', check.names = FALSE)
tbMCF7 <- read.csv('MCF7(Kallisto).csv', check.names = FALSE)
tbPANC <- read.csv('PANC(Kallisto).csv', check.names = FALSE)

tbPC3 <- read.csv('PC3(Kallisto).csv', check.names = FALSE)
tbPC3M <- read.csv('PC3M(Kallisto).csv', check.names = FALSE)
tbRW <- read.csv('RW(Kallisto).csv', check.names = FALSE)
tbTRAMP <- read.csv('TRAMP(Kallisto).csv', check.names = FALSE)

colnames(tb184A1)

AllReadCountTable <- 
  
  tb184A1


