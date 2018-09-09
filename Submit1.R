
#================       Report All Data together      =======================


# 

library(dplyr)

tb184A1 <- read.csv('184A1(Kallisto).csv', check.names = FALSE)
colnames(tb184A1)[1] <- 'Gene'
tbC32 <- read.csv('C32(Kallisto).csv', check.names = FALSE)
colnames(tbC32)[1] <- 'Gene'
tbCCD18 <- read.csv('CCD18(Kallisto).csv', check.names = FALSE)
colnames(tbCCD18)[1] <- 'Gene'
tbHEM <- read.csv('HEM(Kallisto).csv', check.names = FALSE)
colnames(tbHEM)[1] <- 'Gene'
tbHT29 <- read.csv('HT29(Kallisto).csv', check.names = FALSE)
colnames(tbHT29)[1] <- 'Gene'
tbhTERT <- read.csv('hTERT(Kallisto).csv', check.names = FALSE)
colnames(tbhTERT)[1] <- 'Gene'
tbMCF7 <- read.csv('MCF7(Kallisto).csv', check.names = FALSE)
colnames(tbMCF7)[1] <- 'Gene'
tbPANC <- read.csv('PANC(Kallisto).csv', check.names = FALSE)
colnames(tbPANC)[1] <- 'Gene'

tbPC3 <- read.csv('PC3(Kallisto).csv', check.names = FALSE)
colnames(tbPC3)[1] <- 'Gene'
tbPC3M <- read.csv('PC3M(Kallisto).csv', check.names = FALSE)
colnames(tbPC3M)[1] <- 'Gene'
tbRW <- read.csv('RW(Kallisto).csv', check.names = FALSE)
colnames(tbRW)[1] <- 'Gene'
tbTRAMP <- read.csv('TRAMP(Kallisto).csv', check.names = FALSE)
colnames(tbTRAMP)[1] <- 'Gene'



AllReadCountTable <- 
  full_join(
    full_join(
      full_join(
        full_join(
          tb184A1, tbC32, by='Gene' 
        ),
        full_join(
          tbCCD18, tbHEM, by='Gene' 
        ),
        by='Gene'
      ),
      full_join(
        full_join(
          tbHT29, tbhTERT, by='Gene' 
        ),
        full_join(
          tbMCF7, tbPANC, by='Gene' 
        ),
        by='Gene'
      ),
      by='Gene'
    ),
    full_join(
      full_join(
        tbPC3, tbPC3M, by='Gene' 
      ),
      full_join(
        tbRW, tbTRAMP, by='Gene' 
      ),
      by='Gene'
    ),
    by='Gene'
  )
  
write.csv(AllReadCountTable, file='Phase_1_Result/AllReadCount.csv')

AllProjectXut <- AllReadCountTable[,1:49]
write.csv(AllProjectXut, file='Phase_1_Result/AllProjectXut.csv')

AllMissouri <- AllReadCountTable[,c(1, 50:73)]
write.csv(AllMissouri, file='Phase_1_Result/AllMissouri.csv')


## COmbine all FC and FDR table
tbDEG184A1 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_184A1.csv', check.names = FALSE)
colnames(tbDEG184A1)[1] <- 'Gene'
colnames(tbDEG184A1)[2] <- 'Symbol_184A1'
tbDEGC32 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_C32.csv', check.names = FALSE)
colnames(tbDEGC32)[1] <- 'Gene'
colnames(tbDEGC32)[2] <- 'Symbol_C32'
tbDEGCCD18 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_CCD18.csv', check.names = FALSE)
colnames(tbDEGCCD18)[1] <- 'Gene'
colnames(tbDEGCCD18)[2] <- 'Symbol_CCD18'
tbDEGHEM <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_HEM.csv', check.names = FALSE)
colnames(tbDEGHEM)[1] <- 'Gene'
colnames(tbDEGHEM)[2] <- 'Symbol_HEM'
tbDEGHT29 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_HT29.csv', check.names = FALSE)
colnames(tbDEGHT29)[1] <- 'Gene'
colnames(tbDEGHT29)[2] <- 'Symbol_HT29'
tbDEGhTERT <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_hTERT.csv', check.names = FALSE)
colnames(tbDEGhTERT)[1] <- 'Gene'
colnames(tbDEGhTERT)[2] <- 'Symbol_hTERT'
tbDEGMCF7 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_MCF7.csv', check.names = FALSE)
colnames(tbDEGMCF7)[1] <- 'Gene'
colnames(tbDEGMCF7)[2] <- 'Symbol_MCF7'
tbDEGPANC <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_PANC.csv', check.names = FALSE)
colnames(tbDEGPANC)[1] <- 'Gene'
colnames(tbDEGPANC)[2] <- 'Symbol_PANC'

tbDEGPC3 <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_PC3.csv', check.names = FALSE)
colnames(tbDEGPC3)[1] <- 'Gene'
colnames(tbDEGPC3)[2] <- 'Symbol_PC3'
tbDEGPC3M <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_PC3M.csv', check.names = FALSE)
colnames(tbDEGPC3M)[1] <- 'Gene'
colnames(tbDEGPC3M)[2] <- 'Symbol_PC3M'
tbDEGRW <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_RW.csv', check.names = FALSE)
colnames(tbDEGRW)[1] <- 'Gene'
colnames(tbDEGRW)[2] <- 'Symbol_RW'
tbDEGTRAMP <- read.csv('Phase_1_Result/FC FDR of each cell line/DEGResult_TRAMP.csv', check.names = FALSE)
colnames(tbDEGTRAMP)[1] <- 'Gene'
colnames(tbDEGTRAMP)[2] <- 'Symbol_TRAMP'


AllDEGResult <- 
  full_join(
    full_join(
      full_join(
        full_join(
          tbDEG184A1, tbDEGC32, by='Gene' 
        ),
        full_join(
          tbDEGCCD18, tbDEGHEM, by='Gene' 
        ),
        by='Gene'
      ),
      full_join(
        full_join(
          tbDEGHT29, tbDEGhTERT, by='Gene' 
        ),
        full_join(
          tbDEGMCF7, tbDEGPANC, by='Gene' 
        ),
        by='Gene'
      ),
      by='Gene'
    ),
    full_join(
      full_join(
        tbDEGPC3, tbDEGPC3M, by='Gene' 
      ),
      full_join(
        tbDEGRW, tbDEGTRAMP, by='Gene' 
      ),
      by='Gene'
    ),
    by='Gene'
  )

write.csv(AllReadCountTable, file='Phase_1_Result/FC FDR of each cell line/AllDEGResult.csv')


