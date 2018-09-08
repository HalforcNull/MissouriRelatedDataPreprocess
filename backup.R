
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("tximport")
biocLite("EnsDb.Hsapiens.v86")
biocLite("DESeq2")
biocLite("tximport")
biocLite("tximportData")
biocLite("ensembldb")
library(biomaRt)
library(dplyr)
library('DESeq2')
library("tximport")
library("readr")
library("tximportData")
library("ensembldb")
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(dplyr)

# tophat+cufflink solution

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)


gene.count <- read.table("Missouri Data Result/genes.count_table", header = TRUE)
gene.attr <- read.table("Missouri Data Result/genes.attr_table", header = TRUE)

gene.tracking.name <- gene.attr %>% select(one_of(c("tracking_id", "gene_short_name")))


rawresult_1 <- dplyr::left_join(gene.count, gene.tracking.name, by="tracking_id") 
names(rawresult_1)[26] <- 'hgnc_symbol'
rawresult<- dplyr::inner_join(rawresult_1, hgnc_swissprot, by = 'hgnc_symbol')

PC3M<-rawresult %>%
  select(one_of(c('ensembl_gene_id', 'q7_0', 'q8_0', 'q9_0', 'q10_0', 'q11_0', 'q12_0') ))
colnames(PC3M) <- c('ensembl_gene_id', 'PC3M_NT1', 'PC3M_NT2', 'PC3M_NT3', 'PC3M_TR1', 'PC3M_TR2', 'PC3M_TR3')


PC3<-rawresult %>%
  select(one_of(c('ensembl_gene_id', 'q13_0', 'q14_0', 'q15_0', 'q16_0', 'q17_0', 'q18_0') ))
colnames(PC3) <- c('ensembl_gene_id', 'PC3_NT1', 'PC3_NT2', 'PC3_NT3', 'PC3_TR1', 'PC3_TR2', 'PC3_TR3')



write.csv(PC3M, file="pc3m.csv", row.names = FALSE)
write.csv(PC3, file="pc3.csv", row.names = FALSE)



gene.count <- read.table("Yale Data Result/genes.count_table", header = TRUE)
gene.attr <- read.table("Yale Data Result/genes.attr_table", header = TRUE)

gene.tracking.name <- gene.attr %>% select(one_of(c("tracking_id", "gene_short_name")))


rawresult_1 <- dplyr::left_join(gene.count, gene.tracking.name, by="tracking_id") 
names(rawresult_1)[8] <- 'hgnc_symbol'
rawresult<- dplyr::inner_join(rawresult_1, hgnc_swissprot, by = 'hgnc_symbol')

MCF7<-rawresult %>%
  select(one_of(c('ensembl_gene_id', 'q1_0', 'q2_0', 'q3_0', 'q4_0', 'q5_0', 'q6_0') ))
colnames(MCF7) <- c('ensembl_gene_id', 'MCF7_NT1', 'MCF7_NT2', 'MCF7_NT3', 'MCF7_TR1', 'MCF7_TR2', 'MCF7_TR3')

write.csv(MCF7, file="MCF7.csv", row.names = FALSE)


# kallisto solution PC3 Verify
PC3NT1 <- read.table("VerifyPC3/kallistoResult2_NT1/abundance.tsv", header = TRUE)
PC3NT2 <- read.table("VerifyPC3/kallistoResult2_NT2/abundance.tsv", header = TRUE)
PC3NT3 <- read.table("VerifyPC3/kallistoResult2_NT3/abundance.tsv", header = TRUE)
PC3TR1 <- read.table("VerifyPC3/kallistoResult2_TR1/abundance.tsv", header = TRUE)
PC3TR2 <- read.table("VerifyPC3/kallistoResult2_TR2/abundance.tsv", header = TRUE)
PC3TR3 <- read.table("VerifyPC3/kallistoResult2_TR3/abundance.tsv", header = TRUE)


PC3NT1 <- PC3NT1 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3NT1)[2] <- 'PC3NT1'
PC3NT2 <- PC3NT2 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3NT2)[2] <- 'PC3NT2'
PC3NT3 <- PC3NT3 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3NT3)[2] <- 'PC3NT3'
PC3TR1 <- PC3TR1 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3TR1)[2] <- 'PC3TR1'
PC3TR2 <- PC3TR2 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3TR2)[2] <- 'PC3TR2'
PC3TR3 <- PC3TR3 %>%
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(PC3TR3)[2] <- 'PC3TR3'

PC3Result <- full_join(
                full_join(  full_join(PC3NT1,PC3NT2, by='target_id'),
                            full_join(PC3NT3,PC3TR1, by='target_id'),
                            by='target_id'
                          ),
                full_join(  PC3TR2, PC3TR3, by='target_id' ),
                by='target_id')

PC3Result$target_id <- sub("[.].*", "", as.character( PC3Result$target_id) )

length(unique(PC3Result$target_id))

write.csv(PC3Result, file="PC3Verify.csv", row.names = FALSE)


# kallisto solution MCF7 Verify
MCF7NT1 <- read.table("VerifyMCF7/kallistoResult_NT1/abundance.tsv", header = TRUE)
MCF7NT2 <- read.table("VerifyMCF7/kallistoResult_NT2/abundance.tsv", header = TRUE)
MCF7NT3 <- read.table("VerifyMCF7/kallistoResult_NT3/abundance.tsv", header = TRUE)
MCF7TR1 <- read.table("VerifyMCF7/kallistoResult_TR1/abundance.tsv", header = TRUE)
MCF7TR2 <- read.table("VerifyMCF7/kallistoResult_TR2/abundance.tsv", header = TRUE)
MCF7TR3 <- read.table("VerifyMCF7/kallistoResult_TR3/abundance.tsv", header = TRUE)


MCF7NT1 <- MCF7NT1 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7NT1)[2] <- 'MCF7NT1'
MCF7NT2 <- MCF7NT2 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7NT2)[2] <- 'MCF7NT2'
MCF7NT3 <- MCF7NT3 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7NT3)[2] <- 'MCF7NT3'
MCF7TR1 <- MCF7TR1 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7TR1)[2] <- 'MCF7TR1'
MCF7TR2 <- MCF7TR2 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7TR2)[2] <- 'MCF7TR2'
MCF7TR3 <- MCF7TR3 %>% 
  select(one_of('target_id', 'tpm'))  ## NOTE: we are using transcript per million (TPM) here.
names(MCF7TR3)[2] <- 'MCF7TR3'


MCF7Result <- full_join(
  full_join(  full_join(MCF7NT1,MCF7NT2, by='target_id'),
              full_join(MCF7NT3,MCF7TR1, by='target_id'),
              by='target_id'
  ),
  full_join(  MCF7TR2, MCF7TR3, by='target_id' ),
  by='target_id')

MCF7Result$target_id <- sub("[.].*", "", as.character( MCF7Result$target_id) )


write.csv(MCF7Result, file="MCF7Verify.csv", row.names = FALSE)



# DESeq2



MCF7NT1 <- read.table("VerifyMCF7/kallistoResult_NT1/abundance.tsv", header = TRUE)
MCF7NT2 <- read.table("VerifyMCF7/kallistoResult_NT2/abundance.tsv", header = TRUE)
MCF7NT3 <- read.table("VerifyMCF7/kallistoResult_NT3/abundance.tsv", header = TRUE)
MCF7TR1 <- read.table("VerifyMCF7/kallistoResult_TR1/abundance.tsv", header = TRUE)
MCF7TR2 <- read.table("VerifyMCF7/kallistoResult_TR2/abundance.tsv", header = TRUE)
MCF7TR3 <- read.table("VerifyMCF7/kallistoResult_TR3/abundance.tsv", header = TRUE)


# dir <- system.file("extdata", package="tximportData")
# samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# 
# ensembldb::m
# 
# 
# library(EnsDb.Hsapiens.v86)
# transcripts(x, columns = listColumns(x, "tx"),
#             filter = AnnotationFilterList(),
#             return.type = "data.frame")
# 
# 
# 
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, "EnsDb.Hsapiens")
# 
# edb <-ah[['AH60977']]
# txs <- transcripts(edb, return.type = "DataFrame")
# k <- keys(edb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txs, k, 'gene_id', 'tx_name')
# 
# head(txs)


esdb <- EnsDb.Hsapiens.v86
newtxs <- transcripts(esdb, return.type = 'data.frame')
k <- keys(esdb, keytype = "TXNAME")
tx2gene <- dplyr::select(newtxs, one_of(c('tx_name', 'gene_id')))
colnames(tx2gene) <- c('TXNAME', 'GENEID')
# 
# txi.kallisto.tsv <- tximport("VerifyMCF7/kallistoResult_TR3/abundance.tsv", type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

files <- c(
  'VerifyMCF7/kallistoResult_NT1/abundance.tsv',
  'VerifyMCF7/kallistoResult_NT2/abundance.tsv',
  'VerifyMCF7/kallistoResult_NT3/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR1/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR2/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR3/abundance.tsv')
names(files) <- c('MCF7NT1','MCF7NT2','MCF7NT3','MCF7TR1','MCF7TR2','MCF7TR3')


txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="MCF7(Kallisto).csv", row.names = TRUE)


files <-c(
  'VerifyPC3/kallistoResult2_NT1/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT2/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT3/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR1/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR2/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR3/abundance.tsv')
names(files) <- c('PC3NT1','PC3NT2','PC3NT3','PC3TR1','PC3TR2','PC3TR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="PC3(Kallisto).csv", row.names = TRUE)



files <-c('kallistoPC3M_NT1/abundance.tsv',
'kallistoPC3M_NT2/abundance.tsv',
'kallistoPC3M_NT3/abundance.tsv',
'kallistoPC3M_TR1/abundance.tsv',
'kallistoPC3M_TR2/abundance.tsv',
'kallistoPC3M_TR3/abundance.tsv')
names(files) <- c('PC3MNT1','PC3MNT2','PC3MNT3','PC3MTR1','PC3MTR2','PC3MTR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="PC3M(Kallisto).csv", row.names = TRUE)


files <-c('kallistoRW_NT1/abundance.tsv',
          'kallistoRW_NT2/abundance.tsv',
          'kallistoRW_NT3/abundance.tsv',
          'kallistoRW_TR1/abundance.tsv',
          'kallistoRW_TR2/abundance.tsv',
          'kallistoRW_TR3/abundance.tsv')
names(files) <- c('RWNT1','RWNT2','RWNT3','RWTR1','RWTR2','RWTR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="RW(Kallisto).csv", row.names = TRUE)


files <-c('kallistoTRAMP_NT1/abundance.tsv',
          'kallistoTRAMP_NT2/abundance.tsv',
          'kallistoTRAMP_NT3/abundance.tsv',
          'kallistoTRAMP_TR1/abundance.tsv',
          'kallistoTRAMP_TR2/abundance.tsv',
          'kallistoTRAMP_TR3/abundance.tsv')
names(files) <- c('TRAMPNT1','TRAMPNT2','TRAMPNT3','TRAMPTR1','TRAMPTR2','TRAMPTR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="TRAMP(Kallisto).csv", row.names = TRUE)

files <-c('kallistoPC3M_NT1/abundance.tsv',
          'kallistoPC3M_NT2/abundance.tsv',
          'kallistoPC3M_NT3/abundance.tsv',
          'kallistoPC3M_TR1/abundance.tsv',
          'kallistoPC3M_TR2/abundance.tsv',
          'kallistoPC3M_TR3/abundance.tsv')
names(files) <- c('PC3MNT1','PC3MNT2','PC3MNT3','PC3MTR1','PC3MTR2','PC3MTR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="PC3M(Kallisto).csv", row.names = TRUE)



# hTERT 13-18
SampleName = 'hTERT'
files <- NULL
for(i in 13:18){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)


# PANC 19-24
SampleName = 'PANC'
files <- NULL
for(i in 19:24){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)

# HT29 25-30
SampleName = 'HT29'
files <- NULL
for(i in 25:30){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)

# CCD18 31-36
SampleName = 'CCD18'
files <- NULL
for(i in 31:36){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)


# HEM 37-42
SampleName = 'HEM'
files <- NULL
for(i in 37:42){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)

# C32 43-48
SampleName = 'C32'
files <- NULL
for(i in 43:48){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)

# 184A1 1-6
SampleName = '184A1'
files <- NULL
for(i in 1:6){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
}
names(files) <- c(paste0(SampleName, 'NT1'), paste0(SampleName, 'NT2'), paste0(SampleName, 'NT3'),
                  paste0(SampleName, 'TR1'), paste0(SampleName, 'TR2'), paste0(SampleName, 'TR3'))
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)




#################==============             targets vs. non-targets               ===================
###       This would mean finding DEG/pathways that are co-deregulated in Tramp, PC3, PC3M, and MCF cells but not in the rest of the samples

######  Missouri Targeting Sample (None targeting sample only contains RW)

###### PC3, PC3M, Tramp vs RW (PC3S3 are excluded)
###### PC3M, Tramp vs RW 

files <-c(
  'VerifyPC3/kallistoResult2_NT1/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT2/abundance.tsv',
  'kallistoPC3M_NT1/abundance.tsv',
  'kallistoPC3M_NT2/abundance.tsv',
  'kallistoPC3M_NT3/abundance.tsv',
  'kallistoTRAMP_NT1/abundance.tsv',
  'kallistoTRAMP_NT2/abundance.tsv',
  'kallistoTRAMP_NT3/abundance.tsv',
  
  
  'VerifyPC3/kallistoResult2_TR1/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR2/abundance.tsv',
  'kallistoPC3M_TR1/abundance.tsv',
  'kallistoPC3M_TR2/abundance.tsv',
  'kallistoPC3M_TR3/abundance.tsv',
  'kallistoTRAMP_TR1/abundance.tsv',
  'kallistoTRAMP_TR2/abundance.tsv',
  'kallistoTRAMP_TR3/abundance.tsv')
names(files) <- c('TPP_NT1','TPP_NT2','TPP_NT3','TPP_NT4','TPP_NT5','TPP_NT6','TPP_NT7','TPP_NT8',
                  'TPP_TR1','TPP_TR2','TPP_TR3','TPP_TR4','TPP_TR5','TPP_TR6','TPP_TR7','TPP_TR8')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="MissouriTarget(Kallisto).csv", row.names = TRUE)


######  Project_Xut None Targeting Sample (Targeting sample only contains MCF7)
SampleName = 'ProjectXut_None_Targeting'
files <- NULL
fnames <- NULL
ntIndex <- 1
trIndex <- 1
for(i in 1:6){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
  if((i-1)%%6<3){
    fnames <- c(fnames, paste0(SampleName, 'NT', as.character(ntIndex)))
    ntIndex=ntIndex+1
  }else{
    fnames <- c(fnames, paste0(SampleName, 'TR', as.character(trIndex)))
    trIndex=trIndex+1
  }
}
for(i in 13:48){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
  if((i-1)%%6<3){
    fnames <- c(fnames, paste0(SampleName, 'NT', as.character(ntIndex)))
    ntIndex=ntIndex+1
  }else{
    fnames <- c(fnames, paste0(SampleName, 'TR', as.character(trIndex)))
    trIndex=trIndex+1
  }
}
names(files) <- fnames
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)


########### ALL Targeting
SampleName = 'All_Targeting'
files <- NULL
fnames <- NULL
ntIndex <- 1
trIndex <- 1


##MCF
files <- c(
  'VerifyMCF7/kallistoResult_NT1/abundance.tsv',
  'VerifyMCF7/kallistoResult_NT2/abundance.tsv',
  'VerifyMCF7/kallistoResult_NT3/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR1/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR2/abundance.tsv',
  'VerifyMCF7/kallistoResult_TR3/abundance.tsv')
## PC3M
files <-c(files,
          'kallistoPC3M_NT1/abundance.tsv',
          'kallistoPC3M_NT2/abundance.tsv',
          'kallistoPC3M_NT3/abundance.tsv',
          'kallistoPC3M_TR1/abundance.tsv',
          'kallistoPC3M_TR2/abundance.tsv',
          'kallistoPC3M_TR3/abundance.tsv')
## TRAMP
files <-c(files,
          'kallistoTRAMP_NT1/abundance.tsv',
          'kallistoTRAMP_NT2/abundance.tsv',
          'kallistoTRAMP_NT3/abundance.tsv',
          'kallistoTRAMP_TR1/abundance.tsv',
          'kallistoTRAMP_TR2/abundance.tsv',
          'kallistoTRAMP_TR3/abundance.tsv')
## PC3 (only s1 and s2)
files <-c(files,
  'VerifyPC3/kallistoResult2_NT1/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT2/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR1/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR2/abundance.tsv')


for(i in 1:18){
  if((i-1)%%6<3){
    fnames <- c(fnames, paste0(SampleName, 'NT', as.character(ntIndex)))
    ntIndex=ntIndex+1
  }else{
    fnames <- c(fnames, paste0(SampleName, 'TR', as.character(trIndex)))
    trIndex=trIndex+1
  }
}
fnames <- c(fnames, paste0(SampleName, 'NT10'), paste0(SampleName, 'NT11'),
            paste0(SampleName, 'TR10'), paste0(SampleName, 'TR11'))

names(files) <- fnames
txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)


#### All None Targeting
SampleName = 'All_Non_Targeting'
files <- NULL
fnames <- NULL
ntIndex <- 1
trIndex <- 1

for(i in 1:6){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
  if((i-1)%%6<3){
    fnames <- c(fnames, paste0(SampleName, 'NT', as.character(ntIndex)))
    ntIndex=ntIndex+1
  }else{
    fnames <- c(fnames, paste0(SampleName, 'TR', as.character(trIndex)))
    trIndex=trIndex+1
  }
}
for(i in 13:48){
  files <- c(files, paste0('kallistoResult_', as.character(i), '/abundance.tsv'))
  if((i-1)%%6<3){
    fnames <- c(fnames, paste0(SampleName, 'NT', as.character(ntIndex)))
    ntIndex=ntIndex+1
  }else{
    fnames <- c(fnames, paste0(SampleName, 'TR', as.character(trIndex)))
    trIndex=trIndex+1
  }
}


files <-c(files,
          'kallistoRW_NT2/abundance.tsv',
          'kallistoRW_NT3/abundance.tsv',
          'kallistoRW_TR2/abundance.tsv',
          'kallistoRW_TR3/abundance.tsv')

fnames <- c(fnames,paste0(SampleName, 'NT', as.character(ntIndex+1)), paste0(SampleName, 'NT',as.character(ntIndex+2)),
            paste0(SampleName, 'TR',as.character(ntIndex+1)), paste0(SampleName, 'TR', as.character(ntIndex+2)))

names(files) <- fnames

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file=paste0(SampleName, "(Kallisto).csv"), row.names = TRUE)




#names(files) <- c('PC3NT1','PC3NT2','PC3TR1','PC3TR2','PC3MNT1','PC3MNT2','PC3MNT3','PC3MTR1','PC3MTR2','PC3MTR3',
#                  'TRAMPNT1','TRAMPNT2','TRAMPNT3','TRAMPTR1','TRAMPTR2','TRAMPTR3')

{
upListInTarget <- c("ENSG00000100906",
                    "ENSG00000144802",
                    "ENSG00000163874",
                    "ENSG00000169242",
                    "ENSG00000115963",
                    "ENSG00000118503",
                    "ENSG00000111859",
                    "ENSG00000145632",
                    "ENSG00000112096",
                    "ENSG00000104312",
                    "ENSG00000111912",
                    "ENSG00000177606",
                    "ENSG00000163739",
                    "ENSG00000171223",
                    "ENSG00000132510",
                    "ENSG00000067082",
                    "ENSG00000125347",
                    "ENSG00000221869",
                    "ENSG00000107968",
                    "ENSG00000109320",
                    "ENSG00000123358",
                    "ENSG00000185215",
                    "ENSG00000128016",
                    "ENSG00000077150",
                    "ENSG00000185022",
                    "ENSG00000120129",
                    "ENSG00000157557",
                    "ENSG00000081041",
                    "ENSG00000162772",
                    "ENSG00000111266",
                    "ENSG00000280962",
                    "ENSG00000095951",
                    "ENSG00000148339",
                    "ENSG00000010818",
                    "ENSG00000163659",
                    "ENSG00000104856",
                    "ENSG00000124391",
                    "ENSG00000137193",
                    "ENSG00000163545",
                    "ENSG00000144655",
                    "ENSG00000167604",
                    "ENSG00000163734",
                    "ENSG00000171522",
                    "ENSG00000151014",
                    "ENSG00000136244",
                    "ENSG00000115009",
                    "ENSG00000166920",
                    "ENSG00000175505",
                    "ENSG00000166016",
                    "ENSG00000173846",
                    "ENSG00000198535",
                    "ENSG00000108551",
                    "ENSG00000135604",
                    "ENSG00000140465",
                    "ENSG00000165474",
                    "ENSG00000131459",
                    "ENSG00000132003",
                    "ENSG00000164949",
                    "ENSG00000119508",
                    "ENSG00000187479",
                    "ENSG00000078401",
                    "ENSG00000204569",
                    "ENSG00000227804",
                    "ENSG00000230995",
                    "ENSG00000231737",
                    "ENSG00000235291",
                    "ENSG00000238104",
                    "ENSG00000206489",
                    "ENSG00000164400",
                    "ENSG00000102962",
                    "ENSG00000006210",
                    "ENSG00000125740",
                    "ENSG00000204490",
                    "ENSG00000206439",
                    "ENSG00000223952",
                    "ENSG00000228321",
                    "ENSG00000228849",
                    "ENSG00000230108",
                    "ENSG00000232810",
                    "ENSG00000179674",
                    "ENSG00000169245",
                    "ENSG00000128271",
                    "ENSG00000267607",
                    "ENSG00000241253",
                    "ENSG00000127528",
                    "ENSG00000007908",
                    "ENSG00000268812",
                    "ENSG00000108342",
                    "ENSG00000167874",
                    "ENSG00000186994",
                    "ENSG00000261618",
                    "ENSG00000176907",
                    "ENSG00000242335",
                    "ENSG00000243570",
                    "ENSG00000241534",
                    "ENSG00000239754",
                    "ENSG00000089692")
upListInNonTarget <- c("ENSG00000124145",
                       "ENSG00000135821",
                       "ENSG00000112096",
                       "ENSG00000115008",
                       "ENSG00000169429",
                       "ENSG00000143878",
                       "ENSG00000124882",
                       "ENSG00000185215",
                       "ENSG00000163739",
                       "ENSG00000115758",
                       "ENSG00000118503",
                       "ENSG00000186480",
                       "ENSG00000185650",
                       "ENSG00000128342",
                       "ENSG00000134954",
                       "ENSG00000100906",
                       "ENSG00000064651",
                       "ENSG00000023445",
                       "ENSG00000177606",
                       "ENSG00000172403",
                       "ENSG00000090339",
                       "ENSG00000067082",
                       "ENSG00000073756",
                       "ENSG00000065534",
                       "ENSG00000168906",
                       "ENSG00000149428",
                       "ENSG00000280682",
                       "ENSG00000115009",
                       "ENSG00000120738",
                       "ENSG00000175592",
                       "ENSG00000149289",
                       "ENSG00000110047",
                       "ENSG00000087074",
                       "ENSG00000117152",
                       "ENSG00000142627",
                       "ENSG00000123358",
                       "ENSG00000120129",
                       "ENSG00000125347",
                       "ENSG00000163661",
                       "ENSG00000128016",
                       "ENSG00000171223",
                       "ENSG00000144802",
                       "ENSG00000143322",
                       "ENSG00000144136",
                       "ENSG00000112972",
                       "ENSG00000265190",
                       "ENSG00000163734",
                       "ENSG00000163874",
                       "ENSG00000182326",
                       "ENSG00000119508",
                       "ENSG00000170345",
                       "ENSG00000197632",
                       "ENSG00000169242",
                       "ENSG00000136244",
                       "ENSG00000160888",
                       "ENSG00000081041",
                       "ENSG00000146278",
                       "ENSG00000164161",
                       "ENSG00000162772",
                       "ENSG00000172216",
                       "ENSG00000221869",
                       "ENSG00000187678",
                       "ENSG00000227231",
                       "ENSG00000230128",
                       "ENSG00000235030",
                       "ENSG00000237155",
                       "ENSG00000196352",
                       "ENSG00000134107",
                       "ENSG00000173391",
                       "ENSG00000184164",
                       "ENSG00000185022",
                       "ENSG00000162413",
                       "ENSG00000115221",
                       "ENSG00000163659",
                       "ENSG00000112715",
                       "ENSG00000159388",
                       "ENSG00000175602",
                       "ENSG00000109320",
                       "ENSG00000163347",
                       "ENSG00000116717",
                       "ENSG00000136158",
                       "ENSG00000132510",
                       "ENSG00000148841",
                       "ENSG00000125538",
                       "ENSG00000159216",
                       "ENSG00000155090",
                       "ENSG00000083799",
                       "ENSG00000127528",
                       "ENSG00000069020",
                       "ENSG00000113070",
                       "ENSG00000095752",
                       "ENSG00000169991",
                       "ENSG00000179388",
                       "ENSG00000211455",
                       "ENSG00000187193",
                       "ENSG00000126003",
                       "ENSG00000166401",
                       "ENSG00000116514",
                       "ENSG00000112658",
                       "ENSG00000143067",
                       "ENSG00000198142",
                       "ENSG00000145779",
                       "ENSG00000170385",
                       "ENSG00000077150",
                       "ENSG00000164949",
                       "ENSG00000136826",
                       "ENSG00000130522",
                       "ENSG00000148339",
                       "ENSG00000011422",
                       "ENSG00000088826",
                       "ENSG00000198355",
                       "ENSG00000119862",
                       "ENSG00000138166",
                       "ENSG00000108342",
                       "ENSG00000198517",
                       "ENSG00000113742",
                       "ENSG00000141682",
                       "ENSG00000135636",
                       "ENSG00000113916",
                       "ENSG00000170006",
                       "ENSG00000134070",
                       "ENSG00000125740",
                       "ENSG00000167695",
                       "ENSG00000167034",
                       "ENSG00000168994",
                       "ENSG00000165891",
                       "ENSG00000177283",
                       "ENSG00000230439",
                       "ENSG00000095951",
                       "ENSG00000069399",
                       "ENSG00000151014",
                       "ENSG00000112149",
                       "ENSG00000162924",
                       "ENSG00000120217",
                       "ENSG00000144655",
                       "ENSG00000111912",
                       "ENSG00000188522",
                       "ENSG00000131979",
                       "ENSG00000105327",
                       "ENSG00000163545",
                       "ENSG00000154914",
                       "ENSG00000153234",
                       "ENSG00000163660",
                       "ENSG00000145244",
                       "ENSG00000162344",
                       "ENSG00000188211",
                       "ENSG00000173918",
                       "ENSG00000173846",
                       "ENSG00000175505",
                       "ENSG00000197019",
                       "ENSG00000175197",
                       "ENSG00000127666",
                       "ENSG00000013441",
                       "ENSG00000111266",
                       "ENSG00000280962",
                       "ENSG00000143333",
                       "ENSG00000146232",
                       "ENSG00000273793",
                       "ENSG00000145365",
                       "ENSG00000196449",
                       "ENSG00000101665",
                       "ENSG00000104856",
                       "ENSG00000115844",
                       "ENSG00000213859",
                       "ENSG00000056558",
                       "ENSG00000163376",
                       "ENSG00000146592",
                       "ENSG00000183092",
                       "ENSG00000249673",
                       "ENSG00000186594",
                       "ENSG00000282800",
                       "ENSG00000277117",
                       "ENSG00000054967",
                       "ENSG00000179859",
                       "ENSG00000107968",
                       "ENSG00000166592",
                       "ENSG00000178607",
                       "ENSG00000006459",
                       "ENSG00000057657",
                       "ENSG00000174945",
                       "ENSG00000160326",
                       "ENSG00000281165",
                       "ENSG00000173110",
                       "ENSG00000109906",
                       "ENSG00000267520",
                       "ENSG00000111012",
                       "ENSG00000173451",
                       "ENSG00000267519",
                       "ENSG00000167604",
                       "ENSG00000226380",
                       "ENSG00000145911",
                       "ENSG00000122877",
                       "ENSG00000270225",
                       "ENSG00000073737",
                       "ENSG00000158050",
                       "ENSG00000100079",
                       "ENSG00000223799",
                       "ENSG00000171453",
                       "ENSG00000125657",
                       "ENSG00000128594",
                       "ENSG00000168334",
                       "ENSG00000123572",
                       "ENSG00000135604",
                       "ENSG00000260196",
                       "ENSG00000269028",
                       "ENSG00000130487",
                       "ENSG00000267607",
                       "ENSG00000275410",
                       "ENSG00000164400",
                       "ENSG00000256671",
                       "ENSG00000172602",
                       "ENSG00000162896",
                       "ENSG00000135625",
                       "ENSG00000124875",
                       "ENSG00000262902",
                       "ENSG00000155530",
                       "ENSG00000141668",
                       "ENSG00000281931",
                       "ENSG00000272114",
                       "ENSG00000198574",
                       "ENSG00000172738",
                       "ENSG00000274979",
                       "ENSG00000123977",
                       "ENSG00000170075",
                       "ENSG00000261618",
                       "ENSG00000099960",
                       "ENSG00000251893",
                       "ENSG00000259884",
                       "ENSG00000049249",
                       "ENSG00000182393",
                       "ENSG00000268812",
                       "ENSG00000244230",
                       "ENSG00000089692",
                       "ENSG00000204490",
                       "ENSG00000206439",
                       "ENSG00000223952",
                       "ENSG00000228321",
                       "ENSG00000228849",
                       "ENSG00000230108",
                       "ENSG00000232810",
                       "ENSG00000281207",
                       "ENSG00000254842",
                       "ENSG00000255443",
                       "ENSG00000008516",
                       "ENSG00000206342",
                       "ENSG00000261604",
                       "ENSG00000234431",
                       "ENSG00000185291",
                       "ENSG00000221949",
                       "ENSG00000163735",
                       "ENSG00000167874",
                       "ENSG00000251279",
                       "ENSG00000120337",
                       "ENSG00000234667",
                       "ENSG00000260464",
                       "ENSG00000244731",
                       "ENSG00000007908",
                       "ENSG00000255521",
                       "ENSG00000258602",
                       "ENSG00000255282",
                       "ENSG00000226738",
                       "ENSG00000159261",
                       "ENSG00000223935",
                       "ENSG00000276462",
                       "ENSG00000118094",
                       "ENSG00000236581",
                       "ENSG00000227507",
                       "ENSG00000156427",
                       "ENSG00000260360",
                       "ENSG00000184545",
                       "ENSG00000278165",
                       "ENSG00000183709",
                       "ENSG00000213344",
                       "ENSG00000200534",
                       "ENSG00000157765",
                       "ENSG00000168703",
                       "ENSG00000283213",
                       "ENSG00000269951",
                       "ENSG00000261114",
                       "ENSG00000237892",
                       "ENSG00000232043",
                       "ENSG00000254746",
                       "ENSG00000275894",
                       "ENSG00000232618",
                       "ENSG00000268568",
                       "ENSG00000207808",
                       "ENSG00000277895",
                       "ENSG00000197110",
                       "ENSG00000261026",
                       "ENSG00000102962",
                       "ENSG00000274976",
                       "ENSG00000221539",
                       "ENSG00000271784",
                       "ENSG00000168928",
                       "ENSG00000230176",
                       "ENSG00000266651",
                       "ENSG00000136286",
                       "ENSG00000206754",
                       "ENSG00000236782",
                       "ENSG00000208028",
                       "ENSG00000221500",
                       "ENSG00000207980",
                       "ENSG00000276216",
                       "ENSG00000263624",
                       "ENSG00000235989",
                       "ENSG00000227908",
                       "ENSG00000265112",
                       "ENSG00000234789",
                       "ENSG00000274322",
                       "ENSG00000200879",
                       "ENSG00000265452",
                       "ENSG00000232530",
                       "ENSG00000261630",
                       "ENSG00000268366",
                       "ENSG00000176907",
                       "ENSG00000244036",
                       "ENSG00000232303",
                       "ENSG00000265096",
                       "ENSG00000200913",
                       "ENSG00000282977",
                       "ENSG00000196664")

DownListInTarget <- c("ENSG00000265972",
                      "ENSG00000125968",
                      "ENSG00000106031",
                      "ENSG00000196659",
                      "ENSG00000180884",
                      "ENSG00000187626",
                      "ENSG00000115507",
                      "ENSG00000251495")
DownListInNonTarget <- c("ENSG00000136111",
                         "ENSG00000168481",
                         "ENSG00000172005",
                         "ENSG00000112984",
                         "ENSG00000166501",
                         "ENSG00000112742",
                         "ENSG00000135914",
                         "ENSG00000134215",
                         "ENSG00000167183",
                         "ENSG00000086548",
                         "ENSG00000197852",
                         "ENSG00000139973",
                         "ENSG00000139263",
                         "ENSG00000134222",
                         "ENSG00000161888",
                         "ENSG00000181274",
                         "ENSG00000066382",
                         "ENSG00000091879",
                         "ENSG00000089101",
                         "ENSG00000168517",
                         "ENSG00000167800",
                         "ENSG00000095932",
                         "ENSG00000165879",
                         "ENSG00000224981",
                         "ENSG00000118946",
                         "ENSG00000133256",
                         "ENSG00000136866",
                         "ENSG00000198416",
                         "ENSG00000161653",
                         "ENSG00000143365",
                         "ENSG00000159905",
                         "ENSG00000123500",
                         "ENSG00000005981",
                         "ENSG00000137397",
                         "ENSG00000186777",
                         "ENSG00000259488",
                         "ENSG00000171116",
                         "ENSG00000197054",
                         "ENSG00000123405",
                         "ENSG00000251246",
                         "ENSG00000267500",
                         "ENSG00000131400",
                         "ENSG00000171847",
                         "ENSG00000255337",
                         "ENSG00000147059",
                         "ENSG00000213057",
                         "ENSG00000277806",
                         "ENSG00000136457",
                         "ENSG00000278206",
                         "ENSG00000266371",
                         "ENSG00000125804",
                         "ENSG00000261360",
                         "ENSG00000225920")
}

upGeneOnlyInTarget <- upListInTarget[!upListInTarget%in%upListInNonTarget]
downGeneOnlyInTarget <- DownListInTarget[!DownListInTarget%in%DownListInNonTarget]

upGeneOnlyInNonTarget <- upListInNonTarget[!upListInNonTarget%in%upListInTarget]
downGeneOnlyInNonTarget <- DownListInNonTarget[!DownListInNonTarget%in%DownListInTarget]



writeClipboard(paste(upGeneOnlyInTarget, sep = '\t'))
writeClipboard(paste(downGeneOnlyInTarget, sep = '\t'))


#================       Report All Data together      =======================


# 
tb184A1 <- read.csv('184A1(Kallisto).csv')
tbC32 <- read.csv('C32(Kallisto).csv')
tbCCD18 <- read.csv('CCD18(Kallisto).csv')
tbHEM <- read.csv('HEM(Kallisto).csv')
tbHT29 <- read.csv('HT29(Kallisto).csv')
tbhTERT <- read.csv('hTERT(Kallisto).csv')
tbMCF7 <- read.csv('MCF7(Kallisto).csv')
PANC <- read.csv('PANC(Kallisto).csv')

PC3 <- read.csv('PC3(Kallisto).csv')
PC3M <- read.csv('PC3M(Kallisto).csv')
RW <- read.csv('RW(Kallisto).csv')
TRAMP <- read.csv('TRAMP(Kallisto).csv')

































################################################################################################
#########################           Paper reading         ######################################
################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite("debrowser")
library('debrowser')



runDE()


