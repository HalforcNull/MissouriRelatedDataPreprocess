
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
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
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("tximport")
biocLite("tximportData")
biocLite("ensembldb")
library('DESeq2')
library("tximport")
library("readr")
library("tximportData")
library("ensembldb")


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
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(dplyr)
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
write.csv(txi.kallisto.tsv$counts, file="MCF7Verify_2.csv", row.names = TRUE)


files <-c(
  'VerifyPC3/kallistoResult2_NT1/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT2/abundance.tsv',
  'VerifyPC3/kallistoResult2_NT3/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR1/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR2/abundance.tsv',
  'VerifyPC3/kallistoResult2_TR3/abundance.tsv')
names(files) <- c('PC3NT1','PC3NT2','PC3NT3','PC3TR1','PC3TR2','PC3TR3')

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion=TRUE )
write.csv(txi.kallisto.tsv$counts, file="PC3Verify_2.csv", row.names = TRUE)




################################################################################################
#########################           Paper reading         ######################################
################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite("debrowser")
library('debrowser')



runDE()


