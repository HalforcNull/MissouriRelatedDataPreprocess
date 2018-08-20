
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(dplyr)


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
colnames(PC3M) <- c('ensembl_gene_id', 'PC3MNT1', 'PC3MNT2', 'PC3MNT3', 'PC3MTR1', 'PC3MTR2', 'PC3MTR3')


PC3<-rawresult %>%
  select(one_of(c('ensembl_gene_id', 'q13_0', 'q14_0', 'q15_0', 'q16_0', 'q17_0', 'q18_0') ))
colnames(PC3) <- c('ensembl_gene_id', 'PC3_NT1', 'PC3_NT2', 'PC3_NT3', 'PC3_TR1', 'PC3_TR2', 'PC3_TR3')



write.csv(PC3M, file="pc3m.csv", row.names = FALSE)
write.csv(PC3, file="pc3.csv", row.names = FALSE)


