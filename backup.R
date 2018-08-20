library(dplyr)

gene.count <- read.table("Missouri Data Result/genes.count_table", header = TRUE)
gene.attr <- read.table("Missouri Data Result/genes.attr_table", header = TRUE)

gene.tracking.name <- gene.attr %>% select(one_of(c("tracking_id", "gene_short_name")))


result <- dplyr::left_join(gene.count, gene.tracking.name, by="tracking_id")



