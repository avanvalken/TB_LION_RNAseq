library(tidyverse)
library(data.table)
library(stringr)
library(EnsDb.Hsapiens.v86)

# set wd to project loacation 
setwd("../")

# Read in files
x <- list.files(path="data/counts", full.names=T)
z <- lapply(x, read.table, sep="\t", skip=4)
id <- list.files(path="data/counts", full.names=F)
id <- gsub("_ReadsPerGene.out.tab", "", id)
names(z) <- id

# make dataframe of counts
z <- lapply(z, function(x) x[,c(1,4)])
z <- lapply(z, column_to_rownames, var='V1')
z <- mapply(function(x,y) setnames(x, "V4", y), z, names(z), SIMPLIFY = F)
df <- bind_cols(z) 
df <- rownames_to_column(df, var="ensemble_gene_id")

# strip off suffixes
df$ensemble_gene_id <- str_replace(df$ensemble_gene_id,
                        pattern = ".[0-9]+$",
                        replacement = "")


rm(z)

# change ensemble to hgnc_symbol
# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- df$ensemble_gene_id
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
length(unique(geneIDs1$SYMBOL)) #55166

# 2. match gene symbols to ensembl
df.1 <- merge(geneIDs1, df, by.x="GENEID", by.y="ensemble_gene_id")

gns <- genes(EnsDb.Hsapiens.v86)
gns2 <- gns[gns$gene_biotype %in% "protein_coding"]

df.2 <- dplyr::filter(df.1, SYMBOL %in% gns2$gene_name)

# save counts file
write_tsv(df.2, "data/counts.txt")



# library('org.Hs.eg.db')
# 
# genes.1 <- data.frame("ensemble_gene_id"=df$ensemble_gene_id)
# genes.1$symbol <- mapIds(org.Hs.eg.db, keys = genes.1$ensemble_gene_id, keytype = "ENSEMBL", column="SYMBOL")
# 
# sum(is.na(genes.1$symbol)==F) #37064
