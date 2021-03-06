#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")
library("taxonomizr")

# Read *.m8 file with best hit NCBI accession number
if (file.info(snakemake@input[["fhtbl"]])$size == 0){
    file.create(snakemake@output[["linout"]])
    quit()
}
m8 <- read_tsv(snakemake@input[["fhtbl"]], col_names = FALSE)
colnames(m8) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

# Create NCBI accession number vector
acc <- m8$target

# Convert accessions to taxonomic IDs
ids <- accessionToTaxa(acc, snakemake@params[["taxtax"]])

# Get taxonomic lineage
ncbi_tax <- as_tibble(getTaxonomy(ids, snakemake@params[["taxtax"]]))

# Bind m8 file to lineage
seqids <- select(m8, "query")
mmseqs_pviral_nt_lineage <- cbind(seqids, ncbi_tax)

# Write results to table
write_tsv(mmseqs_pviral_nt_lineage, path = snakemake@output[["linout"]])
