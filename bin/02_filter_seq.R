#!/usr/bin/env Rscript

# Remove sequences from downstream analyses if they have more then 10% ambiguous bases

library(ape)
library(magrittr)
library(readxl)
library(parallel)
library(dplyr)

args <- commandArgs(TRUE)

rdata <- args[1]
metadata <- args[2]
sequences <- args[3]
outgroup_fastas <- args[4]
args1 <- args

load(rdata)

# read in GISAID fastas and metadata
fastas <- read.dna(sequences, "fasta")
fasta_ids <- names(fastas)


# read in outgroup fasta
outgroup_id <- "MG772933.1"
outgroup_fasta <- read.dna(outgroup_fastas, "fasta")

# filter sequences with >10% ambiguous bases
base_freq <- lapply(1:length(fastas), function (i) base.freq(fastas[i], all=TRUE))
ambiguous_prop <- sapply(base_freq, function (x) sum(x[!(names(x) %in% c("a", "c", "t", "g"))]))
seq_selection_criteria <- (sapply(fastas, length)>=29000)&(ambiguous_prop<=.1)
output_fasta <- fastas[seq_selection_criteria]
output_fasta$outgroup <- outgroup_fasta
class(output_fasta) <- "DNAbin"

print(head(seq_metadata))

print(colnames(seq_metadata))

seq_metadata <- filter(seq_metadata, host=="Human")

# split sequences by location
output_fasta_provinces <- filter(seq_metadata, division_exposure %in% provinces_to_model) %>%
  split(., .$division_exposure) %>%
  lapply(function (x) output_fasta[names(output_fasta) %in% x$strain]) %>%
  `[`(., names(.) %in% provinces_to_model)

output_fasta_countries <- with(seq_metadata, coalesce(country_exposure, country)) %>%
  split(seq_metadata, .) %>%
  lapply(function (x) output_fasta[names(output_fasta) %in% x$strain]) %>%
  `[`(., names(.) %in% countries_to_model)

output_fasta_location <- c(output_fasta_provinces, output_fasta_countries)

# output files

time_str <- Sys.time() %>% as.character() %>% gsub("[^[:alnum:]]", "", .)

mclapply(names(output_fasta_location), function (loc_name) {
  output_fasta <- output_fasta_location[[loc_name]]
  loc_name %<>% gsub(" ", "", .) %>% tolower()
  msa_dir <- file.path("msa", loc_name)
  if (!dir.exists(msa_dir)) {
    dir.create(msa_dir, recursive=TRUE)
  }
  existing_msa <- list.files(msa_dir, "msa_mafft_", full.names=TRUE)
  if (length(existing_msa)>0) {
    # if there is already an existing MSA, only output sequences that have not already been processed
    processed_seq <- tail(existing_msa, 1) %>% read.dna("fasta") %>% rownames()
    output_fasta <- output_fasta[!(names(output_fasta) %in% processed_seq)]
  }
  write.dna(output_fasta, 
            file.path(msa_dir, paste0("input_mafft_", time_str, ".fasta")),
            format="fasta", nbcol=1, colw = 80)
}, mc.cores=detectCores())


save.image(paste0("02_filter_seq_", time_str, ".RData"))
