#!/usr/bin/env Rscript

library(dplyr)
library(EpiGenR)
library(magrittr)
library(readr)
library(phytools)
library(lubridate)

args <- commandArgs(TRUE)
timetree <- args[1] # ex: "~/code/nf-core-epigen/results/newicks/tree_california.nwk"
tip_dates_file <- args[2] # ex: "~/code/nf-core-epigen/results/gisaid/metadata.tsv"
cleaned_timeseries <- file.path(args[3]) # ex: "~/code/nf-core-epigen/results/timeseries_cleaned/summary_california_timeseries_new_cases.tsv"


# impute dates of infections
# -----------------------------------------------------
infer_dates_from_timeseries <- function (full_timeseries, lognormal_mean, lognormal_sd) {
  apply(full_timeseries, 1, function (x) {
    decimal_date(as.Date(x[1])) %>%
      `-`(rlnorm(x[2], lognormal_mean, lognormal_sd)) %>%
      date_decimal() %>%
      as.Date()
  }) %>%
    do.call(what=c)
}


timeseries <- cleaned_timeseries %>% 
  read_tsv %>%
  filter(new_cases>0)

incub_pars <- EpiGenR::get_lognormal_params(5.5/365, 2.1/365)
shape_param <- incub_pars[1]
scale_param <- incub_pars[2]
set.seed(2342342)
inferred_dates <- infer_dates_from_timeseries(timeseries, shape_param, scale_param)
output_timeseries <- get_data(epi=inferred_dates) %>% 
  setNames(c("date", "new_cases")) %>%
  mutate(date=as.Date(date_decimal(date)))

basename <- sub('\\..*$', '', basename(cleaned_timeseries))
write_csv(output_timeseries, paste0(basename, "_value_counts", ".txt"))
