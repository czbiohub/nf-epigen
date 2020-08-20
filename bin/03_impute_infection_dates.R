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


# reformat labels on trees to include date of sampling
# -----------------------------------------------------
trs <- mapply(function (x, y) {
  tr <- read.newick(x)
  node_dates <- read_tsv(y)
  tr <- keep.tip(tr, which(tr$tip.label %in% node_dates$gisaid_epi_isl))
  tr$tip.label %<>% paste(node_dates$date[match(., node_dates$gisaid_epi_isl)], sep="_")
  tr
}, timetree, tip_dates_file, SIMPLIFY=FALSE) %>%
  setNames(., strsplit(names(.), "/") %>% sapply(`[`, 2))

# Exclude because there is no dominant lineage and the TMRCA is too far back: Ontario
# Exclude because of insufficient time-series data: Austria, Belgium
# Check TMRCA (should be <0.5 years)
# lapply(skylines, `[[`, "time") %>% sapply(max)
# Select dominant lineage: California, France, Guangdong, Hong Kong, Hubei, Italy, 
# Japan, Shanghai

select_nodes <- c("california"="NODE_0000002",
                  "france"="NODE_0000017",
                  "guangdong"="NODE_0000012",
                  "hongkong"="NODE_0000001",
                  "hubei"="NODE_0000005",
                  "italy"="NODE_0000008",
                  "japan"="NODE_0000005",
                  "shanghai"="NODE_0000001"
)

selected_trs <- trs[!(names(trs) %in% c("austria", "belgium", "ontario"))] 


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
  lapply(read_tsv) %>%
  lapply(filter, new_cases>0) %>%
  setNames(., names(selected_trs))

incub_pars <- EpiGenR::get_lognormal_params(5.5/365, 2.1/365)
shape_param <- incub_pars[1]
scale_param <- incub_pars[2]
set.seed(2342342)
inferred_dates <- table(lapply(timeseries, infer_dates_from_timeseries, shape_param, scale_param))
output_timeseries <- get_data(epi=inferred_dates) %>% 
  setNames(c("date", "new_cases")) %>%
  mutate(date=as.Date(date_decimal(date)))

basename <- sub('\\..*$', '', basename(cleaned_timeseries))
write_csv(output_timeseries, paste0(basename, "_value_counts", ".txt"))
