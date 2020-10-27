#!/usr/bin/env python3

import re
import datetime
import urllib.request
import io
import math
from datetime import date
import json
from datetime import datetime as dt
import time
import argparse

import torch
from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
from augur.utils import json_to_tree

import pyro
import pyro.distributions as dist
from pyro.contrib.epidemiology import CompartmentalModel
from pyro.contrib.epidemiology.models import SuperspreadingSEIRModel
from pyro.contrib.epidemiology.distributions import binomial_dist, infection_dist

pyro.enable_validation(True)
torch.set_default_dtype(torch.double)
torch.set_printoptions(precision=2)


def get_phylogeny(tree_path, tree_type="newick"):
    if tree_type == "newick":
        with open(tree_path) as f:
            for phylogeny in Phylo.parse(f, "newick"):
                break
            return phylogeny

    elif tree_type == "auspice_json":
        with open(tree_path) as f:
            tree = json.load(f)

        phylogeny = json_to_tree(tree)
        return phylogeny

    else:
        raise ValueError(f"current tree_type {tree_type} not supported")


def fix_parsing_annotations(phylogeny):
    # Fix a parsing error for whereby internal nodes interpret .name as .confidence
    for clade in phylogeny.find_clades():
        if clade.confidence:
            clade.name = clade.confidence
            clade.confidence = None


def toYearFraction(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year + 1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def main():

    parser = argparse.ArgumentParser(description=f"SEIR model for specific region")
    parser.add_argument("--tree", help="tree newick file or auspice json")
    parser.add_argument(
        "--infection_dates", help="timeseries of new cases value counts"
    )
    parser.add_argument("--metadata", help="metadata.tsv file from gisaid")

    args = parser.parse_args()

    phylogeny_newick = get_phylogeny(args.tree, tree_type="newick")
    infection_dates = pd.read_csv(args.infection_dates, sep=",")
    full_metadata = pd.read_csv(args.metadata, sep="\t")

    fix_parsing_annotations(phylogeny_newick)
    new_cases = list(infection_dates["new_cases"])
    if new_cases[-1] == 0:
        new_cases.pop(-1)
    new_cases = torch.tensor(new_cases, dtype=torch.double)

    metadata = full_metadata[
        full_metadata["strain"].isin([x.name for x in phylogeny_newick.get_terminals()])
    ]
    metadata["date"] = pd.to_datetime(metadata["date"])
    metadata["decimal_date"] = metadata["date"].apply(toYearFraction)

    last_tip_date = metadata["decimal_date"].max()
    leaf_times, coal_times = dist.coalescent.bio_phylo_to_times(phylogeny_newick)
    shift = last_tip_date - max(leaf_times)
    first_timeseries_date = pd.to_datetime(infection_dates["date"]).apply(toYearFraction).min()
    leaf_times = (leaf_times + shift - first_timeseries_date)*365.25
    coal_times = (coal_times + shift - first_timeseries_date)*365.25

    model = SuperspreadingSEIRModel(
        population=int(4e7),
        incubation_time=5.5,
        recovery_time=14.0,
        data=new_cases,
        leaf_times=leaf_times,
        coal_times=coal_times,
    )

    mcmc = model.fit_mcmc(num_samples=200, haar_full_mass=7)

    mcmc.summary()


if __name__ == "__main__":
    main()
