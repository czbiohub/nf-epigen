#!/usr/bin/env python3

import json
from datetime import datetime as dt
import logging
import time
import argparse

import torch
from Bio import Phylo
import pandas as pd
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


def prune_tree(phylogeny, metadata, region):

    # get strain names from metadata
    strains = metadata.query(
        f"country == '{region.title()}' | \
        division == '{region.title()}'"
    )["strain"].values

    logging.info(f"beginning terminal nodes {phylogeny.count_terminals()}")

    # collapse internal nodes that have data not from region
    internals = phylogeny.get_nonterminals()
    for node in internals:
        if (not(node.name.startswith("NODE_"))&(node.name not in strains)):
            phylogeny.collapse(node)

    # prune terminal nodes
    terminals = phylogeny.get_terminals()
    for node in terminals:
        if (node.name not in strains):
            phylogeny.prune(node)

    logging.info(f"ending terminal nodes {phylogeny.count_terminals()}")

    for node in phylogeny.find_clades():
        if node.is_terminal():
            logging.info("node: {node.name} retained in tree")


def get_new_cases(infection_dates):
    new_cases = list(infection_dates["new_cases"])
    if new_cases[-1] == 0:
        new_cases.pop(-1)
    return torch.tensor(new_cases, dtype=torch.double)


def get_leaf_time_coal_times(phylogeny, last_tip_date, infection_dates):
    leaf_times, coal_times = dist.coalescent.bio_phylo_to_times(
        phylogeny
    )
    shift = last_tip_date - max(leaf_times)
    first_timeseries_date = pd.to_datetime(
        infection_dates["date"]
    ).apply(toYearFraction).min()

    leaf_times = (leaf_times + shift - first_timeseries_date)*365.25
    coal_times = (coal_times + shift - first_timeseries_date)*365.25

    return leaf_times, coal_times

def main():

    parser = argparse.ArgumentParser(description="SEIR model for specific region")
    parser.add_argument("--region", help="region of interest (division or country)")
    parser.add_argument("--tree", help="tree newick file")
    parser.add_argument(
        "--infection_dates", help="timeseries of new cases value counts"
    )
    parser.add_argument("--metadata", help="metadata.tsv file from gisaid", required=True)

    parser.add_argument(
        "--model_type",
        help="stochastic discrete-time discrete-state model",
        choices=[
            "SimpleSIRModel",
            "SimpleSEIRModel",
            "SimpleSEIRDModel",
            "OverdispersedSIRModel",
            "OverdispersedSEIRModel",
            "SuperspreadingSIRModel",
            "SuperspreadingSEIRModel",
            "HeterogeneousSIRModel"
        ],
        required=True,
        type=str
    )
    parser.add_argument(
        "--population",
        help="the total population of a single-region (S + I + R)",
        required=True,
        type=float
    )
    parser.add_argument(
        "--num_samples",
        help="Number of posterior samples to draw via mcmc",
        required=True,
        type=int
    )
    parser.add_argument(
        "--haar_full_mass",
        help="Number of low frequency Haar components to include in the full mass matrix.",
        required=True,
        type=int
    )
    parser.add_argument(
        "--incubation_time",
        help="Mean incubation time (duration in state E)",
        required=False,
        type=float
    )
    parser.add_argument(
        "--recovery_time",
        help="Mean recovery time (duration in state I)",
        required=True,
        type=float
    )
    parser.add_argument(
        "--mortality_rate",
        help="Mean mortality rate",
        type=float
    )

    args = parser.parse_args()

    phylogeny = get_phylogeny(args.tree, tree_type="newick")
    infection_dates = pd.read_csv(args.infection_dates, sep=",")
    full_metadata = pd.read_csv(args.metadata, sep="\t")

    fix_parsing_annotations(phylogeny)
    prune_tree(phylogeny, full_metadata, args.region)

    new_cases = get_new_cases(infection_dates)

    metadata = full_metadata[
        full_metadata["strain"].isin([x.name for x in phylogeny.get_terminals()])
    ]
    metadata["date"] = pd.to_datetime(metadata["date"])
    metadata["decimal_date"] = metadata["date"].apply(toYearFraction)

    compartment_model = getattr(pyro.contrib.epidemiology.models, args.model_type)
    compartment_model_params = {
        "population": int(args.population),
        "recovery_time": args.recovery_time,
        "data": new_cases,
    }

    # add optional kwargs
    if args.incubation_time:
        compartment_model_params["incubation_time"] = args.incubation_time

    if args.mortality_rate:
        compartment_model_params["mortality_rate"] = args.mortality_rate

    if args.model_type == "SuperspreadingSEIRModel":
        last_tip_date = metadata["decimal_date"].max()
        leaf_times, coal_times = get_leaf_time_coal_times(
            phylogeny, last_tip_date, infection_dates
        )
        compartment_model_params["leaf_times"] = leaf_times
        compartment_model_params["coal_times"] = coal_times

    model = compartment_model(**compartment_model_params)

    mcmc_fit_params = {
        "num_samples": args.num_samples,
        "haar_full_mass": args.haar_full_mass
    }
    mcmc = model.fit_mcmc(**mcmc_fit_params)

    mcmc.summary()


if __name__ == "__main__":
    main()
