#!/usr/bin/bash

mkdir newicks
mkdir gisaid

# download newicks
aws s3 cp s3://czb-covid-results/gisaid/results/ newicks/ --recursive --exclude '*' --include '*.nwk'

# download gisaid sequences and metadata
aws s3 cp s3://czb-covid-results/gisaid/metadata.tsv.gz - | gunzip -cq >gisaid/metadata.tsv
aws s3 cp s3://czb-covid-results/gisaid/sequences.fasta.gz - | gunzip -cq >gisaid/sequences.fasta
