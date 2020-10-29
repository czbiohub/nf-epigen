#NEWICKS="timetree_error/*.nwk"
NEWICKS="timetrees/*.nwk"
AUSPICE_JSONS=jsons/
GISAID_METADATA=gisaid/metadata.tsv
GISAID_SEQUENCES=gisaid/sequences.tsv
POPULATION_CSV=model_params/population.csv
NUM_SAMPLES_CSV=model_params/num_samples.csv

run-nf-epigen-newicks:
	nextflow run main.nf \
	-profile docker \
	-resume \
	--newicks ${NEWICKS} \
	--gisaid_metadata ${GISAID_METADATA} \
	--gisaid_sequences ${GISAID_SEQUENCES} \
	--population_csv ${POPULATION_CSV} \
	--num_samples_csv ${NUM_SAMPLES_CSV}

