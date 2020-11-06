#NEWICKS="timetree_error/*.nwk"
NEWICKS="timetrees/*.nwk"
AUSPICE_JSONS=jsons/
GISAID_METADATA=gisaid/metadata.tsv
GISAID_SEQUENCES=gisaid/sequences.tsv
POPULATION_CSV=model_params/population.csv
NUM_SAMPLES_CSV=model_params/num_samples.csv

test-python:
	python bin/seir_model.py \
	--metadata ${GISAID_METADATA} \
	--infection_dates results/imputed_infection_dates/summary_california_timeseries_new_cases_value_counts.txt \
	--tree results/imputed_infection_dates/tree_california.nwk \
	--model_type SuperspreadingSEIRModel \
	--population 39e7 \
	--num_samples 100 \
	--incubation_time 5.5 \
	--recovery_time 14 \
	--haar_full_mass 7

run-nf-epigen-superspreadingseirmodel:
	nextflow run main.nf \
	-profile docker \
	-resume \
	--newicks ${NEWICKS} \
	--gisaid_metadata ${GISAID_METADATA} \
	--gisaid_sequences ${GISAID_SEQUENCES} \
	--model_type SuperspreadingSEIRModel \
	--population_csv ${POPULATION_CSV} \
	--num_samples_csv ${NUM_SAMPLES_CSV} \
	--incubation_time 5.5 \
	--recovery_time 14 \
	--haar_full_mass 7


run-nf-epigen-simpleseirdmodel:
	nextflow run main.nf \
	-profile docker \
	-resume \
	--newicks ${NEWICKS} \
	--gisaid_metadata ${GISAID_METADATA} \
	--gisaid_sequences ${GISAID_SEQUENCES} \
	--model_type SimpleSEIRDModel \
	--population_csv ${POPULATION_CSV} \
	--num_samples_csv ${NUM_SAMPLES_CSV} \
	--mortality_rate 0.01
	--incubation_time 5.5 \
	--recovery_time 14 \
	--haar_full_mass 7
