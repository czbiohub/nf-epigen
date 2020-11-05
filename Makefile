#NEWICKS="timetree_error/*.nwk"
NEWICKS="timetrees/*.nwk"
AUSPICE_JSONS=jsons/
GISAID_METADATA=gisaid/metadata.tsv
GISAID_SEQUENCES=gisaid/sequences.tsv
POPULATION_CSV=model_params/population.csv
NUM_SAMPLES_CSV=model_params/num_samples.csv

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
	--mortality_rate 2
	--incubation_time 5.5 \
	--recovery_time 14 \
	--haar_full_mass 7
