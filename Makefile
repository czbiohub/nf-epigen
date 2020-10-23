NEWICKS="timetrees/*.nwk"
AUSPICE_JSONS=jsons/
GISAID_METADATA=gisaid/metadata.tsv
GISAID_SEQUENCES= gisaid/sequences.tsv


run-nf-epigen-newicks:
	nextflow run main.nf \
	-profile docker \
	-resume \
	--newicks ${NEWICKS} \
	--gisaid_metadata ${GISAID_METADATA} \
	--gisaid_sequences ${GISAID_SEQUENCES}


