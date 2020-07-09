
process download_gisaid_metadata_and_sequences {
    publishDir "${params.outdir}/gisaid", mode: 'copy'

    output:
    file "*metadata.tsv" into ch_metadata
    file "*sequences.fasta" into ch_sequences

    shell:
    """
    aws s3 cp s3://czb-covid-results/gisaid/metadata.tsv.gz - | gunzip -cq >metadata.tsv
    aws s3 cp s3://czb-covid-results/gisaid/sequences.fasta.gz - | gunzip -cq >sequences.fasta
    """
}

process download_timeseries {
    publishDir "${params.outdir}/timeseries", mode: 'copy'

    output:
    file "*time_series_covid19_deaths_US.csv" into ch_timeseries_deaths_us
    file "*time_series_covid19_confirmed_US.csv" into ch_timeseries_confirmed_us
    file "*time_series_covid19_recovered_global.csv" into ch_timeseries_recovered_global
    file "*time_series_covid19_deaths_global.csv" into ch_timeseries_deaths_global
    file "*time_series_covid19_confirmed_global.csv" into ch_timeseries_confirmed_global

    shell:
    """
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv
    """
}

process clean_and_transform_timeseries {
    publishDir "${params.outdir}/timeseries_cleaned", mode: 'copy'

    input:
    file deaths_us from ch_timeseries_deaths_us
    file confirmed_us from ch_timeseries_confirmed_us
    file recovered_global from ch_timeseries_recovered_global
    file deaths_global from ch_timeseries_deaths_global
    file confirmed_global from ch_timeseries_confirmed_global
    file metadata from ch_metadata
    file who_sitreps from file("$baseDir/datasets/WHO_sitreps_20200121-20200122.tsv")
    file wuhan_incidence from file("$baseDir/datasets/li2020nejm_wuhan_incidence.tsv")

    output:
    file "*.tsv" into ch_cleaned_timeseries
    file "*.RData" into ch_rdata_cleaned_timeseries

    script:
    """
    01_curate_time_series.R \
        ${metadata} \
        ${deaths_us} \
        ${confirmed_us} \
        ${recovered_global} \
        ${deaths_global} \
        ${confirmed_global} \
        ${who_sitreps} \
        ${wuhan_incidence}
    """
}


process rename_sequences_to_include_collection_dates {
    publishDir "${params.outdir}/renamed_sequences", mode: 'copy'

    input:
    file metadata from ch_metadata
    file sequences from ch_sequences
    file rdata from ch_rdata_cleaned_timeseries
    file outgroup_fasta from file("$baseDir/datasets/MG772933.1.fasta")

    output:
    file "*.fasta" into ch_renamed_fastas
    file "*.RData" into ch_rdata_renamed_sequences

    script:
    """
    02_filter_seq.R \
        ${rdata} \
        ${metadata} \
        ${sequences} \
        ${outgroup_fasta}
    """
}

ch_renamed_fastas.view()


