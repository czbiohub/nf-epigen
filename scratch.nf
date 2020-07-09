
process download_metadata {
    publishDir "${params.outdir}/gisaid", mode: 'copy'

    output:
    file "*metadata.tsv" into ch_metadata

    shell:
    """
    aws s3 cp s3://czb-covid-results/gisaid/metadata.tsv.gz - | gunzip -cq >metadata.tsv
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

    // output:
    // stdout "x" into ch_strings

    script:
    """
    01_curate_time_series.R \
        ${metadata} \
        ${deaths_us} \
        ${confirmed_us} \
        ${recovered_global} \
        ${deaths_global} \
        ${confirmed_global} \
        ${who_sitreps}
        ${wuhan_incidence}
    """
    // """
    // practice.R \
    //     ${metadata} \
    //     ${deaths_us} \
    //     ${confirmed_us} \
    //     ${recovered_global} \
    //     ${deaths_global} \
    //     ${confirmed_global} \
    //     ${who_sitreps} \
    //     ${wuhan_incidence}
    // """    
}


ch_cleaned_timeseries.view()
// ch_cleaned_timeseries.view()

