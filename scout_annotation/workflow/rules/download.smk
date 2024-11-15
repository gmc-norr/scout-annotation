include: "common.smk"

# Define Base URL for civic downloads
CIVIC_BASE_URL  = "https://civicdb.org/downloads/" 

# Define Base URL for Cancer Biomarkers Database https://www.cancergenomeinterpreter.org/data/biomarkers/cgi_biomarkers_latest.tsv
CGI_BASE_URL = "https://www.cancergenomeinterpreter.org/data/biomarkers/"

# Function to generate the full URL for downloading the CIVIC file
def get_civic_url():
    civic_file = get_resource_filename('civic')
    civic_version = config['annotation_references']['civic']['version']
    return f"{CIVIC_BASE_URL}{civic_version}/{civic_file}"

def get_cgi_url():
    cgi_file = get_resource_filename('cgi')
    cgi_version = config['annotation_references']['cgi']['version']
    return f"{CGI_BASE_URL}/{cgi_file}"

rule download_civic:
    output:
        civic_variants_tsv=f"{config['downloads_dir']}/{get_resource_filename('civic')}",
    run:
        civic_url_str = get_civic_url()
        shell(
            f"""
                wget -O {output.civic_variants_tsv} {civic_url_str}
            """
        )

rule download_cgi:
    output:
        cgi_variants_tsv=f"{config['downloads_dir']}/{get_resource_filename('cgi')}",
    run:
        cgi_url_str = get_cgi_url()
        shell(
            f"""
                wget -O {output.cgi_variants_tsv} {cgi_url_str}
            """
        )

