# This script annotates a vcf against a civic variant tsv file
# The matching is done by either rsID or protein change

from cyvcf2 import VCF, Writer
import pandas as pd
import re 
#from snakemake import snakemake

def annotate_vcf(vcf_path, tsv_path, output_path):
    # Load the VCF file
    vcf = VCF(vcf_path)

    # Load Civic TSV file
    civic_df = pd.read_csv(tsv_path, sep='\t', header=0)

    # Add new INFO fields to header
    vcf.add_info_to_header(
        {
            'ID': 'CIVIC_URL',
            'Description': 'URL for variant in CIVic',
            'Type': 'String',
            'Number': '.'
        }
    )

    writer = Writer(output_path, vcf)

    for variant in vcf:
        CSQ_fields = variant.INFO.get("CSQ").split(',')
        matching_civic = pd.DataFrame(columns=civic_df.columns)
        for CSQ_field in CSQ_fields:
            CSQ_dict = {key: value for key,value in zip(format_list, CSQ_field.split('|'))}
            hgvsc = CSQ_dict.get('HGVSc')
            matching = civic_df[
                civic_df['hgvs_descriptions'].notna() &  # Exclude NaN
                civic_df['hgvs_descriptions'].str.split(',').apply(lambda ids: hgvsc in ids if isinstance(ids, list) else False)
            ]

            if not matching.empty:
                matching_civic = pd.concat([matching_civic, matching],axis=0)
    
        if not matching_civic.empty:
            print(matching_civic['variant_civic_url'])
            variant.INFO['CIVIC_URL'] = matching_civic.iloc[0]['variant_civic_url']
        writer.write_record(variant)

    writer.close()
    vcf.close()

def main():
    annotate_vcf(
        snakemake.input.vcf, 
        snakemake.input.civic_variants_tsv, 
        snakemake.output.vcf
    ) 

if __name__ == "__main__":
    main()
