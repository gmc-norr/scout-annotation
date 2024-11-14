# This script annotates a vcf against a civic variant tsv file
# The matching is done by either rsID or protein change

from cyvcf2 import VCF, Writer
import pandas as pd
from snakemake import snakemake

def annotate_vcf(vcf_path, tsv_path, output_path):
    # Load the VCF file
    vcf = VCF(vcf_path)

    # Load Civic TSV file
    tsv = pd.read_csv(tsv_path, sep='\t')

    # Add new INFO fields to header
    vcf.add_info_to_header(
        {
            'ID': 'CIVIC_EVIDENCE',
            'Description': 'CIVIC evidence annotations',
            'Type': 'String',
            'Number': '.'
        }
    )

    writer = Writer(output_path, vcf)

    for variant in vcf:
        rs_id = variant.ID
        matches = tsv[tsv['rsID'] == rs_id]
        if not matches.empty:
            evidence = ','.join(matches['EVIDENCE'].astype(str))
            variant.INFO['CIVIC_EVIDENCE'] = evidence
        writer.write_record(variant)

    writer.close()
    vcf.close()

def main():
    vcf = snakemake.input.vcf
    civic_tsv = snakemake.input.civic_tsv
    annotate_vcf(
        snakemake.input.vcf, 
        snakemake.input.civic_variants_tsv, 
        snakemake.output.vcf
    )

if __name__ == "__main__":
    main()
