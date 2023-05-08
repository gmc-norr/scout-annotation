import pathlib
from typing import List, Union
from pysam import VariantFile, VariantRecord, VariantHeader
import sys

def fix_vcf_af(vcf: VariantFile) -> List[VariantRecord]:
    header = vcf.header

    # Check for AF, DP, and AO in the header
    if not "AF" in vcf.header.formats:
        header.formats.add("AF", "A", "Float", "Alternative allele frequency")

    if not "AO" in vcf.header.formats:
        header.formats.add("AO", "A", "Integer", "Read depth for each alternative allele")

    modified_variants = []

    for v in vcf:
        if "AO" not in v.format:
            # Cannot tell with the information at hand, set to missing
            v.samples[0]["AO"] = None
        if "DP" not in v.format:
            # Cannot tell with the information at hand, set to missing
            v.samples[0]["DP"] = None
        if "AF" not in v.format:
            # Calculate the alternative allele frequency as the number of
            # reads representing the allele divided by the total read depth
            afs = []
            if v.samples[0]["AO"][0] is None:
                afs.append(None)
            else:
                for ao in v.samples[0]["AO"]:
                    afs.append(ao / v.samples[0]["DP"])
            if len(afs) > 1:
                v.samples[0]["AF"] = None
            else:
                v.samples[0]["AF"] = afs
        
        modified_variants.append(v)

    return modified_variants


def write_variants(filename: Union[str, pathlib.Path], variants: List[VariantRecord], header: VariantHeader):
    out_vcf = VariantFile(filename, "w", header=header)
    for v in variants:
        out_vcf.write(v)
    out_vcf.close()


def main():
    input_vcf = snakemake.input.vcf
    output_vcf = snakemake.output.vcf
    log = snakemake.log

    vcf = VariantFile(input_vcf)

    variants = fix_vcf_af(vcf)
    write_variants(output_vcf, variants, vcf.header)

if __name__ == "__main__":
    main()
