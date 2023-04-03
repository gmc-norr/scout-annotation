from pysam import VariantFile
import sys


def add_ad(variant):
    pass


def add_af(variant):
    pass


def main():
    input_vcf = snakemake.input.vcf
    output_vcf = snakemake.output.vcf
    log = snakemake.log

    vcf = VariantFile(input_vcf)

    header = vcf.header

    # Check for AF, DP, and AO in the header
    if not "AF" in vcf.header.formats:
        header.formats.add("AF", "A", "Float", "Alternative allele frequency")

    if not "AO" in vcf.header.formats:
        header.formats.add("AO", "A", "Integer", "Read depth for each alternative allele")

    ovcf = VariantFile(output_vcf, "w", header=header)

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
            v.samples[0]["AF"] = afs
        
        ovcf.write(v)

    ovcf.close()


if __name__ == "__main__":
    main()
