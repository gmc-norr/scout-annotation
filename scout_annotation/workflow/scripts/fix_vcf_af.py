import pathlib
from typing import List, Union
from pysam import VariantFile, VariantRecord, VariantHeader


def fix_vcf_af(vcf: VariantFile) -> List[VariantRecord]:
    header = vcf.header

    af_ao_number = "A"

    # Set the AF/AO headers if they don't exist, and set the
    # Number parameter accordingly.
    if not "AF" in vcf.header.formats and "AO" not in vcf.header.formats:
        header.formats.add(
            "AF", af_ao_number, "Float", "Allele frequency of each alternative allele"
        )
        header.formats.add(
            "AO", af_ao_number, "Integer", "Read depth for each alternative allele"
        )

    if "AF" in vcf.header.formats and "AO" in vcf.header.formats:
        af_ao_number = str(header.formats["AF"].number)
    elif "AF" in vcf.header.formats:
        af_ao_number = str(header.formats["AF"].number)
        header.formats.add(
            "AO", af_ao_number, "Integer", "Alternative allele observations"
        )
    elif "AO" in vcf.header.formats:
        af_ao_number = str(header.formats["AO"].number)
        header.formats.add(
            "AF", af_ao_number, "Float", "Allele frequency of the alternative allele"
        )

    modified_variants = []

    for v in vcf:
        n_alternate = len(v.alts) if v.alts is not None else None

        if "AO" not in v.format:
            # Cannot tell with the information at hand, set to missing
            if n_alternate is None and af_ao_number == "1":
                v.samples[0]["AO"] = None
            elif n_alternate is None:
                v.samples[0]["AO"] = tuple()
            else:
                v.samples[0]["AO"] = tuple(None for _ in range(n_alternate))

        if "DP" not in v.format:
            # Cannot tell with the information at hand, set to missing
            v.samples[0]["DP"] = None

        if (
            "AF" not in v.format
            and v.samples[0]["DP"] is not None
            and v.samples[0]["AO"] is not None
        ):
            # Calculate the alternative allele frequency as the number of
            # reads representing the allele divided by the total read depth
            afs = []

            if af_ao_number == "A":
                for ao in v.samples[0]["AO"]:
                    if ao is None:
                        afs.append(None)
                    else:
                        afs.append(ao / v.samples[0]["DP"])
            elif af_ao_number == "1":
                if isinstance(v.samples[0]["AO"], tuple):
                    # Weird cases where the AO number is supposed to be 1, but
                    # multiple values still occur.
                    v.samples[0]["AO"] = None
                    afs.append(None)
                else:
                    afs.append(v.samples[0]["AO"] / v.samples[0]["DP"])
            else:
                raise TypeError(f"invalid Number for AF/AO: {af_ao_number}")

            if n_alternate == 0:
                assert len(afs) == 0
                v.samples[0]["AF"] = tuple()
            else:
                v.samples[0]["AF"] = afs
        elif "AF" not in v.format or v.samples[0]["AF"] is None:
            v.samples[0]["AF"] = None
        else:
            # AF is already defined, add variant as-is
            pass

        modified_variants.append(v)

    return modified_variants


def write_variants(
    filename: Union[str, pathlib.Path],
    variants: List[VariantRecord],
    header: VariantHeader,
):
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
