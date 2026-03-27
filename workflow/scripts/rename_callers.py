import pysam
import logging

logger = logging.getLogger(__name__)

CALLERS_MAP = {"gatk_mutect2": "mutect"}

CALLERS_FIELD = "FOUND_IN"


def _rename_callers_in_rec(
    rec: pysam.VariantRecord, callers_map: dict, callers_field: str
) -> None:
    try:
        callers = list(rec.info[callers_field])
    except KeyError:
        return
    new_callers = list()
    for caller in callers:
        if caller in callers_map:
            new_callers.append(callers_map[caller])
        else:
            new_callers.append(caller)
    rec.info[callers_field] = new_callers


def rename_callers(
    vcf_in: str,
    vcf_out: str,
    callers_map: dict = CALLERS_MAP,
    callers_field: str = CALLERS_FIELD,
) -> None:
    logger.info(f"renaming {callers_field} field in INFO of {vcf_in}")
    for key, value in callers_map.items():
        logger.info(f"  {key} -> {value}")
    logger.info(f"and writing to {vcf_out}")

    with pysam.VariantFile(vcf_in, "r") as in_vcf:
        with pysam.VariantFile(vcf_out, "w", header=in_vcf.header) as out_vcf:
            for record in in_vcf:
                _rename_callers_in_rec(record, callers_map, callers_field)
                out_vcf.write(record)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    try:
        log = snakemake.log_fmt_shell(stdout=False, stderr=True)
        rename_callers(
            snakemake.input.vcf,
            snakemake.output.vcf,
            snakemake.params.callers_map,
            snakemake.params.callers_field,
        )
    except NameError:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("-i", "--input-vcf", required=True)
        p.add_argument("-o", "--output-vcf", required=True)
        a = p.parse_args()

        rename_callers(a.input_vcf, a.output_vcf)
