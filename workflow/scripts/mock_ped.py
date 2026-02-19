import gzip
import logging
import pathlib


def get_sample(vcf_filename):
    openfun = open
    if vcf_filename.suffix == ".gz":
        openfun = gzip.open

    with openfun(vcf_filename, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.strip().split()[9:]
                assert len(samples) == 1, "only single-sample vcf files supported"
                return samples[0]

    raise ValueError(f"no header line found: {vcf_filename}")


def generate_ped_line(family, sample):
    ped_line = [
        family,
        sample,
        0,  # paternal id
        0,  # maternal id
        0,  # sex
        2,  # affected
    ]

    return "\t".join(map(str, ped_line)) + "\n"


def write_ped(line, ped_filename):
    with open(ped_filename, "w") as f:
        f.write(line)


def main():
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])

    vcf_filename = pathlib.Path(snakemake.input.vcf)
    ped_filename = pathlib.Path(snakemake.output.ped)

    logging.info(f"Creating mock ped for {vcf_filename}")

    family_id = snakemake.wildcards.sample
    logging.info(f"Using family id '{family_id}'")

    sample = get_sample(vcf_filename)
    logging.info(f"sample: '{sample}'")

    ped_line = generate_ped_line(family_id, sample)
    write_ped(ped_line, ped_filename)

    logging.info(f"ped written to {ped_filename}")


if __name__ == "__main__":
    main()
