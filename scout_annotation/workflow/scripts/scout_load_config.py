from pathlib import Path
import yaml


def get_track_name(track):
    match track:
        case "cancer":
            return "cancer"
        case "rare_disease" | "rare":
            return "rare"


def save_config(config, path):
    with open(path, "wt") as f:
        f.write(yaml.dump(config))


def phenotype(ped_path, sample):
    phenotypes = {
        -9: "missing",
        0: "missing",
        1: "unaffected",
        2: "affected",
    }

    with open(ped_path) as ped:
        for line in ped:
            _, s, _, _, _, p = line.strip().split()
            if s == sample:
                return phenotypes[int(p)]

    raise KeyError(f"no phenotype information found for sample: {sample}")


def generate_sample_config(sample):
    ped_file = snakemake.input.ped
    include_bam = snakemake.params.include_bam
    sex = snakemake.params.sex
    analysis_type = snakemake.params.analysis_type
    msi_score = snakemake.params.msi_score
    hrd_score = snakemake.params.hrd_score
    tmb_score = snakemake.params.tmb_score

    sample_config = {
        "sample_id": sample,
        "phenotype": phenotype(ped_file, sample),
        "sex": sex,
        "analysis_type": analysis_type,
    }

    if msi_score is not None:
        sample_config["msi"] = int(msi_score)
    if hrd_score is not None:
        sample_config["hrd"] = int(hrd_score)
    if tmb_score is not None:
        sample_config["tmb"] = int(tmb_score)

    if include_bam:
        sample_config["alignment_path"] = "{}.bam".format(snakemake.wildcards.sample)

    return sample_config


def generate_family_config(family, sample_config_files):
    track = get_track_name(snakemake.params.track)
    madeline2_svg = snakemake.input.madeline2_svg
    owner = snakemake.params.owner
    peddy_ped = snakemake.input.peddy_ped
    peddy_ped_check = snakemake.input.peddy_ped_check
    peddy_sex_check = snakemake.input.peddy_sex_check
    vcf = Path(snakemake.input.vcf)
    gene_panels = snakemake.params.panels
    rank_model_version = snakemake.params.rank_model_version
    rank_score_threshold = snakemake.params.rank_score_threshold

    sample_configs = []
    for fn in sample_config_files:
        with open(fn) as f:
            sample_configs.append(yaml.safe_load(f))

    family_config = {
        "family": family,
        "track": track,
        "owner": owner,
        "human_genome_build": "37",
        "gene_panels": gene_panels,
        "samples": sample_configs,
        "rank_model_version": rank_model_version,
        "rank_score_threshold": rank_score_threshold,
    }

    match track:
        case "rare":
            family_config["vcf_snv"] = vcf.name
        case "cancer":
            family_config["vcf_cancer"] = vcf.name

    if len(peddy_ped) > 0:
        family_config["peddy_ped"] = Path(peddy_ped).name
    if len(peddy_ped_check) > 0:
        family_config["peddy_check"] = Path(peddy_ped_check).name
    if len(peddy_sex_check) > 0:
        family_config["peddy_sex"] = Path(peddy_sex_check).name
    if len(madeline2_svg) > 0:
        family_config["madeline"] = Path(madeline2_svg).name

    return family_config


def main():
    config_type = snakemake.params.type

    output_filename = snakemake.output.yaml

    match config_type:
        case "sample":
            sample = snakemake.wildcards.sample
            config = generate_sample_config(sample)
        case "family":
            family = snakemake.wildcards.family
            sample_configs = snakemake.input.sample_configs
            config = generate_family_config(family, sample_configs)
        case _:
            raise Exception(f"invalid config type: {config_type}")

    save_config(config, output_filename)


if __name__ == "__main__":
    main()
