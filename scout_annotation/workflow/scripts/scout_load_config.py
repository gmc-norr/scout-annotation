import cyvcf2
import pathlib
import yaml


def get_track_name(track):
    match track:
        case "cancer":
            return "cancer"
        case "rare_disease" | "rare":
            return "rare"


def generate_load_config(vcf, ped):
    genome_build = snakemake.config["genome_build"]
    rank_model_version = snakemake.params["rank_model_version"]
    owner = snakemake.params["owner"]

    sample_ids = cyvcf2.VCF(vcf).samples
    assert len(sample_ids) == 1, "multiple samples found in VCF"

    load_config = dict(
        family=snakemake.params["sample_name"],
        genome_build=genome_build,
        rank_model_version=rank_model_version,
        rank_score_threshold=-1000,
        owner=owner,
        track=get_track_name(snakemake.params["track"]),
        gene_panels=snakemake.params["panels"],
        samples=[
            dict(
                sample_id=sample_ids[0],
                analysis_type=snakemake.params["analysis_type"],
                phenotype=snakemake.params["phenotype"],
                sex=snakemake.params["sex"],
            ),
        ],
    )

    if snakemake.params["track"] == "cancer":
        load_config["vcf_cancer"] = vcf.name
    else:
        load_config["vcf_snv"] = vcf.name

    if snakemake.params["include_bam"]:
        bam = f"{snakemake.params['sample_name']}.bam"
        load_config["samples"][0]["alignment_path"] = bam

    return load_config


def save_config(config, path):
    with open(path, "wt") as f:
        f.write(yaml.dump(config))


def main():
    vcf_filename = pathlib.Path(snakemake.input["vcf"])
    ped_filename = pathlib.Path(snakemake.input["ped"])
    yaml_filename = pathlib.Path(snakemake.output["yaml"])

    load_config = generate_load_config(vcf_filename, ped_filename)

    save_config(load_config, yaml_filename)


if __name__ == "__main__":
    main()