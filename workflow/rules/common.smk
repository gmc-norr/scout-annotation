import cyvcf2
import numpy as np
import pandas as pd
from pathlib import Path
from snakemake.utils import validate
import yaml

configfile: "config/config.yaml"
validate(config, "../schema/config.schema.yaml")

with open(config["resources"]) as f:
    resources = yaml.load(f, Loader=yaml.FullLoader)
validate(resources, "../schema/resources.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", comment="#")
validate(samples, "../schema/samples.schema.yaml")

wildcard_constraints:
    ext=r"vcf(\.gz)?$",
    track=r"(rare_disease|cancer)"

def _get_sample_row(wildcards):
    return samples[samples["sample"] == wildcards.sample]

def get_vcf_file(wildcards):
    vcf_filename = _get_sample_row(wildcards)["vcf"].values
    assert len(vcf_filename) == 1, "duplicate sample IDs found"
    return vcf_filename[0]

def get_ped(wildcards):
    ped = _get_sample_row(wildcards)["ped"].values
    assert len(ped) == 1
    if not ped[0]:
        return ped
    return "mock_ped/{sample}.ped".format(**wildcards)

def get_sex(wildcards):
    sex = _get_sample_row(wildcards)["sex"]
    assert len(sex) == 1
    return sex.values[0]

def get_analysis_type(wildcards):
    analysis_type = _get_sample_row(wildcards)["type"]
    assert len(analysis_type) == 1
    return analysis_type.values[0]

def get_track(wildcards):
    track = _get_sample_row(wildcards)["track"]
    assert len(track) == 1
    return track.values[0]

def get_vcf_samples(wildcards):
    vcf_filename = _get_sample_row(wildcards)["vcf"]
    assert len(vcf_filename) == 1
    vcf = cyvcf2.VCF(vcf_filename.values[0])
    samples = vcf.samples
    assert len(samples) == 1
    return samples[0]

def get_rank_model_version(wildcards):
    sample_track = get_track(wildcards)
    # TODO: set this up for structural variants
    return config["genmod"]["rank_model_version"][sample_track]

def get_rank_model(wildcards):
    sample_track = get_track(wildcards)
    # TODO: set this up for structural variants
    version = get_rank_model_version(wildcards)
    return f"rank_model/{sample_track}_rank_model_{version}.ini"

def get_vcfanno_config(wildcards):
    sample_track = get_track(wildcards)
    genome_build = config["genome_build"]
    # TODO: set this up for structural variants
    version = config["vcfanno"]["config_version"][sample_track]
    return f"rank_model/grch{genome_build}_{sample_track}_vcfanno_config_{version}.toml"

def get_result_files():
    infiles = []
    outfiles = []
    load_configs = []
    outdir = config.get("output_directory", "results")
    for s, p in zip(samples["sample"], samples["ped"]):
        infiles.append(f"annotation/{s}/{s}.decomposed.vep.most_severe_csq.vcfanno.genmod.vcf.gz")
        outfiles.append(f"{outdir}/{s}/{s}.scout-annotated.vcf.gz")

        infiles.append(f"annotation/{s}/{s}.decomposed.vep.most_severe_csq.vcfanno.genmod.vcf.gz.tbi")
        outfiles.append(f"{outdir}/{s}/{s}.scout-annotated.vcf.gz.tbi")

        if isinstance(p, str) and len(p) > 0:
            infiles.append(p)
        else:
            infiles.append(f"mock_ped/{s}.ped")
        outfiles.append(f"{outdir}/{s}/{s}.ped")

        load_configs.append(f"{outdir}/{s}/{s}.load_config.yaml")

    return infiles, outfiles, load_configs

infiles, outfiles, load_configs = get_result_files()
