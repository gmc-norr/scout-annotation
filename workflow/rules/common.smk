import cyvcf2
import numpy as np
import pandas as pd
from pathlib import Path
from snakemake.utils import validate
import sys
import yaml

ruleorder: copy_results > tabix

if not workflow.overwrite_configfiles:
    print("error: ", file=sys.stderr)

validate(config, "../schema/config.schema.yaml")

with open(config["resources"]) as f:
    resources = yaml.load(f, Loader=yaml.FullLoader)
validate(resources, "../schema/resources.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", comment="#")
validate(samples, "../schema/samples.schema.yaml")

panel_path = Path(config["panel_filtering"]["panel_directory"])
if not panel_path.exists() and not panel_path.is_absolute():
    panel_path = Path(snakemake.workflow.srcdir("../.."), panel_path).resolve()

wildcard_constraints:
    ext=r"vcf(\.gz)?$",
    track=r"(rare_disease|cancer)"

def _get_sample_row(wildcards):
    return samples[samples["sample"] == wildcards.sample]

def get_bai_file(wildcards):
    bam = get_bam_file(wildcards)
    if len(bam) == 0:
        return []
    return f"{get_bam_file(wildcards)}.bai"

def get_bam_file(wildcards):
    sample_row = _get_sample_row(wildcards)
    if "bam" not in sample_row:
        return []
    bam = _get_sample_row(wildcards)["bam"].values
    assert len(bam) == 1
    if pd.isnull(bam[0]):
        return []
    return bam[0]

def get_filter_tag(wildcards):
    track = _get_sample_row(wildcards)["track"].values
    assert len(track) == 1
    if track == "rare_disease":
        return "rare_disease"
    elif track == "cancer":
        return "somatic"
    else:
        raise ValueError(f"no filter for track {track}")

def get_panel_dict():
    panels = {}
    if not panel_path.exists():
        raise FileNotFoundError(f"directory not found: {panel_path}")
    for p in panel_path.glob("*.tsv"):
        panels[p.stem] = p
    return panels

def get_sample_panels(wildcards):
    panels = _get_sample_row(wildcards)["panels"].values
    assert len(panels) == 1
    if pd.isnull(panels[0]):
        return []
    return panels[0].split(",")

def get_panel_files(wildcards):
    panel_dict = get_panel_dict()
    try:
        panel_files = [panel_dict[p] for p in get_sample_panels(wildcards)]
    except KeyError as ie:
        raise KeyError(f"panel not found: {ie}")
    return panel_files

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

def get_vcf_filter(wildcards):
    filter_definition = config.get("vcf_filter", {}).get(wildcards.tag)
    if not filter_definition:
        raise KeyError(f"no such vcf filter definition: {wildcards.tag}")
    return filter_definition

def get_vcfanno_config(wildcards):
    sample_track = get_track(wildcards)
    genome_build = config["genome_build"]
    # TODO: set this up for structural variants
    version = config["vcfanno"]["config_version"][sample_track]
    return f"rank_model/grch{genome_build}_{sample_track}_vcfanno_config_{version}.toml"

def get_vembrane_expression(wildcards):
    filter_file = get_vcf_filter(wildcards)
    with open(filter_file) as f:
        definition = yaml.load(f, Loader=yaml.FullLoader)
        expressions = []
        for filter in definition["filters"]:
            expressions.append(" ".join([x.strip() for x in filter["expression"].splitlines()]))
    return "(" + ") and (".join(expressions) + ")"

def get_output_files():
    outfiles = []
    load_configs = []
    outdir = config.get("output_directory", "results")
    for s, p in zip(samples["sample"], samples["ped"]):
        outfiles.append(f"{outdir}/{s}/{s}.scout-annotated.vcf.gz")
        outfiles.append(f"{outdir}/{s}/{s}.scout-annotated.vcf.gz.tbi")
        outfiles.append(f"{outdir}/{s}/{s}.ped")
        if len(get_bam_file(snakemake.io.Namedlist(fromdict={"sample": s}))) != 0:
            outfiles.append(f"{outdir}/{s}/{s}.bam")
            outfiles.append(f"{outdir}/{s}/{s}.bam.bai")
        load_configs.append(f"{outdir}/{s}/{s}.load_config.yaml")
    return outfiles, load_configs

outfiles, load_configs = get_output_files()
