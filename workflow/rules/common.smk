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

def get_annotated_vcf(wildcards):
    snv_filter = get_sample_snv_filter_tag(wildcards)
    sample_panels = get_sample_panels(wildcards)

    if len(sample_panels) > 0:
        # Check the number of variants after panel filtering
        panel_filtering = checkpoints.panel_filtering.get(**wildcards)
        with panel_filtering.output["vcf"].open() as f:
            n_variants = sum(1 for line in f if not line.startswith("#"))
        if n_variants == 0:
            # No variants left for SNV filtering and ranking
            return f'{panel_filtering.output["vcf"]}.gz'

    if snv_filter is not None:
        # Check the number of variants after snv filtering
        vcf_filename = "annotation/{sample}/{sample}.annotated".format(**wildcards)
        if len(sample_panels) > 0:
            vcf_filename = "annotation/{sample}/{sample}.annotated.panel_filtered".format(**wildcards)
        vcf_filtering = checkpoints.vcf_filtering.get(
            file=vcf_filename,
            tag=get_sample_snv_filter_tag(wildcards),
        )
        with vcf_filtering.output["vcf"].open() as f:
            n_variants = sum(1 for line in f if not line.startswith("#"))
        if n_variants == 0:
            # No variants left for ranking
            return f'{vcf_filtering.output["vcf"]}.gz'

    # genmod-annotated VCF
    return "annotation/{sample}/{sample}.annotated.genmod.vcf.gz".format(**wildcards)

def get_annotated_vcf_index(wildcards):
    return f"{get_annotated_vcf(wildcards)}.tbi"

def get_filtered_vcf(wildcards):
    snv_filter = get_sample_snv_filter_tag(wildcards)
    sample_panels = get_sample_panels(wildcards)

    if snv_filter is None and len(sample_panels) == 0:
        # no filtering at all, return annotated VCF
        return "annotation/{sample}/{sample}.annotated.vcf"

    if snv_filter is None:
        # only panel filtering
        return "annotation/{sample}/{sample}.annotated.panel_filtered.vcf".format(**wildcards)

    # both panel and snv filtering
    return "annotation/{sample}/{sample}.annotated.panel_filtered.filter-{snv_filter}.vcf".format(**wildcards, snv_filter=snv_filter)


def get_bai_file(wildcards):
    bam = get_bam_file(wildcards)
    if not isinstance(bam, Path):
        return []
    return get_bam_file(wildcards).with_suffix(".bam.bai")

def get_bam_file(wildcards):
    sample_row = _get_sample_row(wildcards)
    if "bam" not in sample_row:
        return []
    bam = _get_sample_row(wildcards)["bam"].values
    assert len(bam) == 1
    if pd.isnull(bam[0]):
        return []
    return Path(bam[0]).resolve()

def get_panel_dict():
    panels = {}
    if not panel_path.exists():
        raise FileNotFoundError(f"directory not found: {panel_path}")
    for p in panel_path.glob("*.tsv"):
        panels[p.stem] = p
    return panels

def get_sample_panels(wildcards):
    panels = _get_sample_row(wildcards)["panels"].values
    if len(panels) == 0:
        return []
    assert len(panels) == 1
    if pd.isnull(panels[0]):
        return []
    return panels[0].split(",")

def get_sample_snv_filter_tag(wildcards):
    filtering = _get_sample_row(wildcards)["filtering"]
    assert len(filtering) == 1
    if pd.isnull(filtering.values[0]) or filtering.values[0] == "":
        return None
    return filtering.values[0]

def get_panel_files(wildcards):
    panel_dict = get_panel_dict()
    try:
        panel_files = [panel_dict[p] for p in get_sample_panels(wildcards) if p != ""]
    except KeyError as ie:
        raise KeyError(f"panel not found: {ie}")
    return panel_files

def get_vcf_file(wildcards):
    vcf_filename = _get_sample_row(wildcards)["vcf"].values
    if len(vcf_filename) == 0:
        raise RuntimeError("sample not found in sample file")
    if len(vcf_filename) != 1:
        raise RuntimeError("duplicate sample IDs found in sample file")
    return vcf_filename[0]

def get_reheadered_vcf_file(wildcards):
    vcf_sample_name = get_vcf_samples(wildcards)
    if vcf_sample_name != wildcards.sample:
        return rules.bcftools_reheader.output.vcf
    return get_vcf_file(wildcards)

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
        if isinstance(get_bam_file(snakemake.io.Namedlist(fromdict={"sample": s})), Path):
            outfiles.append(f"{outdir}/{s}/{s}.bam")
            outfiles.append(f"{outdir}/{s}/{s}.bam.bai")
        load_configs.append(f"{outdir}/{s}/{s}.load_config.yaml")
    return outfiles, load_configs

outfiles, load_configs = get_output_files()
