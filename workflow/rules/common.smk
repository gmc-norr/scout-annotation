import cyvcf2
import numpy as np
import pandas as pd
from pathlib import Path
import yaml

configfile: "config/config.yaml"
with open(config["resources"]) as f:
    resources = yaml.load(f, Loader=yaml.FullLoader)

samples = pd.read_csv(config["samples"], sep="\t")

wildcard_constraints:
    ext=r"vcf(\.gz)?$"


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
    return sex[0]

def get_analysis_type(wildcards):
    analysis_type = _get_sample_row(wildcards)["type"]
    assert len(analysis_type) == 1
    return analysis_type[0]

def get_track(wildcards):
    track = _get_sample_row(wildcards)["track"]
    assert len(track) == 1
    return track[0]

def get_vcf_samples(wildcards):
    vcf_filename = _get_sample_row(wildcards)["vcf"]
    assert len(vcf_filename) == 1
    vcf = cyvcf2.VCF(vcf_filename[0])
    samples = vcf.samples
    assert len(samples) == 1
    return samples[0]

def get_rank_model(wildcards):
    sample_track = get_track(wildcards)
    # TODO: set this up for cancer and structural variants
    version = config["genmod"]["rd_rank_model_version"]
    return f"rank_model/rd_rank_model_v{version}.ini"

def get_result_files():
    infiles = []
    outfiles = []
    load_configs = []
    for s, p in zip(samples["sample"], samples["ped"]):
        infiles.append(f"annotation/{s}/{s}.decomposed.vep.annovar.genmod.vcf.gz")
        outfiles.append(f"results/{s}/{s}.scout-annotated.vcf.gz")

        infiles.append(f"annotation/{s}/{s}.decomposed.vep.annovar.genmod.vcf.gz.tbi")
        outfiles.append(f"results/{s}/{s}.scout-annotated.vcf.gz.tbi")

        if isinstance(p, str) and len(p) > 0:
            infiles.append(p)
        else:
            infiles.append(f"mock_ped/{s}.ped")
        outfiles.append(f"results/{s}/{s}.ped")

        load_configs.append(f"results/{s}/{s}.load_config.yaml")

    return infiles, outfiles, load_configs

infiles, outfiles, load_configs = get_result_files()

# for inf, outf in zip(infiles, outfiles):
#     print(inf, ":", outf)