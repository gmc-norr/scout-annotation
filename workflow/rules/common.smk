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

def get_vcf_file(wildcards):
    vcf_filename = samples[samples["sample"] == wildcards.sample]["vcf"].values
    assert len(vcf_filename) == 1, "duplicate sample IDs found"
    return vcf_filename[0]

def get_ped(wildcards):
    ped = samples[samples["sample"] == wildcards.sample]["ped"].values
    assert len(ped) == 1
    if not ped[0]:
        return ped
    return "mock_ped/{sample}.ped".format(**wildcards)

def get_result_files():
    infiles = []
    outfiles = []
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
    return infiles, outfiles

infiles, outfiles = get_result_files()

# for inf, outf in zip(infiles, outfiles):
#     print(inf, ":", outf)