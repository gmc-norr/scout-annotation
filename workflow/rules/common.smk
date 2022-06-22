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
    print(ped, wildcards)
    assert len(ped) == 1
    if not ped[0]:
        print(ped)
        return ped
    return "results/{sample}/{sample}.ped".format(**wildcards)