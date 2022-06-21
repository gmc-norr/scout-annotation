import pandas as pd
from pathlib import Path
import yaml

configfile: "config/config.yaml"
with open(config["resources"]) as f:
    resources = yaml.load(f, Loader=yaml.FullLoader)

samples = pd.read_csv(config["samples"], sep="\t")

wildcard_constraints:
    ext=r"vcf(\.gz)?$"

def get_vcf_file(sample):
    return samples[samples["sample"] == sample]["vcf"].values

def get_ped(wildcards):
    return samples[samples["sample"] == wildcards.sample]["ped"].values