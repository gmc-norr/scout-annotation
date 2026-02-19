import pandas as pd
from pathlib import Path
from snakemake.utils import validate
import sys

if not workflow.overwrite_configfiles:
    print("error: ", file=sys.stderr)

validate(config, "../schema/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", comment="#")
validate(samples, "../schema/samples.schema.yaml")

sample_dict = samples.set_index("sample").to_dict(orient="index")
family_dict = samples.set_index("family").to_dict(orient="index")

# Create output directories
results_dir = Path(config.get('results_dir_name', "scout_anno_results"))   
results_dir.mkdir(exist_ok=True)

annotation_dir = results_dir / "annotation"
annotation_dir.mkdir(exist_ok=True)
annotation_dir = str(annotation_dir)

coverage_dir = results_dir / "coverage"
coverage_dir.mkdir(exist_ok=True)

decompose_dir = results_dir / "decompose"
decompose_dir.mkdir(exist_ok=True)

output_vcfs = [f"{annotation_dir}/{family}/{family}.annotated.genmod.vcf" for family in family_dict.keys()]
output_d4s = []
for sample, sample_data in sample_dict.items():
    if sample_data.get("bam") is not None:
        output_d4s.append(f"{coverage_dir}/{sample_data['family']}/{sample}.coverage.d4")

def get_madeline2_svg(wildcards):
    samples = get_family_samples(wildcards)
    outdir = config.get("output_directory")
    if len(samples) < 2:
        return []
    return f"{outdir}/{{family}}/{{family}}.pedigree.svg"

def get_initial_vcf_file(wildcards):
    return family_dict[wildcards.family]['vcf']

def get_bam_file(wildcards):
    return family_dict[wildcards.family]['bam']


