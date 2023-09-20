import cyvcf2
import numpy as np
import pandas as pd
from pathlib import Path
from snakemake.utils import validate
import sys
import yaml

if not workflow.overwrite_configfiles:
    print("error: ", file=sys.stderr)

validate(config, "../schema/config.schema.yaml")

# Resolve relative paths to other config files
paths_to_check = zip(
    [
        config.get("panel_filtering", {}).get("panel_directory"),
        *config.get("vcf_filter", {}).values(),
        config.get("peddy", {}).get("sites"),
        config.get("resources"),
    ],
    [
        ("panel_filtering", "panel_directory"),
        *[("vcf_filter", x) for x in config.get("vcf_filter", {}).keys()],
        ("peddy", "sites"),
        ("resources",),
    ],
)

for filepath, config_key in paths_to_check:
    if filepath is None:
        continue

    filepath = Path(filepath)
    if filepath.is_absolute() or filepath.exists():
        continue

    # Try to resolve relative to the location of any of the config files,
    # starting from the last one supplied
    path_found = False
    for config_path in workflow.configfiles[::-1]:
        new_path = (Path(config_path).parent / filepath).resolve()
        if new_path.exists():
            path_found = True
            config_item = config
            for ck in config_key[:-1]:
                config_item = config_item.get(ck, {})
            config_item[config_key[-1]] = str(new_path)
            break
    if not path_found:
        raise IOError(f"file or directory not found for filter '{key}': {filepath}")

with open(config["resources"]) as f:
    resources = yaml.load(f, Loader=yaml.FullLoader)
validate(resources, "../schema/resources.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", comment="#")
validate(samples, "../schema/samples.schema.yaml")

# add empty family column if it doesn't exist
if samples.get("family", pd.Series(dtype="str")).empty:
    samples["family"] = None

# if family ID is missing, set it to the same as the sample ID
samples["family"] = np.where(
    samples["family"].isnull(), samples["sample"], samples["family"]
)


wildcard_constraints:
    ext=r"vcf(\.gz)?$",
    family=r"[a-zA-Z][-a-zA-Z0-9]+",
    sample=r"[a-zA-Z][-a-zA-Z0-9]+",
    track=r"(rare_disease|cancer)",


def _get_sample_row(wildcards):
    sample_row = samples[samples["sample"] == wildcards.sample]
    assert (
        len(sample_row) == 1
    ), f"expected 1 sample, found {len(sample_row)} for {wildcards.sample}"
    return sample_row


def _get_family_rows(wildcards):
    family_rows = samples[samples["family"] == wildcards.family]
    assert len(family_rows) > 0, f"no samples found for family: {wildcards.family}"
    return family_rows


def get_annotated_vcf(wildcards):
    snv_filter = get_family_snv_filter_tag(wildcards)
    sample_panels = get_family_panels(wildcards)

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
        vcf_filename = "annotation/{family}/{family}.annotated".format(**wildcards)
        if len(sample_panels) > 0:
            vcf_filename = (
                "annotation/{family}/{family}.annotated.panel_filtered".format(
                    **wildcards
                )
            )
        vcf_filtering = checkpoints.vcf_filtering.get(
            file=vcf_filename,
            tag=get_family_snv_filter_tag(wildcards),
        )
        with vcf_filtering.output["vcf"].open() as f:
            n_variants = sum(1 for line in f if not line.startswith("#"))
        if n_variants == 0:
            # No variants left for ranking
            return f'{vcf_filtering.output["vcf"]}.gz'

    # genmod-annotated VCF
    return "annotation/{family}/{family}.annotated.genmod.vcf.gz".format(**wildcards)


def get_annotated_vcf_index(wildcards):
    return f"{get_annotated_vcf(wildcards)}.tbi"


def get_filtered_vcf(wildcards):
    snv_filter = get_family_snv_filter_tag(wildcards)
    sample_panels = get_family_panels(wildcards)

    if snv_filter is None and len(sample_panels) == 0:
        # no filtering at all, return annotated VCF
        return "annotation/{family}/{family}.annotated.vcf"

    if snv_filter is None:
        # only panel filtering
        return "annotation/{family}/{family}.annotated.panel_filtered.vcf".format(
            **wildcards
        )

    # both panel and snv filtering
    return "annotation/{family}/{family}.annotated.panel_filtered.filter-{snv_filter}.vcf".format(
        **wildcards, snv_filter=snv_filter
    )


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


def get_case_owner(wildcards):
    case_owner = _get_family_rows(wildcards)["owner"].unique()
    assert len(case_owner) == 1, "all samples in a family must have the same owner"
    return case_owner[0]


def get_family_id(wildcards):
    sample_row = _get_sample_row(wildcards)
    family_id = sample_row["family"].values[0]
    if not pd.isnull(family_id):
        return family_id
    return sample_row["sample"].values[0]


def get_panel_dict():
    panels = {}
    panel_path = Path(config.get("panel_filtering", {}).get("panel_directory"))
    if not panel_path.exists():
        raise FileNotFoundError(f"directory not found: {panel_path}")
    for p in panel_path.glob("*.tsv"):
        panels[p.stem] = p
    return panels


def get_peddy_file(wildcards, filetype):
    samples = get_family_samples(wildcards)
    if len(samples) < 2:
        return []
    peddy_files = {
        "ped": "peddy/{family}/{family}.peddy.ped",
        "het_check": "peddy/{family}/{family}.het_check.csv",
        "ped_check": "peddy/{family}/{family}.ped_check.csv",
        "sex_check": "peddy/{family}/{family}.sex_check.csv",
        "html": "peddy/{family}/{family}.peddy.ped",
    }
    return peddy_files[filetype]


def get_preprocessed_vcf(wildcards):
    # if trio, get merged vcf, otherwise get sample vcf
    family_samples = samples[samples["family"] == wildcards["family"]]["sample"].values
    assert len(family_samples) == 1 or len(family_samples) == 3
    if len(family_samples) == 1:
        return "decompose/{sample}/{sample}.decomposed.normalized.uniq.fix-af.vcf.gz".format(
            sample=family_samples[0]
        )
    return "merge_trio/{family}_trio.vcf.gz".format(**wildcards)


def get_preprocessed_vcf_index(wildcards):
    vcf = get_preprocessed_vcf(wildcards)
    return f"{vcf}.tbi"


def get_family_panels(wildcards):
    panels = samples[samples["family"] == wildcards["family"]]["panels"].unique()
    assert len(panels) == 1, "all samples in a family need to have the same panels"
    if pd.isnull(panels[0]):
        return []
    return panels[0].split(",")


def get_family_snv_filter_tag(wildcards):
    filtering = samples[samples["family"] == wildcards.family]["filtering"].unique()
    assert len(filtering) == 1, "all samples in a family must use the same filters"
    if pd.isnull(filtering[0]) or filtering[0] == "":
        return None
    return filtering[0]


def get_panel_files(wildcards):
    panel_dict = get_panel_dict()
    try:
        panel_files = [panel_dict[p] for p in get_family_panels(wildcards) if p != ""]
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
    vcf_sample_name = get_vcf_samples(get_vcf_file(wildcards))
    if vcf_sample_name != wildcards.sample:
        return rules.bcftools_reheader.output.vcf
    return get_vcf_file(wildcards)


def get_family_ped(wildcards):
    peds = _get_family_rows(wildcards)["ped"].unique()
    assert len(peds) == 1, "samples within a family must have the same ped file"
    if not pd.isnull(peds[0]):
        return peds[0]
    return "mock_ped/{family}.ped".format(**wildcards)


def get_family_samples(wildcards):
    return _get_family_rows(wildcards)["sample"]


def get_family_vcfs(wildcards):
    """Get the sample VCFs from a family ID"""
    # TODO: support family structures other than trios
    trio = samples[samples["family"] == wildcards.family]["sample"].values
    assert len(trio) == 3, f"found {len(trio)} samples, expected 3"
    return [
        f"decompose/{sample}/{sample}.decomposed.normalized.uniq.fix-af.vcf.gz"
        for sample in trio
    ]


def get_family_vcf_index(wildcards):
    vcfs = get_family_vcfs(wildcards)
    return [f"{vcf}.tbi" for vcf in vcfs]


def get_sample_sex(wildcards):
    sex = _get_sample_row(wildcards)["sex"]
    assert len(sex) == 1
    return sex.values[0]


def get_analysis_type(wildcards):
    analysis_type = _get_family_rows(wildcards)["type"].unique()
    assert (
        len(analysis_type) == 1
    ), "all samples within a family must have the same analysis type"
    return analysis_type[0]


def get_track(wildcards):
    # track = _get_sample_row(wildcards)["track"]
    tracks = samples[samples["family"] == wildcards.family]["track"].unique()
    assert len(tracks) == 1, "all samples in a family need to belong to the same track"
    return tracks[0]


def get_vcf_samples(vcf_filename):
    vcf = cyvcf2.VCF(vcf_filename)
    assert len(vcf.samples) > 0, f"no samples found in VCF: {vcf_filename}"
    return vcf.samples


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
            expressions.append(
                " ".join([x.strip() for x in filter["expression"].splitlines()])
            )
    return "(" + ") and (".join(expressions) + ")"


def get_output_files():
    outfiles = []
    load_configs = []
    outdir = config.get("output_directory", "results")
    for f in set(samples["family"].values):
        outfiles.append(f"{outdir}/{f}/{f}.scout-annotated.vcf.gz")
        outfiles.append(f"{outdir}/{f}/{f}.scout-annotated.vcf.gz.tbi")
        outfiles.append(f"{outdir}/{f}/{f}.ped")
        family_wildcard = snakemake.io.Namedlist(fromdict={"family": f})
        family_samples = get_family_samples(family_wildcard)
        assert len(family_samples) == 1 or len(family_samples) == 3
        for s in family_samples:
            sample_wildcard = snakemake.io.Namedlist(fromdict={"sample": s})
            if isinstance(get_bam_file(sample_wildcard), Path):
                outfiles.append(f"{outdir}/{s}/{s}.bam")
                outfiles.append(f"{outdir}/{s}/{s}.bam.bai")
        if len(family_samples) == 3:
            outfiles.append(f"{outdir}/{f}/{f}.peddy.ped")
        load_configs.append(f"{outdir}/{f}/{f}.load_config.yaml")
    return outfiles, load_configs


outfiles, load_configs = get_output_files()
