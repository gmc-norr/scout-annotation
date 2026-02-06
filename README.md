# Scout annotation workflow

Snakemake workflow for annotating VCFs prior loading into Scout.

## Setup

The package is managed using uv, and this can be used to install it as a tool, or just run it directly in the project directory:

```bash
uv sync
uv run scout-annotation
# or
uv tool install --python 3.12 .
```

> [!NOTE]
> Due to a [current limitation in uv](https://github.com/astral-sh/uv/issues/11624), the Python version restriction in `pyproject.toml` isn't fully respected.
> Therefore the version must be passed in when installing using `uv tool install`.

It can also be installed using pip. Just activate the environment you want to use and run:

```bash
python -m pip install .
```

Currently only Python 3.12 is supported.

## Running 

```bash
scout-annotation --help
# Usage: scout-annotation [OPTIONS] COMMAND [ARGS]...
#
# Options:
#   -c, --config TEXT               config file used for overwriting defaults
#   --cores INTEGER                 number of cores available for snakemake
#   --use-apptainer, --use-singularity
#                                   use apptainer as executor
#   --apptainer-args, --singularity-args TEXT
#                                   arguments for apptainer
#   --apptainer-prefix, --singularity-prefix TEXT
#                                   path to cached apptainer containers
#   --loglevel [DEBUG|INFO|WARNING|ERROR]
#                                   set logging level
#   --version                       Show the version and exit.
#   -h, --help                      Show this message and exit.
#
# Commands:
#   batch   Annotate a batch of samples.
#   panels  List available gene panels
#   single  Annotate a single sample.
#   trio    Annotate a trio of samples
```

## Testing

This package is managed by uv, and tests are implemented using pytest.

```bash
uv run pytest
```

Running the above will run both unit tests and an integration test.
Currently, the integration test is only expected to run properly on one of the RV servers where `/storage` is mounted.
If `/storage` is not found, the integration tests will be skipped.

## Rule graph

```mermaid
flowchart TB
    1[rename_samples] -->
        2[bcftools_reheader] -->
        3[decompose] -->
        4[normalize] -->
        5[vt_sort] -->
        6[vt_uniq] -->
        7[fix_vcf_af] --> 30{Part of trio?}

    30 -- yes -->
        31[merge_trio] -->  9

    31 --> 70[peddy]
    9 --> 70
    10 --> 70

    30 -- no --> 9

    9[bgzip] -->
        11[vep] -->
        12[vcfanno] -->
        14[most_severe_consequence]

    9 --> 10[tabix] --> 11[vep]

    2 -->
        29{PED file given?} -- no -->
        8[mock_ped]

    8 --> 11
    8 --> 21
    8 --> 23

    13[vcfanno_config] --> 12
    14[most_severe_consequence] -->
        15{Panels and/or filters?}
    15 -- only filter --> 17[vcf_filtering]
    15 -- panels --> 18[panel_filtering]
    18 --> 19{Filters too?}
    19 -- yes --> 17[vcf_filtering]
    19 -- no --> 20{Any variants left?}
    15 -- none --> 20

    17 --> 20

    20 -- no --> 21[copy_results]

    20 -- yes -->
        16[genmod_annotate] -->
        22[genmod_models] -->
        23[genmod_score] -->
        24[genmod_compound] -->
        26[bgzip] -->
        21 -->
        40[sample_config] -->
        41[family_config] -->
        28[all]

    26 --> 27[tabix] --> 21

    50{BAM files given?} -- yes -->
        51[link_bam] -->
        28

    51 --> 40
    70 --> 21

    30 -- yes -->
        80[madeline2_data] -->
        81[madeline2_svg] -->
        21
    29 -- yes --> 80

    60[genmod_rank_model] --> 23
```
