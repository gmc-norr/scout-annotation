# Scout annotation workflow

Snakemake workflow for annotating VCFs prior loading into Scout.

## Setup

### uv

The package is uv-compatible, and can be installed into a virtual environment in the repo by running:

```bash
git clone https://github.com/gmc-norr/scout-annotation
cd scout-annotation
uv sync
```

### pip

It can also be installed using pip. Just activate the environment you want to use and run:

```bash
git clone https://github.com/gmc-norr/scout-annotation
cd scout-annotation
python -m pip install .
```

For development, the dev dependencies are useful (requires pip >=25.1):

```bash
python -m pip install --group dev
```

Currently only Python 3.12 is supported.

## Example Usage

```bash
snakemake \
    -s scout-annotation/workflow/Snakefile \
    --config samples=samples_minimal.tsv output_directory=results_test \
    --configfiles config.yaml \
    --debug-dag \
    --cores 1
```

## Testing

The simplest way to run the tests is using uv:

```bash
uv run pytest
```

Running the above will run both unit tests and an integration test.
Currently, the integration test is only expected to run properly on one of the RV servers where `/storage` is mounted.
If `/storage` is not found, the integration tests will be skipped.
