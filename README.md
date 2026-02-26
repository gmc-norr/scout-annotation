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

## Example Usage 

```bash
snakemake -s scout-annotation/workflow/Snakefile --config samples=samples_minimal.tsv output_directory=results_test --configfiles config.yaml --debug-dag --cores 1

## Testing

This package is managed by uv, and tests are implemented using pytest.

```bash
uv run pytest
```

Running the above will run both unit tests and an integration test.
Currently, the integration test is only expected to run properly on one of the RV servers where `/storage` is mounted.
If `/storage` is not found, the integration tests will be skipped.

