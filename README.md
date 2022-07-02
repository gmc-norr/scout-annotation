# Scout annotation workflow

Snakemake workflow for annotating VCFs prior loading into Scout.

## Setup

```bash
python -m virtualenv -p python3.9 venv
. ./venv/bin/activate
```

## Integration test

On vs478:

```bash
cd .test
snakemake \
    -s ../workflow/Snakefile \
    -c 1 \
    --use-singularity \
    --singularity-args "--bind /storage" \
    --use-conda \
    --configfile config.yaml
```
