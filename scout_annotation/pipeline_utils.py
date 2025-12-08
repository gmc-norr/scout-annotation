import json
from pathlib import Path
import re
from typing import Tuple
import yaml


def detect_pipeline(path: str | Path) -> (str, str):
    snakemake_dir = Path(path) / ".snakemake"
    nextflow_dir = Path(path) / ".nextflow"
    if snakemake_dir.exists() and snakemake_dir.is_dir():
        return _identify_snakemake_pipeline(path)
    elif nextflow_dir.exists() and nextflow_dir.is_dir():
        return _identify_nextflow_pipeline(path)
    raise ValueError("could not identify path as snakemake or nextflow pipeline results")


def is_version(s: str) -> bool:
    """Check if a string is a three-digit version number.
    Surrounding whitespace is allowed, as well as a v-prefix.
    """
    parts = s.lower().strip().lstrip("v").split(".")
    return len(parts) == 3 and all(x.isdigit() for x in parts)


def _identify_snakemake_pipeline(path: str | Path) -> Tuple[str, str]:
    name, version = general_report_pipeline_version(path)
    if name is not None and version is not None:
        return name, version
    name, version = multiqc_pipeline_version(path)
    if name is not None and version is not None:
        return name, version
    else:
        raise ValueError("could not find snakemake pipeline version")


def general_report_pipeline_version(path: str | Path) -> Tuple[str | None, str | None]:
    general_report_json = Path(path).glob("**/*.general.json")
    for json_path in general_report_json:
        with open(json_path) as f:
            d = json.load(f)
            name = d.get("pipeline", {}).get("name")
            version = d.get("pipeline", {}).get("version")
            if name is None or version is None:
                continue
            return name, version
    else:
        return None, None


def multiqc_pipeline_version(path: str | Path) -> Tuple[str | None, str | None]:
    version_d = Path(path) / "results" / "versions"
    if not version_d.exists():
        raise ValueError("potentially unsupported snakemake pipeline")
    for yaml_file in version_d.glob("**/*.yaml"):
        with open(yaml_file) as f:
            d = yaml.safe_load(f)
            if d and len(d) == 1:
                name = list(d.keys())[0]
                version = list(d.values())[0]
                if not is_version(version):
                    continue
                return name, version
                break
    else:
        return None, None


def _identify_nextflow_pipeline(path: str | Path) -> Tuple[str, str]:
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    logfile = Path(path) / ".nextflow.log"
    if not logfile.exists():
        raise ValueError("could not find nextflow log")
    with open(logfile) as f:
        for line in f:
            line = ansi_escape.sub("", line).strip().split()
            if len(line) != 2 or not is_version(line[1]):
                continue
            return (line[0], line[1])
    raise ValueError("could not find nextflow pipeline version")
