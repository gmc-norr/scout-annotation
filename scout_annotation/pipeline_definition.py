from collections import defaultdict
from pathlib import Path
from pydantic import BaseModel, model_validator
from typing import TextIO
import yaml

from scout_annotation.path import WildcardPath


class FileModel(BaseModel):
    path: str | Path
    required: bool = True


class PipelineFiles(BaseModel):
    bam: FileModel | None = None
    bai: FileModel | None = None
    snv_vcf: FileModel
    report: FileModel | None = None
    cnv_report: FileModel | None = None
    hrd: FileModel | None = None
    tmb: FileModel | None = None
    msi: FileModel | None = None
    ped: FileModel | None = None

    @model_validator(mode="after")
    def check_bai_if_bam(self) -> "PipelineFiles":
        if self.bam and not self.bai:
            raise ValueError("bai is required if bam is supplied")
        return self

class PipelineDefinition(BaseModel):
    name: str
    version: str
    track: str
    owner: str
    files: PipelineFiles

    def resolve_paths(self, path: str | Path) -> dict[str, PipelineFiles]:
        samples = defaultdict(dict)
        for ftype, fd in self.files.model_dump(exclude_none=True).items():
            p = WildcardPath(path / fd["path"])
            fpaths = p.expand()
            if not fpaths:
                if fd["required"]:
                    raise ValueError(f"required file not found: {ftype}")
                continue
            for fpath in p.expand():
                if "sample" not in fpath[1]:
                    continue
                sample = fpath[1]["sample"]
                samples[sample][ftype] = {"path": fpath[0]}
        plfiles = {}
        for sample, d in samples.items():
            plfiles[sample] = PipelineFiles(**d)
        return plfiles

class PipelineDefinitions(BaseModel):
    pipelines: list[PipelineDefinition]

    def find(self, name: str, version: str) -> PipelineDefinition | None:
        for pd in self.pipelines:
            if pd.name == name and pd.version == version:
                return pd
        return None


def parse(file: TextIO) -> PipelineDefinitions:
    defs = yaml.safe_load(file)
    return PipelineDefinitions(**defs)


def read(path: str | Path) -> PipelineDefinitions:
    with open(path) as f:
        return parse(f)
