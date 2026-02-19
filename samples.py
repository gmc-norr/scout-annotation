import hashlib
import pathlib
import time
from typing import Dict, List


def write_samples(samples: List[Dict], directory: pathlib.Path):
    if not directory.is_dir():
        directory.mkdir()

    sample_md5 = hashlib.md5(
        ",".join(d["sample"] for d in samples).encode("utf8")
    ).hexdigest()[:8]
    filename = pathlib.Path(
        directory, f"{time.strftime('%Y%m%d')}-{sample_md5}_samples.txt"
    )

    # Header
    cols = (
        "sample",
        "family",
        "owner",
        "sex",
        "type",
        "filtering",
        "track",
        "vcf",
        "bam",
        "ped",
        "panels",
        "msi_score",
        "hrd_score",
        "tmb_score"
    )
    with open(filename, "w") as f:
        print("\t".join(cols), file=f)
        for s in samples:
            print("\t".join(str(s[c]) if c in s else "" for c in cols), file=f)

    return filename
