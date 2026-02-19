from pathlib import Path
from snakemake.shell import shell


def common_prefix(paths):
    paths = sorted(list(paths))

    n = len(paths)
    if n == 0:
        return None
    if n == 1:
        return paths[0]

    prefix_len = min(len(paths[0]), len(paths[-1]))

    i = 0
    while i < prefix_len and paths[0][i] == paths[-1][i]:
        i += 1

    return paths[0][0:i].rstrip(".-_")


vcf = Path(snakemake.input.vcf)
ped = Path(snakemake.input.ped)

out_ped = Path(snakemake.output.ped)

prefix = f"--prefix {common_prefix(snakemake.output)}"

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

each = snakemake.params.get("each", None)
plot = snakemake.params.get("plot")
sites = snakemake.params.get("sites")

if sites:
    sites = f"--sites {sites}"
else:
    sites = ""

if each:
    each = f"--each {each}"
else:
    each = ""

if plot:
    plot = "--plot"
else:
    plot = ""

shell(
    """
    (python -m peddy \\
        --procs {snakemake.threads} \\
        {prefix} \\
        {sites} \\
        {each} \\
        {plot} \\
        {extra} \\
        {vcf} \\
        {ped}) {log}
    """
)
