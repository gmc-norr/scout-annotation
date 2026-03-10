import re
import pysam
from typing import List, Tuple, Optional

old_re = re.compile(
    r"^(?P<chrom>.+):(?P<pos>\d+):(?P<ref>[ACGTN]+)\/(?P<alt>[ACGTN]+)$"
)


def parse_old_clumped(v: str) -> Optional[Tuple[str, int, str, str]]:
    """Parse OLD_CLUMPED like chr1:1111111:AA/CC"""
    if v is None:
        return None

    m = old_re.match(str(v))
    if not m:
        return None

    chrom = m.group("chrom")
    pos = int(m.group("pos"))
    ref = m.group("ref")
    alt = m.group("alt")

    return chrom, pos, ref, alt


def copy_same_info(
    member_recs: List[pysam.VariantRecord],
    new_rec: pysam.VariantRecord,
    header: pysam.VariantHeader,
) -> None:
    first = member_recs[0]

    # copy INFO from first record except OLD_CLUMPED, CALLERS and CSQ
    for key in first.info.keys():
        if key in {"OLD_CLUMPED", "CALLERS", "CSQ"}:
            continue
        _set_info_value(new_rec, key, first.info[key], header)

    callers = []
    seen = set()
    for rec in member_recs:
        for c in _as_list(rec.info.get("CALLERS")):
            if c not in seen:
                seen.add(c)
                callers.append(c)
    _set_info_value(new_rec, "CALLERS", callers, header)


def copy_fmt(
    member_recs: List[pysam.VariantRecord],
    new_rec: pysam.VariantRecord,
) -> None:
    first = member_recs[0]

    for sample in first.samples:
        src = first.samples[sample]
        dst = new_rec.samples[sample]

        for fmt_key in src.keys():
            dst[fmt_key] = src[fmt_key]


def make_new_rec(
    recs: List[pysam.VariantRecord],
    out_vcf: pysam.VariantFile,
) -> pysam.VariantRecord:
    first = recs[0]

    old_txt = first.info.get("OLD_CLUMPED")
    old = parse_old_clumped(old_txt)

    chrom, pos, ref, alt = old

    new_rec = out_vcf.new_record(
        contig=chrom,
        start=pos - 1,
        stop=(pos - 1) + len(ref),
        alleles=(ref, alt),
        filter=list(first.filter.keys()),
    )

    copy_same_info(recs, new_rec, out_vcf.header)
    copy_fmt(recs, new_rec)

    return new_rec


def _as_list(x):
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


def _set_info_value(
    rec: pysam.VariantRecord,
    key: str,
    value,
    header: pysam.VariantHeader,
) -> None:
    if value is None:
        return

    if key not in header.info:
        return

    number = header.info[key].number

    if number == 0:
        if value:
            rec.info[key] = True
        return

    if number == 1:
        if isinstance(value, (list, tuple)):
            if len(value) == 0:
                return
            rec.info[key] = value[0]
        else:
            rec.info[key] = value
        return

    if isinstance(value, list):
        rec.info[key] = tuple(value)
    else:
        rec.info[key] = value


def _flush_buffer(
    buffer: List[pysam.VariantRecord],
    out_vcf: pysam.VariantFile,
) -> None:
    if not buffer:
        return
    out_vcf.write(make_new_rec(buffer, out_vcf))
    buffer.clear()


def undecompose(vcf_in: str, vcf_out: str) -> None:
    """Read vcf and merge adjacent decomposed records back using OLD_CLUMPED."""
    with pysam.VariantFile(vcf_in, "r") as in_vcf:
        with pysam.VariantFile(vcf_out, "w", header=in_vcf.header) as out_vcf:
            buffer: List[pysam.VariantRecord] = []
            prev_old_clumped = None

            for record in in_vcf:
                old_val = record.info.get("OLD_CLUMPED")
                old_clumped = parse_old_clumped(old_val)

                if old_clumped is None:
                    _flush_buffer(buffer, out_vcf)
                    prev_old_clumped = None
                    out_vcf.write(record)
                    continue

                if prev_old_clumped is None or old_clumped == prev_old_clumped:
                    buffer.append(record)
                    prev_old_clumped = old_clumped
                    continue

                _flush_buffer(buffer, out_vcf)
                buffer.append(record)
                prev_old_clumped = old_clumped

            _flush_buffer(buffer, out_vcf)


if __name__ == "__main__":
    try:
        log = snakemake.log_fmt_shell(stdout=False, stderr=True)
        undecompose(
            snakemake.input.vcf,
            snakemake.output.vcf,
        )
    except NameError:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("-i", "--input-vcf", required=True)
        p.add_argument("-o", "--output-vcf", required=True)
        a = p.parse_args()

        undecompose(a.input_vcf, a.output_vcf)
