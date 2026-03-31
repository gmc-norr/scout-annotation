import re
import pysam
import logging
from typing import Dict, List, Tuple, Optional, Set

logger = logging.getLogger(__name__)

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


def expected_decomposed_snvs(
    old: Optional[Tuple[str, int, str, str]],
) -> Optional[List[Tuple[str, int, str, str]]]:
    """
    For a same-length MNV in OLD_CLUMPED, return the expected decomposed SNVs.

    Example:
      chr1:1111:ACGT/GGGC
    becomes:
      [
        ("chr1", 1111, "A", "G"),
        ("chr1", 1112, "C", "G"),
        ("chr1", 1114, "T", "C"),
      ]
    """
    if old is None:
        return None

    chrom, pos, ref, alt = old

    if len(ref) != len(alt):
        return None

    expected = []
    for i, (r, a) in enumerate(zip(ref, alt)):
        if r != a:
            expected.append((chrom, pos + i, r, a))

    return expected


def rec_as_simple_snv(
    rec: pysam.VariantRecord,
) -> Optional[Tuple[str, int, str, str]]:
    """Return (chrom, pos, ref, alt) if record is a simple biallelic SNV, else None."""
    if len(rec.alleles) != 2:
        return None

    ref, alt = rec.alleles
    if len(ref) != 1 or len(alt) != 1:
        return None

    return rec.contig, rec.pos, ref, alt


def copy_same_info(
    member_recs: List[pysam.VariantRecord],
    new_rec: pysam.VariantRecord,
    header: pysam.VariantHeader,
) -> None:
    first = member_recs[0]

    # copy INFO from first record except OLD_CLUMPED and CSQ
    for key in first.info.keys():
        if key in {"OLD_CLUMPED", "CSQ"}:
            continue
        _set_info_value(new_rec, key, first.info[key], header)


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
    if old is None:
        raise ValueError(f"Cannot make merged record: invalid OLD_CLUMPED={old_txt!r}")

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


def _collect_old_clumped_status(
    vcf_in: str,
) -> Set[Tuple[str, int, str, str]]:
    """
    collect which OLD_CLUMPED groups have the full expected decomposed SNVs.
    """
    groups: Dict[Tuple[str, int, str, str], dict] = {}

    with pysam.VariantFile(vcf_in, "r") as in_vcf:
        for record in in_vcf:
            old_txt = record.info.get("OLD_CLUMPED")
            old = parse_old_clumped(old_txt)
            if old is None:
                continue

            if old not in groups:
                expected = expected_decomposed_snvs(old)
                groups[old] = {
                    "old_txt": old_txt,
                    "expected": expected,
                    "expected_set": set(expected) if expected is not None else None,
                    "observed": set(),
                    "invalid": expected is None,
                }

            group = groups[old]
            snv = rec_as_simple_snv(record)

            if snv is None:
                group["invalid"] = True
                continue

            if group["expected_set"] is None:
                group["invalid"] = True
                continue

            if snv not in group["expected_set"]:
                group["invalid"] = True
                continue

            group["observed"].add(snv)

    complete_groups: Set[Tuple[str, int, str, str]] = set()

    for old, group in groups.items():
        expected = group["expected"]
        observed = sorted(group["observed"], key=lambda x: x[1])

        if not group["invalid"] and group["expected_set"] == group["observed"]:
            complete_groups.add(old)
            chrom, pos, ref, alt = old
            logger.info(
                "Found complete OLD_CLUMPED group %r -> %s:%d %s>%s from %d decomposed SNVs",
                group["old_txt"],
                chrom,
                pos,
                ref,
                alt,
                len(observed),
            )
        else:
            logger.warning(
                "Incomplete or invalid OLD_CLUMPED group %r; expected=%s observed=%s. Will write original records unchanged.",
                group["old_txt"],
                expected,
                observed,
            )

    return complete_groups


def undecompose(vcf_in: str, vcf_out: str) -> None:
    """
    Read VCF in two passes and merge OLD_CLUMPED groups only when the full
    expected decomposed SNV set is present somewhere in the file.
    """
    complete_groups = _collect_old_clumped_status(vcf_in)
    written_groups: Set[Tuple[str, int, str, str]] = set()

    with pysam.VariantFile(vcf_in, "r") as in_vcf:
        with pysam.VariantFile(vcf_out, "w", header=in_vcf.header) as out_vcf:
            for record in in_vcf:
                old_txt = record.info.get("OLD_CLUMPED")
                old = parse_old_clumped(old_txt)

                out_vcf.write(record)

                if old is None:
                    continue

                if old not in complete_groups:
                    continue

                # complete group: write merged once, skip the rest
                if old in written_groups:
                    continue

                chrom, pos, ref, alt = old
                logger.info(
                    "Writing merged OLD_CLUMPED group %r as %s:%d %s>%s",
                    old_txt,
                    chrom,
                    pos,
                    ref,
                    alt,
                )
                out_vcf.write(make_new_rec([record], out_vcf))
                written_groups.add(old)


if __name__ == "__main__":
    try:
        logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])
        log = snakemake.log_fmt_shell(stdout=False, stderr=True)
        undecompose(
            snakemake.input.vcf,
            snakemake.output.vcf,
        )
    except NameError:
        import argparse

        logging.basicConfig(level=logging.INFO)
        p = argparse.ArgumentParser()
        p.add_argument("-i", "--input-vcf", required=True)
        p.add_argument("-o", "--output-vcf", required=True)
        a = p.parse_args()

        undecompose(a.input_vcf, a.output_vcf)
