from collections import defaultdict
import pathlib
import cyvcf2
import pandas as pd
from typing import AnyStr, ByteString, Dict, List

CSQ_DEF = None


def get_panel_genes(panel_file: AnyStr) -> List[ByteString]:
    d = pd.read_csv(panel_file, sep=r"[,;\t]", engine="python")
    return d["hgnc_symbol"].values


def get_vep_gene(variant: cyvcf2.Variant) -> List[AnyStr]:
    return variant.INFO.get("CSQ").split("|")[CSQ_DEF.index("SYMBOL")]


def filter_vcf(vcf_filename: AnyStr, panels: Dict, output_filename: AnyStr, hard_filter: bool=True) -> cyvcf2.VCF:
    global CSQ_DEF
    vcf = cyvcf2.VCF(vcf_filename)
    CSQ_DEF = vcf.get_header_type("CSQ").get("Description").split()[-1].split("|")

    vcf.add_filter_to_header(dict(
        ID="PANEL",
        Description="Variant not in panel gene"
    ))

    vcf.add_info_to_header(dict(
        ID="PANEL",
        Number=".",
        Type="String",
        Description="Gene panel in which the gene was found, if any"
    ))

    vcf.add_to_header(f"##panels={','.join(panels.keys())}")

    panel_union = set()
    gene_to_panel = defaultdict(list)
    for pname, p in panels.items():
        panel_union.update(set(p))
        for g in p:
            gene_to_panel[g].append(pname)
    
    writer = cyvcf2.Writer(output_filename, vcf, "w")
    writer.write_header()

    for variant in vcf:
        variant_gene = get_vep_gene(variant)
        variant.INFO["PANEL"] = ",".join(gene_to_panel[variant_gene])
        if variant_gene not in panel_union:
            if variant.FILTER is None:
                variant.FILTER = ["PANEL"]
            else:
                variant.FILTER = variant.FILTERS + ["PANEL"]
            if not hard_filter:
                writer.write_record(variant)
        else:
            writer.write_record(variant)

    return vcf


def main():
    vcf_filename = snakemake.input["vcf"]
    panel_files = snakemake.input["panels"]

    output_filename = snakemake.output["vcf"]

    hard_filter = snakemake.params["hard_filter"]

    panels = {}
    for p in panel_files:
        p = pathlib.Path(p)
        panels[p.stem] = get_panel_genes(p) 

    filter_vcf(vcf_filename, panels, output_filename, hard_filter)


if __name__ == "__main__":
    main()
