import cyvcf2


SO_RANK = {
    "transcript_ablation": dict(so_id="SO:0001893", rank=1),
    "splice_acceptor_variant": dict(so_id="SO:0001574", rank=2),
    "splice_donor_variant": dict(so_id="SO:0001575", rank=3),
    "stop_gained": dict(so_id="SO:0001587", rank=4),
    "frameshift_variant": dict(so_id="SO:0001589", rank=5),
    "stop_lost": dict(so_id="SO:0001578", rank=6),
    "start_lost": dict(so_id="SO:0002012", rank=7),
    "transcript_amplification": dict(so_id="SO:0001889", rank=8),
    "inframe_insertion": dict(so_id="SO:0001821", rank=9),
    "inframe_deletion": dict(so_id="SO:0001822", rank=10),
    "missense_variant": dict(so_id="SO:0001583", rank=11),
    "protein_altering_variant": dict(so_id="SO:0001818", rank=12),
    "splice_region_variant": dict(so_id="SO:0001630", rank=13),
    "splice_donor_5th_base_variant": dict(so_id="SO:0001787", rank=14),
    "splice_donor_region_variant": dict(so_id="SO:0002170", rank=15),
    "splice_polypyrimidine_tract_variant": dict(so_id="SO:0002169", rank=16),
    "incomplete_terminal_codon_variant": dict(so_id="SO:0001626", rank=17),
    "start_retained_variant": dict(so_id="SO:0002019", rank=18),
    "stop_retained_variant": dict(so_id="SO:0001567", rank=19),
    "synonymous_variant": dict(so_id="SO:0001819", rank=20),
    "coding_sequence_variant": dict(so_id="SO:0001580", rank=21),
    "mature_miRNA_variant": dict(so_id="SO:0001620", rank=22),
    "5_prime_UTR_variant": dict(so_id="SO:0001623", rank=23),
    "3_prime_UTR_variant": dict(so_id="SO:0001624", rank=24),
    "non_coding_transcript_exon_variant": dict(so_id="SO:0001792", rank=25),
    "intron_variant": dict(so_id="SO:0001627", rank=26),
    "NMD_transcript_variant": dict(so_id="SO:0001621", rank=27),
    "non_coding_transcript_variant": dict(so_id="SO:0001619", rank=28),
    "upstream_gene_variant": dict(so_id="SO:0001631", rank=29),
    "downstream_gene_variant": dict(so_id="SO:0001632", rank=30),
    "TFBS_ablation": dict(so_id="SO:0001895", rank=31),
    "TFBS_amplification": dict(so_id="SO:0001892", rank=32),
    "TF_binding_site_variant": dict(so_id="SO:0001782", rank=33),
    "regulatory_region_ablation": dict(so_id="SO:0001894", rank=34),
    "regulatory_region_amplification": dict(so_id="SO:0001891", rank=35),
    "feature_elongation": dict(so_id="SO:0001907", rank=36),
    "regulatory_region_variant": dict(so_id="SO:0001566", rank=37),
    "feature_truncation": dict(so_id="SO:0001906", rank=38),
    "intergenic_variant": dict(so_id="SO:0001628", rank=39)
}


def main():
    vcf = cyvcf2.VCF(snakemake.input.vcf)

    vcf.add_info_to_header(dict(
        ID="MOST_SEVERE_CONSEQUENCE",
        Number=".",
        Type="String",
        Description="Most severe variant consequence according to https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences",
    ))

    consequence_index = vcf.get_header_type("CSQ")["Description"].split()[-1].split("|").index("Consequence")

    writer = cyvcf2.Writer(snakemake.output.vcf, vcf)
    writer.write_header()

    def get_most_severe_consequence(variant):
        consequences = variant.INFO["CSQ"].split("|")[consequence_index].split("&")
        if len(consequences) == 1:
            return consequences[0]
        most_severe = None
        most_severe_rank = 100
        for csq in consequences:
            if SO_RANK[csq]["rank"] < most_severe_rank:
                most_severe = csq
                most_severe_rank = SO_RANK[csq]["rank"]
        return most_severe

    for v in vcf:
        v.INFO["MOST_SEVERE_CONSEQUENCE"] = get_most_severe_consequence(v)
        writer.write_record(v)


if __name__ == "__main__":
    main()