rule merge_trio:
    input:
        vcf=get_family_vcfs,
        idx=get_family_vcf_index,
    output:
        vcf=temp("merge_trio/{family}_trio.vcf"),
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        "bcftools merge -o {output.vcf} {input.vcf}"
