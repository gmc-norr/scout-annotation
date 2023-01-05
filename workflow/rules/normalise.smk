rule decompose:
    input:
        vcf=get_vcf_file
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.vcf")
    log: "decompose/{sample}/{sample}.decomposed.log"
    container: "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt decompose {input.vcf} -o {output.vcf} 2> {log}
        """
