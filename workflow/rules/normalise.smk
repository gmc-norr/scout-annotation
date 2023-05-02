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

rule normalize:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.vcf",
        fasta=config.get("reference", {}).get("fasta"),
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.vcf"),
    log: "decompose/{sample}/{sample}.decomposed.normalized.log",
    container: "docker://hydragenetics/vt:2015.11.10",
    shell:
        """
        vt normalize -n -r {input.fasta} {input.vcf} -o {output.vcf} 2> {log}
        """

rule vt_uniq:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.normalized.vcf",
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.uniq.vcf"),
    log: "decompose/{sample}/{sample}.decomposed.normalized.uniq.log",
    container: "docker://hydragenetics/vt:2015.11.10",
    shell:
        """
        vt uniq {input.vcf} -o {output.vcf} 2> {log}
        """

rule fix_vcf_af:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.normalized.uniq.vcf"
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.uniq.fix-af.vcf")
    log: "decompose/{sample}/{sample}.decomposed.fix-af.log"
    container: "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    script: "../scripts/fix_vcf_af.py"
