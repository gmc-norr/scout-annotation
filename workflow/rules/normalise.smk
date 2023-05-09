rule sample_name_translation:
    output:
        samples=temp("decompose/{sample}/{sample}.samples.txt"),
    params:
        old_name=lambda wc: get_vcf_samples(get_vcf_file(wc)),
        new_name=lambda wc: wc.sample,
    run:
        with open(output.samples, "w") as f:
            print(f"{params.old_name} {params.new_name}", file=f)


rule bcftools_reheader:
    input:
        vcf=get_vcf_file,
        samples="decompose/{sample}/{sample}.samples.txt",
    output:
        vcf=temp("decompose/{sample}/{sample}.renamed.vcf"),
    log:
        "decompose/{sample}/{sample}.renamed.log",
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        """
        bcftools reheader -s {input.samples} -o {output.vcf} {input.vcf} 2> {log}
        """


rule decompose:
    input:
        vcf=get_reheadered_vcf_file,
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.vcf"),
    log:
        "decompose/{sample}/{sample}.decomposed.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        (vt decompose -s {input.vcf} | vt decompose_blocksub -o {output.vcf} -) 2> {log}
        """


rule normalize:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.vcf",
        fasta=config.get("reference", {}).get("fasta"),
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.vcf"),
    log:
        "decompose/{sample}/{sample}.decomposed.normalized.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt normalize -n -r {input.fasta} {input.vcf} -o {output.vcf} 2> {log}
        """


rule vt_sort:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.normalized.vcf",
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.sort.vcf"),
    log:
        "decompose/{sample}/{sample}.decomposed.normalized.sort.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt sort -o {output.vcf} {input.vcf} 2> {log}
        """


rule vt_uniq:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.normalized.sort.vcf",
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.uniq.vcf"),
    log:
        "decompose/{sample}/{sample}.decomposed.normalized.uniq.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt uniq -o {output.vcf} {input.vcf} 2> {log}
        """


rule fix_vcf_af:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.normalized.uniq.vcf",
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.normalized.uniq.fix-af.vcf"),
    log:
        "decompose/{sample}/{sample}.decomposed.fix-af.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    script:
        "../scripts/fix_vcf_af.py"
