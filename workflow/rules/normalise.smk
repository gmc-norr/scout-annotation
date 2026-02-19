rule sample_name_translation:
    output:
        samples=temp("decompose/{sample}/{sample}.samples.txt"),
    params:
        old_name=lambda wc: get_vcf_samples(get_vcf_file(wc))[0],
        new_name=lambda wc: wc.sample,
    localrule: True
    run:
        with open(output.samples, "w") as f:
            print(f"{params.old_name} {params.new_name}", file=f)


rule bcftools_reheader:
    input:
        vcf=get_initial_vcf_file,
    output:
        sample_conversion=temp(
            f"{decompose_dir}/{{family}}/{{family}}.sample_name_conversion.txt"
        ),
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.reheadered.vcf"),
    params:
        new_name=lambda wc: family_dict[wc.family]['sample'],
    log:
        f"{decompose_dir}/{{family}}/{{family}}.reheadered.log",
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        """
        echo "$(bcftools query -l file.bcf | head -1) {params.new_name}" >> {output.sample_conversion} &&
        bcftools reheader -s {output.sample_conversion} -o {output.vcf} {input.vcf} 2> {log}
        """

rule decompose:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.reheadered.vcf"
    output:
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.vcf"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        (vt decompose -s {input.vcf} | vt decompose_blocksub -o {output.vcf} -) 2> {log}
        """

rule normalize:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.decomposed.vcf",
        fasta=config.get("reference", {}).get("fasta"),
    output:
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.vcf"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt normalize -n -r {input.fasta} {input.vcf} -o {output.vcf} 2> {log}
        """


rule vt_sort:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.vcf",
    output:
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.sort.vcf"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.sort.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt sort -o {output.vcf} {input.vcf} 2> {log}
        """

rule vt_uniq:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.sort.vcf",
    output:
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.vcf"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.log",
    container:
        "docker://hydragenetics/vt:2015.11.10"
    shell:
        """
        vt uniq -o {output.vcf} {input.vcf} 2> {log}

rule fix_vcf_af:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.vcf",
    output:
        vcf=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.fix-af.vcf"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.fix-af.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.15.2--py38h7be0bb8_11"
    script:
        "scripts/fix_vcf_af.py"

rule tabix:
    input:
        vcf=f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.fix-af.vcf"
    output:
        tabix=temp(f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.fix-af.vcf.tbi"),
    log:
        f"{decompose_dir}/{{family}}/{{family}}.decomposed.normalized.uniq.fix-af.tabix.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    localrule: True
    shell:
        """
        tabix {input.vcf}
        """
