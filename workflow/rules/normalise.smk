
rule bcftools_reheader:
    input:
        vcf=get_initial_vcf_file,
    output:
        sample_conversion=temp(
            decompose_dir + "/{family}/{family}.sample_name_conversion.txt"
        ),
        vcf=temp(decompose_dir + "/{family}/{family}.reheadered.vcf"),
    params:
        new_name=lambda wc: family_dict[wc.family]["sample"],
    log:
        decompose_dir + "/{family}/{family}.reheadered.log",
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        """
        echo "$(bcftools query -l {input.vcf} | head -1) {params.new_name}" >> {output.sample_conversion} &&
        bcftools reheader -s {output.sample_conversion} -o {output.vcf} {input.vcf} 2> {log}
        """


rule undecompose:
    input:
        vcf=decompose_dir + "/{family}/{family}.reheadered.vcf",
    output:
        vcf=temp(decompose_dir + "/{family}/{family}.undecomposed.vcf"),
    log:
        decompose_dir + "/{family}/{family}.undecomposed.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.23.3--py312h8f9e533_2"
    script:
        "../scripts/undecompose_vcf.py"


rule bcftools_norm:
    input:
        vcf=decompose_dir + "/{family}/{family}.undecomposed.vcf",
        fasta=config["reference"]["fasta"],
    output:
        vcf=temp(decompose_dir + "/{family}/{family}.normalized.undecomposed.vcf"),
    log:
        decompose_dir + "/{family}/{family}.bcftools-norm.log",
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        """
        bcftools norm -f {input.fasta} --check-ref x -Ou {input.vcf} \\
        | bcftools norm --rm-dup exact -o {output.vcf} 2> {log}
        """


rule rename_info_fields:
    input:
        vcf=decompose_dir + "/{family}/{family}.normalized.undecomposed.vcf",
    output:
        vcf=temp(
            decompose_dir
            + "/{family}/{family}.normalized.undecomposed.renamed_info.vcf"
        ),
    log:
        decompose_dir + "/{family}/{family}.renamed-info.log",
    params:
        rename="INFO/CALLERS FOUND_IN",
    container:
        "docker://hydragenetics/common:0.3.0"
    shell:
        """
        bcftools annotate --rename-annots <(echo "{params.rename}") -O v -o {output.vcf} {input.vcf} 2> {log}
        """


rule rename_callers:
    input:
        vcf=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.vcf",
    output:
        vcf=temp(
            decompose_dir
            + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.vcf"
        ),
    log:
        decompose_dir + "/{family}/{family}.renamed-callers.log",
    params:
        callers_map={"gatk_mutect2": "gatk", "mutect2": "gatk"},
        callers_field="FOUND_IN",
    container:
        "docker://quay.io/biocontainers/pysam:0.23.3--py312h8f9e533_2"
    script:
        "../scripts/rename_callers.py"


rule fix_vcf_af:
    input:
        vcf=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.vcf",
    output:
        vcf=temp(
            decompose_dir
            + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf"
        ),
    log:
        decompose_dir + "/{family}/{family}.fix-af.log",
    container:
        "docker://quay.io/biocontainers/pysam:0.23.3--py312h8f9e533_2"
    script:
        "../scripts/fix_vcf_af.py"


rule bgzip_normalized:
    input:
        vcf=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf",
    output:
        vcfgz=temp(
            decompose_dir
            + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf.gz"
        ),
    log:
        decompose_dir + "/{family}/{family}.bgzip.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcfgz}
        """


rule tabix_normalized:
    input:
        vcf=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf.gz",
    output:
        tabix=temp(
            decompose_dir
            + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf.gz.tbi"
        ),
    log:
        decompose_dir + "/{family}/{family}.fix-af.tabix.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    localrule: True
    shell:
        """
        tabix {input.vcf}
        """
