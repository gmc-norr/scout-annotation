rule tabix:
    input:
        vcf="{filepath}.vcf.gz"
    output:
        tabix="{filepath}.vcf.gz.tbi"
    log: "{filepath}.tabix.log"
    container: config.get("tabix", {}).get("container", config["default_container"])
    shell:
        """
        tabix {input.vcf}
        """

rule copy_results:
    input: "results/{sample}/{sample}.decomposed.vep.annovar.genmod.vcf.gz"
    output: "results/{sample}/{sample}.scout-annotated.vcf.gz"
    shell:
        """
        cp {input} {output}
        """

rule bgzip:
    input:
        vcf="{filepath}.vcf"
    output:
        vcfgz="{filepath}.vcf.gz"
    log: "{filepath}.bgzip.log"
    container: config.get("bgzip", {}).get("container", config["default_container"])
    shell:
        """
        bgzip {input.vcf}
        """
