rule decompose:
    input:
        vcf=get_vcf_file
    output:
        vcf=temp("decompose/{sample}/{sample}.decomposed.vcf")
    log: "decompose/{sample}/{sample}.decomposed.log"
    container: config.get("decompose", {}).get("container", config["default_container"])
    shell:
        """
        vt decompose {input.vcf} -o {output.vcf} 2> {log}
        """