rule decompose:
    input:
        vcf=lambda wc: get_vcf_file(wc.sample)
    output:
        vcf=temp("results/{sample}/{sample}.decomposed.vcf")
    log: "results/{sample}/{sample}.decomposed.log"
    container: config.get("decompose", {}).get("container", config["default_container"])
    shell:
        """
        vt decompose {input.vcf} -o {output.vcf} 2> {log}
        """