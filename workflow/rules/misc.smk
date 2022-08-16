rule copy_results:
    input: infiles
    output: outfiles
    log: f"{config.get('output_directory', 'results')}/copy_results.log"
    run:
        import subprocess

        for infile, outfile in zip(input, output):
            subprocess.run([
                "cp", infile, outfile
            ], shell=False, check=True)

rule tabix:
    input:
        vcf="{filepath}.vcf.gz"
    output:
        tabix=temp("{filepath}.vcf.gz.tbi")
    log: "{filepath}.tabix.log"
    container: config.get("tabix", {}).get("container", config["default_container"])
    shell:
        """
        tabix {input.vcf}
        """

rule bgzip:
    input:
        vcf="{filepath}.vcf"
    output:
        vcfgz=temp("{filepath}.vcf.gz")
    log: "{filepath}.bgzip.log"
    container: config.get("bgzip", {}).get("container", config["default_container"])
    shell:
        """
        bgzip {input.vcf}
        """
