rule copy_results:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcfanno.genmod.vcf.gz",
        tbi="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcfanno.genmod.vcf.gz.tbi",
        ped=get_ped,
    output:
        vcf=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.scout-annotated.vcf.gz",
        tbi=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.scout-annotated.vcf.gz.tbi",
        ped=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.ped",
    log: f"{config.get('output_directory', 'results')}/{{sample}}/copy_results.log"
    run:
        import logging
        import pathlib
        import subprocess

        logging.basicConfig(
            filename=log[0],
            filemode="a",
            format="%(asctime)s:%(levelname)s: %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
            level=0,
        )
        logging.info("copying %d files", len(input))

        for infile, outfile in zip(input, output):
            logging.info("%s --> %s", infile, pathlib.Path(outfile).resolve())
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
