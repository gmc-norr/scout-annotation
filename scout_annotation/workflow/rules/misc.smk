rule link_bam:
    input:
        bam=get_bam_file,
        bai=get_bai_file,
    output:
        bam=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.bam",
        bai=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.bam.bai",
    log:
        f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.link_bam.log",
    run:
        from pathlib import Path

        outbam = Path(output.bam)
        outbai = Path(output.bai)
        if outbam.exists():
            outbam.unlink()
        if outbai.exists():
            outbai.unlink()
        outbam.symlink_to(Path(input.bam).resolve())
        outbai.symlink_to(Path(input.bai).resolve())


rule copy_peddy_files:
    input:
        het_check="peddy/{family}/{family}.het_check.csv",
        ped_check="peddy/{family}/{family}.ped_check.csv",
        sex_check="peddy/{family}/{family}.sex_check.csv",
        html="peddy/{family}/{family}.html",
        ped="peddy/{family}/{family}.peddy.ped",
        pedigree="madeline2/{family}/{family}.pedigree.svg",
    output:
        het_check=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.peddy.het_check.csv",
        ped_check=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.peddy.ped_check.csv",
        sex_check=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.peddy.sex_check.csv",
        html=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.peddy.html",
        ped=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.peddy.ped",
        pedigree=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.pedigree.svg",
    log:
        f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.copy_peddy_files.log",
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
        logging.info("copying peddy files")

        for infile, outfile in zip(input, output):
            logging.info("%s --> %s", infile, pathlib.Path(outfile).resolve())
            subprocess.run(["cp", infile, outfile], shell=False, check=True)


rule copy_results:
    input:
        vcf=get_annotated_vcf,
        tbi=get_annotated_vcf_index,
        ped=get_family_ped,
    output:
        vcf=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.scout-annotated.vcf.gz",
        tbi=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.scout-annotated.vcf.gz.tbi",
        ped=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.ped",
    log:
        f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.copy_results.log",
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
            subprocess.run(["cp", infile, outfile], shell=False, check=True)


rule tabix:
    input:
        vcf="{filepath}.vcf.gz",
    output:
        tabix=temp("{filepath}.vcf.gz.tbi"),
    log:
        "{filepath}.tabix.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    shell:
        """
        tabix {input.vcf}
        """


rule bgzip:
    input:
        vcf="{filepath}.vcf",
    output:
        vcfgz=temp("{filepath}.vcf.gz"),
    log:
        "{filepath}.bgzip.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcfgz}
        """
