from pathlib import Path

rule scout_load_config:
    input:
        vcf=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.scout-annotated.vcf.gz",
        ped=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.ped",
    output:
        yaml=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.load_config.yaml",
    log: f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.load_config.log"
    container: "docker://python:3.10.7-slim"
    params:
        sample_name=lambda wc: wc.sample,
        sex=get_sex,
        phenotype="affected",
        analysis_type=get_analysis_type,
        track=get_track,
        rank_model_version=get_rank_model_version,
        vcf_samples=get_vcf_samples,
        panels=get_sample_panels,
        include_bam=lambda wc: isinstance(get_bam_file(wc), Path),
    script: "../scripts/scout_load_config.py"
