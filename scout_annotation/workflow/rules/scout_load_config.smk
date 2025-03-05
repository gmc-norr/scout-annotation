from pathlib import Path


rule scout_load_config_family:
    input:
        vcf=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.scout-annotated.vcf.gz",
        ped=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.ped",
        sample_configs=lambda wc: expand(
            "load_config/{{family}}_{sample}.load_config.yaml",
            sample=get_family_samples(wc),
        ),
        peddy_ped=lambda wc: get_peddy_file(wc, "ped"),
        peddy_ped_check=lambda wc: get_peddy_file(wc, "ped_check"),
        peddy_sex_check=lambda wc: get_peddy_file(wc, "sex_check"),
        madeline2_svg=lambda wc: get_madeline2_svg(wc),
    output:
        yaml=f"{config.get('output_directory', 'results')}/{{family}}/{{family}}.load_config.yaml",
    params:
        type="family",
        track=get_track,
        msi_score=get_msi(),
        hrd_score=get_hrd(),
        tmb_score=get_tmb(),
        owner=get_case_owner,
        panels=get_family_panels,
        rank_model_version=get_rank_model_version,
        rank_score_threshold=-1000,
    container:
        "docker://python:3.10.7-slim"
    script:
        "../scripts/scout_load_config.py"


rule scout_load_config_sample:
    input:
        ped=get_family_ped,
    output:
        yaml=temp("load_config/{family}_{sample}.load_config.yaml"),
    params:
        type="sample",
        include_bam=lambda wc: isinstance(get_bam_file(wc), Path),
        sex=get_sample_sex,
        analysis_type=get_analysis_type,
    container:
        "docker://python:3.10.7-slim"
    script:
        "../scripts/scout_load_config.py"
