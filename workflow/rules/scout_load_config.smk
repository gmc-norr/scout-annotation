rule scout_load_config:
    input:
        vcf="results/{sample}/{sample}.scout-annotated.vcf.gz",
        ped="results/{sample}/{sample}.ped"
    output:
        yaml="results/{sample}/{sample}.load_config.yaml"
    log: "results/{sample}/{sample}.load_config.log"
    container: config.get("scout_load_config", {}).get("container", config.get("default_container", ""))
    params:
        sample_name=lambda wc: wc.sample,
        sex=get_sex,
        phenotype="affected",
        analysis_type=get_analysis_type,
        track=get_track,
        rank_model_version=get_rank_model_version,
        vcf_samples=get_vcf_samples
    script: "../scripts/scout_load_config.py"