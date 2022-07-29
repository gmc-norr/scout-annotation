rule scout_load_config:
    input:
        vcf="results/{sample}/{sample}.scout-annotated.vcf.gz",
        ped="results/{sample}/{sample}.ped"
    output:
        yaml="results/{sample}/{sample}.load_config.yaml"
    params:
        sample_name=lambda wc: wc.sample,
        sex=get_sex,
        phenotype="affected",
        analysis_type=get_analysis_type,
        track=get_track,
        vcf_samples=get_vcf_samples
    script: "../scripts/scout_load_config.py"