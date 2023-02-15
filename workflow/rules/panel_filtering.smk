rule panel_filtering:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.vcf",
        panels=get_panels,
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.panel_filtered.vcf"),
    log: "annotation/{sample}/{sample}.panel_filtering.log",
    params:
        hard_filter=config.get("panel_filtering", {}).get("hard_filter", True),
    script: "../scripts/panel_filtering.py"
