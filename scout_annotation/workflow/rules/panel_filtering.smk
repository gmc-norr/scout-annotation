checkpoint panel_filtering:
    input:
        vcf="annotation/{sample}/{sample}.annotated.vcf",
        panels=get_panel_files,
    output:
        vcf=temp("annotation/{sample}/{sample}.annotated.panel_filtered.vcf"),
    log: "annotation/{sample}/{sample}.panel_filtering.log",
    params:
        hard_filter=config.get("panel_filtering", {}).get("hard_filter", True),
    conda: "../env/panel_filtering.yaml",
    script: "../scripts/panel_filtering.py"
