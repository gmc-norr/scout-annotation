checkpoint panel_filtering:
    input:
        vcf="annotation/{family}/{family}.annotated.vcf",
        panels=get_panel_files,
    output:
        vcf=temp("annotation/{family}/{family}.annotated.panel_filtered.vcf"),
    log:
        "annotation/{family}/{family}.panel_filtering.log",
    params:
        hard_filter=config.get("panel_filtering", {}).get("hard_filter", True),
    conda:
        "../env/panel_filtering.yaml"
    script:
        "../scripts/panel_filtering.py"
