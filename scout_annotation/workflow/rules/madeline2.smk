rule madeline2_data:
    input:
        ped=get_family_ped,
    output:
        madeline2_data=temp("madeline2/{family}/{family}.madeline2_data.tsv"),
    log:
        "madeline2/{family}/{family}.madeline2_data.log",
    container:
        config.get("madeline2_data", {}).get("container", config["default_container"])
    script:
        "../scripts/madeline2_data.py"


rule madeline2_svg:
    input:
        madeline2_data="madeline2/{family}/{family}.madeline2_data.tsv",
    output:
        svg=temp("madeline2/{family}/{family}.pedigree.svg"),
    log:
        "madeline2/{family}/{family}.madeline2_svg.log",
    container:
        config.get("madeline2_svg", {}).get("container", "docker://maehler/madeline2")
    shell:
        """
        (madeline2 \
            --embedded \
            -o madeline2/{wildcards.family}/{wildcards.family}.pedigree \
            {input.madeline2_data}) \
            > {log} 2>&1
        """
