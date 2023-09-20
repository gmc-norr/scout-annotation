rule peddy:
    input:
        vcf="merge_trio/{family}_trio.vcf.gz",
        tbi="merge_trio/{family}_trio.vcf.gz.tbi",
        ped=get_family_ped,
    output:
        het_check=temp("peddy/{family}/{family}.het_check.csv"),
        html=temp("peddy/{family}/{family}.html"),
        ped=temp("peddy/{family}/{family}.peddy.ped"),
        ped_check=temp("peddy/{family}/{family}.ped_check.csv"),
        sex_check=temp("peddy/{family}/{family}.sex_check.csv"),
    log:
        "peddy/{family}/{family}.peddy.log",
    params:
        each=config.get("peddy", {}).get("each"),
        plot=config.get("peddy", {}).get("plot"),
        sites=config.get("peddy", {}).get("sites"),
    threads: resources.get("peddy", {}).get("threads", resources["default_resources"]["threads"])
    container:
        config.get("peddy", {}).get("container", "docker://hydragenetics/peddy:0.4.8")
    script:
        "../scripts/peddy.py"
