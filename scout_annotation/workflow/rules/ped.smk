rule mock_ped:
    input:
        vcf=get_reheadered_vcf_file
    output:
        ped=temp("mock_ped/{sample}.ped")
    params:
        family_id = get_family_id
    log: "mock_ped/{sample}.mock_ped.log"
    localrule: True
    container: config.get("mock_ped", {}).get("container", config.get("default_container", ""))
    script: "../scripts/mock_ped.py"
