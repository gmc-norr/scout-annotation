rule mock_ped:
    input:
        vcf=get_reheadered_vcf_file
    output:
        ped=temp("mock_ped/{sample}.ped")
    log: "mock_ped/{sample}.mock_ped.log"
    localrule: True
    script: "../scripts/mock_ped.py"
