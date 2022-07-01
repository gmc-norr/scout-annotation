rule mock_ped:
    input:
        vcf=get_vcf_file
    output:
        ped=temp("mock_ped/{sample}.ped")
    log: "mock_ped/{sample}.mock_ped.log"
    script: "../scripts/mock_ped.py"