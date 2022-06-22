rule mock_ped:
    input:
        vcf=get_vcf_file
    output:
        ped="results/{sample}/{sample}.ped"
    log: "results/{sample}/{sample}.mock_ped.log"
    script: "../scripts/mock_ped.py"
