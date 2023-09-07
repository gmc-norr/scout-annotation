checkpoint vcf_filtering:
    input:
        vcf="{file}.vcf",
        filter_definition=get_vcf_filter,
    output:
        vcf=temp("{file}.filter-{tag}.vcf"),
    log:
        "{file}.filter-{tag}.log",
    params:
        filter_expression=get_vembrane_expression,
    container:
        "docker://quay.io/biocontainers/vembrane:0.14.0--pyhdfd78af_0"
    shell:
        """
        vembrane \\
            filter \\
            --annotation-key CSQ \\
            --output-fmt vcf \\
            --output {output.vcf} \\
            '{params.filter_expression}' \\
            {input.vcf} 2> {log}
        """
