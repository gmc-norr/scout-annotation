rule generate_d4:
    input: 
        alignment=get_bam_file,
        fasta_fai=f"{config['reference']['fasta']}.fai",
    output:
        d4_file=f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.d4",
    log: f"{config.get('output_directory', 'results')}/{{sample}}/{{sample}}.d4tools.log"
    threads: resources.get("d4tools", {}).get("threads", resources["default_resources"]["threads"])
    container: 
        "docker://clinicalgenomics/d4tools:2.0"
    shell:
        """
        d4tools create -F '~1024' -Azr {input.fasta_fai} -t {threads} {input.alignment} {output.d4_file} > {log} 2>&1
        """