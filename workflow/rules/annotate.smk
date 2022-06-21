rule vep:
    input:
        vcf="results/{sample}/{sample}.decomposed.vcf.gz",
        tabix="results/{sample}/{sample}.decomposed.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=config["vep"]["cache"],
        plugin=config["vep"]["plugin"],
        plugin_data=config["vep"]["plugin-data"],
        swegen=config["vep"]["swegen"]
    output:
        vcf="results/{sample}/{sample}.decomposed.vep.vcf"
    log: "results/{sample}/{sample}.decomposed.vep.log"
    params:
        mode=config.get("vep", {}).get("mode", ""),
        cache_type=config.get("vep", {}).get("cache_type", "merged")
    threads: resources.get("vep", {}).get("threads", resources["default_resources"]["threads"])
    container: config.get("vep", {}).get("container", config["default_container"])
    shell:
        """
        vep {params.mode} \\
            -i {input.vcf} \\
            -o {output.vcf} \\
            --vcf \\
            --assembly GRCh37 \\
            --dir_plugins {input.plugin} \\
            --plugin CADD,{input.plugin_data}/cadd_data/cadd/hg19/whole_genome_SNVs.tsv.gz,{input.plugin_data}/cadd_data/cadd/hg19/InDels.tsv.gz \\
            --plugin LoFtool,{input.plugin_data}/loftool_data/LoFtool_scores.txt \\
            --plugin MaxEntScan,{input.plugin_data}/maxEntScan,SWA,NCS \\
            --plugin REVEL,{input.plugin_data}/revel/new_tabbed_revel.tsv.gz \\
            --plugin dbNSFP,{input.plugin_data}/dbNSFP/dbNSFP3.5a_hg19.gz,GERP++_RS,phastCons100way_vertebrate,phyloP100way_vertebrate \\
            --custom {input.plugin_data}/clinvar/hg19/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT \\
            --custom {input.swegen},SweGen,vcf,exact,0,AF \\
            --dir_cache {input.cache} \\
            --fork {threads} \\
            --distance 5000 \\
            --buffer_size 2000 \\
            --fasta {input.fasta} \\
            --check_existing \\
            --pick \\
            --sift b \\
            --polyphen b \\
            --ccds \\
            --uniprot \\
            --hgvs \\
            --symbol \\
            --numbers \\
            --domains \\
            --regulatory \\
            --canonical \\
            --protein \\
            --biotype \\
            --uniprot \\
            --tsl \\
            --appris \\
            --gene_phenotype \\
            --af \\
            --af_1kg \\
            --af_gnomad \\
            --max_af \\
            --pubmed \\
            --variant_class \\
            --exclude_predicted \\
            --humdiv \\
            --no_stats \\
            --{params.cache_type} \\
            &> {log}
        """

rule annovar:
    input:
        vcf="results/{sample}/{sample}.decomposed.vep.vcf"
    output:
        vcf="results/{sample}/{sample}.decomposed.vep.annovar.hg19_multianno.vcf",
        txt="results/{sample}/{sample}.decomposed.vep.annovar.hg19_multianno.txt",
        avinput="results/{sample}/{sample}.decomposed.vep.annovar.avinput",
    log: "results/{sample}/{sample}.decomposed.vep.annovar.log"
    params:
        prefix="results/{sample}/{sample}.decomposed.vep.annovar"
    shell:
        """
        /usr/local/bin/annovar/table_annovar.pl \\
            --buildver hg19 \\
            --vcfinput \\
            --out {params.prefix} \\
            --remove \\
            --protocol spidex \\
            --operation f \\
            --nastring . \\
            --polish \\
            {input.vcf} \\
            /usr/local/bin/annovar/humandb &> {log}
        """

rule genmod:
    input:
        vcf="results/{sample}/{sample}.decomposed.vep.annovar.hg19_multianno.vcf",
        ped=get_ped,
        rank_model=config["genmod"]["rank_model"]
    output:
        vcf="results/{sample}/{sample}.decomposed.vep.annovar.genmod.vcf"
    conda: "../env/genmod.yaml"
    shell:
        """
        genmod annotate \\
            --annotate_regions \\
            {input.vcf} | \\
        genmod models \\
            --family_file {input.ped} \\
            --vep \\
            - | \\
        genmod score \\
            --rank_results \\
            --score_config {input.rank_model} \\
            - | \\
        genmod compound \\
            --vep \\
            - > {output.vcf}
        """