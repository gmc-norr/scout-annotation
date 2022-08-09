rule vep:
    input:
        vcf="decompose/{sample}/{sample}.decomposed.vcf.gz",
        tabix="decompose/{sample}/{sample}.decomposed.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=config["vep"]["cache"],
        plugin=config["vep"]["plugin"],
        plugin_data=config["vep"]["plugin-data"],
        swegen=config["vep"]["swegen"],
        clinvar=config["vep"]["clinvar"]
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.vcf")
    log: "annotation/{sample}/{sample}.decomposed.vep.log"
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
            --plugin CADD,{input.plugin_data}/CADD_1.6/whole_genome_SNVs.tsv.gz,{input.plugin_data}/CADD_1.6/InDels.tsv.gz \\
            --plugin LoFtool,{input.plugin_data}/LoFtool/LoFtool_scores.txt \\
            --plugin MaxEntScan,{input.plugin_data}/MaxEntScan,SWA,NCS \\
            --plugin REVEL,{input.plugin_data}/revel_1.3/new_tabbed_revel.tsv.gz \\
            --plugin dbNSFP,{input.plugin_data}/dbNSFP_4.1a/dbNSFP4.1a_grch37.gz,GERP++_RS,phastCons100way_vertebrate,phyloP100way_vertebrate \\
            --custom {input.clinvar},CLINVAR,vcf,exact,0,CLNSIG,CLNREVSTAT \\
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

rule most_severe_consequence:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.vcf"
    output:
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcf"
    script: "../scripts/most_severe_consequence.py"

rule annovar:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcf"
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.hg19_multianno.vcf"),
        txt=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.hg19_multianno.txt"),
        avinput=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.avinput"),
    log: "annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.log"
    params:
        prefix="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar"
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
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.hg19_multianno.vcf",
        ped=get_ped,
        rank_model=get_rank_model
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.annovar.genmod.vcf")
    log: "annotation/{sample}/{sample}.genmod.log"
    container: config.get("genmod", {}).get("container", config["default_container"])
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

rule genmod_rankmodel:
    output:
        rank_model="rank_model/{type}_rank_model_v{version}.ini",
    params:
        uri=lambda wc: config["genmod"][f"{wc.type}_rank_model_uri"].format(version=wc.version),
        filename=lambda wc: f"{wc.type}_rank_model_v{wc.version}.ini"
    log: "rank_model/{type}_rank_model_v{version}.log"
    shell:
        """
        echo "fetching {params.uri}" > {log}
        curl -fsSL {params.uri} > {output.rank_model} 2>> {log}
        """