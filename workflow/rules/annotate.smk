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
            --plugin LoFtool,{input.plugin}/LoFtool_scores.txt \\
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
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcf")
    log: "annotation/{sample}/{sample}.most_severe_consequence.log"
    conda: "../env/most_severe_consequence.yaml"
    script: "../scripts/most_severe_consequence.py"

rule vcfanno:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcf",
        toml=get_vcfanno_config
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcfanno.vcf")
    log: "annotation/{sample}/{sample}.vcfanno.log"
    params:
        base_path=config.get("vcfanno", {}).get("base_path", "")
    container: config.get("vcfanno", {}).get("container", config["default_container"])
    shell:
        """
        vcfanno \\
            -base-path {params.base_path} \\
            {input.toml} \\
            {input.vcf} \\
            > {output.vcf} \\
            2> {log}
        """

rule vcfanno_config:
    output:
        toml="rank_model/grch{build}_{track}_vcfanno_config_{version}.toml"
    log: "rank_model/grch{build}_{track}_vcfanno_config_{version}.log"
    params:
        uri=lambda wc: config["vcfanno"]["config_uri"].format(track=wc.track, version=wc.version),
        extra=config.get("vcfanno_config", {}).get("extra", "")
    container: config.get("vcfanno_config", {}).get("container", config["default_container"])
    shell:
        """
        echo "fetching {params.uri}" > {log}
        curl {params.extra} -fsSL {params.uri} > {output.toml} 2>> {log}
        """


rule genmod_annotate:
    input:
        vcf="annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcfanno.vcf",
    output:
        vcf=temp("annotation/{sample}/{sample}.genmod_annotate.vcf"),
    log: "annotation/{sample}/{sample}.genmod_annotate.log"
    container: config.get("genmod", {}).get("container", config["default_container"])
    shell:
        """
        genmod annotate \\
            --annotate_regions \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_models:
    input:
        vcf="annotation/{sample}/{sample}.genmod_annotate.vcf",
        ped=get_ped,
    output:
        vcf=temp("annotation/{sample}/{sample}.genmod_models.vcf"),
    log: "annotation/{sample}/{sample}.genmod_models.log"
    container: config.get("genmod", {}).get("container", config["default_container"])
    shell:
        """
        genmod models \\
            --family_file {input.ped} \\
            --vep \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_score:
    input:
        vcf="annotation/{sample}/{sample}.genmod_models.vcf",
        rank_model=get_rank_model,
    output:
        vcf=temp("annotation/{sample}/{sample}.genmod_score.vcf"),
    log: "annotation/{sample}/{sample}.genmod_score.log"
    container: config.get("genmod", {}).get("container", config["default_container"])
    shell:
        """
        genmod score \\
            --rank_results \\
            --score_config {input.rank_model} \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_compound:
    input:
        vcf="annotation/{sample}/{sample}.genmod_score.vcf",
    output:
        vcf=temp("annotation/{sample}/{sample}.decomposed.vep.most_severe_csq.vcfanno.genmod.vcf"),
    log: "annotation/{sample}/{sample}.genmod_compound.log"
    container: config.get("genmod", {}).get("container", config["default_container"])
    shell:
        """
        genmod compound \\
            --vep \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_rankmodel:
    output:
        rank_model="rank_model/{track}_rank_model_{version}.ini"
    params:
        uri=lambda wc: config["genmod"]["rank_model_uri"].format(track=wc.track, version=wc.version),
        extra=config.get("genmod_rankmodel", {}).get("extra", "")
    log: "rank_model/{track}_rank_model_{version}.log"
    container: config.get("genmod_rankmodel", {}).get("container", config["default_container"])
    shell:
        """
        echo "fetching {params.uri}" > {log}
        curl {params.extra} -fsSL {params.uri} > {output.rank_model} 2>> {log}
        """
