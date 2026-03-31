rule vep:
    input:
        vcf=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf.gz",
        tabix=decompose_dir
        + "/{family}/{family}.normalized.undecomposed.renamed_info.renamed_callers.fix-af.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=config["vep"]["cache"],
        plugin=config["vep"]["plugin"],
        plugin_data=config["vep"]["plugin-data"],
        swegen=config["vep"]["swegen"],
        clinvar=config["vep"]["clinvar"],
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.vep.vcf"),
    log:
        annotation_dir + "/{family}/{family}.vep.log",
    params:
        mode=config.get("vep", {}).get("mode", ""),
        cache_type=config.get("vep", {}).get("cache_type", "merged"),
    container:
        "docker://hydragenetics/vep:105"
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


rule vcfanno:
    input:
        vcf=annotation_dir + "/{family}/{family}.vep.vcf",
        toml=get_vcfanno_config(),
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.vep.vcfanno.vcf"),
    log:
        annotation_dir + "/{family}/{family}.vcfanno.log",
    params:
        base_path=config.get("vcfanno", {}).get("base_path", ""),
    container:
        "docker://clinicalgenomics/vcfanno:0.3.2"
    shell:
        """
        vcfanno \\
            -base-path {params.base_path} \\
            -p {threads} \\
            {input.toml} \\
            {input.vcf} \\
            > {output.vcf} \\
            2> {log}
        """


rule most_severe_consequence:
    input:
        vcf=annotation_dir + "/{family}/{family}.vep.vcfanno.vcf",
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.annotated.vcf"),
    log:
        annotation_dir + "/{family}/{family}.most_severe_consequence.log",
    container:
        "docker://quay.io/biocontainers/cyvcf2:0.32.1--py311h921ead3_0"
    script:
        "../scripts/most_severe_consequence.py"


rule genmod_annotate:
    input:
        vcf=annotation_dir + "/{family}/{family}.annotated.vcf",
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.genmod_annotate.vcf"),
    log:
        annotation_dir + "/{family}/{family}.genmod_annotate.log",
    container:
        "docker://quay.io/biocontainers/genmod:3.7.4--pyh5e36f6f_0"
    shell:
        """
        genmod annotate \\
            --annotate_regions \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule mock_ped:
    input:
        vcf=decompose_dir + "/{family}/{family}.reheadered.vcf",
    output:
        ped=temp(annotation_dir + "/{family}/{family}.ped"),
    log:
        annotation_dir + "/{family}/{family}.mock_ped.log",
    localrule: True
    container:
        config.get("mock_ped", {}).get("container", config.get("default_container", ""))
    script:
        "../scripts/mock_ped.py"


rule genmod_models:
    input:
        vcf=annotation_dir + "/{family}/{family}.genmod_annotate.vcf",
        ped=annotation_dir + "/{family}/{family}.ped",
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.genmod_models.vcf"),
    log:
        annotation_dir + "/{family}/{family}.genmod_models.log",
    container:
        "docker://quay.io/biocontainers/genmod:3.7.4--pyh5e36f6f_0"
    shell:
        """
        genmod models \\
            --family_file {input.ped} \\
            --vep \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_score:
    input:
        vcf=annotation_dir + "/{family}/{family}.genmod_models.vcf",
        ped=annotation_dir + "/{family}/{family}.ped",
        rank_model=get_genmod_rank_model(),
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.genmod_score.vcf"),
    log:
        annotation_dir + "/{family}/{family}.genmod_score.log",
    container:
        "docker://quay.io/biocontainers/genmod:3.7.4--pyh5e36f6f_0"
    shell:
        """
        genmod score \\
            --rank_results \\
            --family_file {input.ped} \\
            --score_config {input.rank_model} \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule genmod_compound:
    input:
        vcf=annotation_dir + "/{family}/{family}.genmod_score.vcf",
    output:
        vcf=temp(annotation_dir + "/{family}/{family}.annotated.genmod.vcf"),
    log:
        annotation_dir + "/{family}/{family}.genmod_compound.log",
    container:
        "docker://quay.io/biocontainers/genmod:3.7.4--pyh5e36f6f_0"
    shell:
        """
        genmod compound \\
            --vep \\
            {input.vcf} > {output.vcf} 2> {log}
        """


rule bgzip_annotated:
    input:
        vcf=annotation_dir + "/{family}/{family}.annotated.genmod.vcf",
    output:
        vcfgz=temp(annotation_dir + "/{family}/{family}.annotated.genmod.vcf.gz"),
    log:
        annotation_dir + "/{family}/{family}.annotated.genmod.bgzip.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcfgz}
        """


rule tabix_annotated:
    input:
        vcf=annotation_dir + "/{family}/{family}.annotated.genmod.vcf.gz",
    output:
        tabix=temp(annotation_dir + "/{family}/{family}.annotated.genmod.vcf.gz.tbi"),
    log:
        annotation_dir + "/{family}/{family}.annotated.genmod.tabix.log",
    container:
        "docker://hydragenetics/common:0.1.1"
    localrule: True
    shell:
        """
        tabix {input.vcf}
        """
