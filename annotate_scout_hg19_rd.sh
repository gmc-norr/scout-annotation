for sample in <sample>; do

echo "##  INFO  ###    Splitting multiallelic"

# split multiallelic and index
bcftools norm -m -both ${sample}.ensemble.vcf \
| bgzip > ${sample}.mono.vcf.gz
tabix ${sample}.mono.vcf.gz

echo "##  INFO  ###    Annotating with vep"

# annotate with vep
/usr/local/share/ensembl-vep/vep \
--fork 4 \
--distance 5000 \
--buffer_size 20000 \
--assembly GRCh37 \
--fasta /storage/userdata/references/vep_data/hg19/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--cache \
--dir_cache /storage/userdata/references/vep_data/hg19 \
--format vcf \
--vcf \
--dir_plugins /storage/userdata/references/vep_data/hg19/Plugins \
--plugin CADD,/storage/userdata/references/vep_data/hg19/cadd_data/cadd/hg19/whole_genome_SNVs.tsv.gz,/storage/userdata/references/vep_data/hg19/cadd_data/cadd/hg19/InDels.tsv.gz \
--plugin LoFtool,/storage/userdata/references/vep_data/hg19/loftool_data/LoFtool_scores.txt \
--plugin MaxEntScan,/storage/userdata/references/vep_data/hg19/maxEntScan,SWA,NCSS \
--plugin REVEL,/storage/userdata/references/vep_data/hg19/revel/new_tabbed_revel.tsv.gz \
--plugin dbNSFP,/storage/userdata/references/vep_data/hg19/dbNSFP/dbNSFP3.5a_hg19.gz,GERP++_RS,phastCons100way_vertebrate,phyloP100way_vertebrate \
--custom /storage/userdata/references/vep_data/hg19/clinvar/hg19/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT \
--custom /storage/userdata/references/swegen/hg19/swegen_20180409/anon-SweGen_STR_NSPHS_1000samples_SNV_hg19.vcf.gz,SweGen_AF,vcf,exact,0,AF \
--appris \
--biotype \
--canonical \
--ccds \
--domains \
--exclude_predicted \
--force_overwrite \
--hgvs \
--humdiv \
--no_progress \
--no_stats \
--numbers \
--merged \
--polyphen p \
--protein \
--offline \
--regulatory \
--sift p \
--symbol \
--tsl \
--uniprot \
--af_gnomad \
--max_af \
--input_file ${sample}.mono.vcf.gz \
--output_file ${sample}.mono.vep_hg19.vcf \

#echo "##  INFO  ###    Add to locusdb"

# add to loqusdb
#/usr/bin/loqusdb load \
#--variant-file ${sample}.mono.vep_hg19.vcf \
#-f ${sample}.ped \

echo "##  INFO  ###    Annotating with annovar"

# annotate with annovar
/usr/local/bin/annovar/table_annovar.pl ${sample}.mono.vep_hg19.vcf /usr/local/bin/annovar/humandb/ \
-buildver hg19 \
-vcfinput \
-out ${sample}.mono.vep \
-remove \
-protocol spidex -operation f \
-nastring . \
-polish


mv ${sample}.mono.vep.hg19_multianno.vcf ${sample}.mono.vep.anno_hg19.vcf

echo "##  INFO  ###    Removing transcripts not in genepanel"

# remove transcripts not in genepanel
python3.6 /home/liko06/bin/misc-scripts/update_vcf.py --vcf ${sample}.mono.vep.anno_hg19.vcf --panel DCM_genelist.txt --out ${sample}.mono.vep.anno_hg19.genepanel.vcf

echo "##  INFO  ###    Annotating with genmod"

# annotate with models
cat ${sample}.mono.vep.anno_hg19.genepanel.vcf | \
genmod annotate - --annotate_regions | \
genmod models - --family_file ${sample}.ped --vep | \
genmod score - -r -c ../rank_score_model_rd.ini --family_file ${sample}.ped | \
genmod compound - --vep > ${sample}.mono.vep.anno_hg19.genepanel.genmod.vcf

echo "##  INFO  ###    Reformatting vcf to fit Scout"

# reformat GT and caller info to fit with scout
python3.6 /home/liko06/bin/misc-scripts/reformat_vcf_for_scout.py --vcf ${sample}.mono.vep.anno_hg19.genepanel.genmod.vcf --out ${sample}.ann_hg19_scout.vcf


rm ${sample}.mono.vep.hg19_multianno.txt
rm ${sample}.mono.vcf.gz
rm ${sample}.mono.vcf.gz.tbi
rm ${sample}.mono.vep.avinput
rm ${sample}.mono.vep_hg19.vcf
rm ${sample}.mono.vep.anno_hg19.vcf
rm ${sample}.mono.vep.anno_hg19.genepanel.vcf
rm ${sample}.mono.vep.anno_hg19.genepanel.genmod.vcf

echo "##  INFO  ###    Done annotating" 


done
