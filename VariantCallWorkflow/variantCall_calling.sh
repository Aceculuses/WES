echo Variant Calling part2 starting at:
date

gatk MarkDuplicates -I $1.sort.bam -M $1.dupMetrix.txt -O $1.rmdup.bam --REMOVE_DUPLICATES true

#addReadGroup
gatk AddOrReplaceReadGroups -I $1.rmdup.bam -O $1.rmdup.rg.bam -LB Solexa-272222 -PL illumina -PU H0164ALXX140820.2 -SM $1

#CreateSequenceDict
#gatk CreateSequenceDictionary -R hg19_v0_Homo_sapiens_assembly19.fasta -O hg19_v0_Homo_sapiens_assembly19.dict

#BSQR
gatk BaseRecalibrator \
  -I $1.rmdup.rg.bam \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  --known-sites /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_dbsnp_138.b37.vcf.gz  \
  --known-sites /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  -O $1.bqsr.grp


gatk ApplyBQSR \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -I $1.rmdup.rg.bam \
  --bqsr-recal-file $1.bqsr.grp \
  -O $1.bqsr.bam


#HaplotypeCaller
gatk HaplotypeCaller \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -L exon_hg19_v0.interval_list \
  -I $1.bqsr.bam \
  -O $1.vcf.gz
  #  -ERC GVCF \


#gunzip $1.vcf.gz


#conda active gatkk
#
#
#***********************Genotype refinement workflow



#**********************CNNscore variant scoring*******************
gatk CNNScoreVariants \
   -I $1.bqsr.bam \
   -V $1.vcf.gz \
   -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
   -O $1.CNN.vcf \
   --tensor-type read_tensor

#**********************FilterVariantTranches variant filtering******************
gatk FilterVariantTranches \
   -V $1.CNN.vcf \
   --resource /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_omni2.5.b37.vcf.gz \
   --resource /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz \
   --info-key CNN_2D \
   --snp-tranche 99.95 \
   --indel-tranche 99.4 \
   --invalidate-previous-filters \
   -O $1.CNNfiltered.vcf


#CalculateGenotypePosteriors
gatk CalculateGenotypePosteriors \
    -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
    -V $1.CNNfiltered.vcf \
    -supporting /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_omni2.5.b37.vcf.gz \
    -supporting /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_hapmap_3.3.b37.vcf.gz \
    -O $1.genotypeCorrected.vcf

#VariantFilter
#if you use genotype filter you will get FT values
#if you use filter, there is no FT values

gatk VariantFiltration \
   -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
   -V $1.genotypeCorrected.vcf \
   -O $1.prefiltered.vcf \
   -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name 'FAIL' \
   --genotype-filter-expression 'GQ <20' \
   --genotype-filter-name 'FAIL'

#Get ONLY PASS FILTER
awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' $1.prefiltered.vcf > $1.preFTfiltered.vcf

echo exit status: $?

echo Variant Calling part2 ending at:
date

