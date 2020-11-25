#!/bin/bash
echo Variant Calling starting at:
date
path=$1
name=$2
infotable=$3
bed=$4


gatk MarkDuplicates -I $path$name.sort.bam -M $path$name.dupMetrix.txt -O $path$name.rmdup.bam --REMOVE_DUPLICATES true

#addReadGroup
gatk AddOrReplaceReadGroups -I $path$name.rmdup.bam -O $path$name.rmdup.rg.bam -LB Solexa-272222 -PL illumina -PU H0164ALXX140820.2 -SM $name

#CreateSequenceDict
#gatk CreateSequenceDictionary -R hg19_v0_Homo_sapiens_assembly19.fasta -O hg19_v0_Homo_sapiens_assembly19.dict

#BSQR
gatk BaseRecalibrator \
  -I $path$name.rmdup.rg.bam \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  --known-sites /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_dbsnp_138.b37.vcf.gz  \
  --known-sites /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  -O $path$name.bqsr.grp


gatk ApplyBQSR \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -I $path$name.rmdup.rg.bam \
  --bqsr-recal-file $path$name.bqsr.grp \
  -O $path$name.bqsr.bam


#HaplotypeCaller
gatk HaplotypeCaller \
  -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -L /home/guoping/bed_file/$bed \
  -I $path$name.bqsr.bam \
  -O $path$name.vcf.gz
  #  -ERC GVCF \


#gunzip $1.vcf.gz


#conda active gatkk
#
#
#***********************Genotype refinement workflow



#**********************CNNscore variant scoring*******************
gatk CNNScoreVariants \
   -I $path$name.bqsr.bam \
   -V $path$name.vcf.gz \
   -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
   -O $path$name.CNN.vcf \
   --tensor-type read_tensor

#**********************FilterVariantTranches variant filtering******************
gatk FilterVariantTranches \
   -V $path$name.CNN.vcf \
   --resource /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_omni2.5.b37.vcf.gz \
   --resource /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz \
   --resource /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
   --info-key CNN_2D \
   --snp-tranche 99.95 \
   --indel-tranche 99.4 \
   --invalidate-previous-filters \
   -O $path$name.CNNfiltered.vcf


#CalculateGenotypePosteriors
gatk CalculateGenotypePosteriors \
    -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
    -V $path$name.CNNfiltered.vcf \
    -supporting /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_1000G_omni2.5.b37.vcf.gz \
    -supporting /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_hapmap_3.3.b37.vcf.gz \
    -O $path$name.genotypeCorrected.vcf

#VariantFilter
#if you use genotype filter you will get FT values
#if you use filter, there is no FT values

gatk VariantFiltration \
   -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
   -V $path$name.genotypeCorrected.vcf \
   -O $path$name.prefiltered.vcf \
   -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name 'FAIL' \
   --genotype-filter-expression 'GQ <20' \
   --genotype-filter-name 'FAIL'

#Get ONLY PASS FILTER
awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' $path$name.prefiltered.vcf > $path$name.final.vcf

convert2annovar.pl -format vcf4 $path$name.final.vcf -outfile $path$name.avinput --withfreq -includeinfo

table_annovar.pl $path$name.avinput /data2/guoping/software/annovar/humandb -buildver hg19 -out $path$name \
                         -remove --protocol refGene,exac03,clinvar_20200602,cytoBand,esp6500siv2_all,genomicSuperDups,gwasCatalog,snp142 \
                         -operation g,f,f,r,f,r,r,f \
                         -nastring . \
                         --polish \
                         -xreffile ~/software/annovar/humandb/hg19_refGene.txt \
                         --otherinfo \
                         --thread 32

variantCall_reformat.py $path$infotable

echo exit status: $?

echo Variant Calling ending at:
date
