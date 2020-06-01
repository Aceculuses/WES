#trimmomatic PE -phred33 -threads 6 -quiet -validatePairs Ca1232_1.fq Ca1232_2.fq Ca1232_1.paired.fq Ca1232_1.unpaired.fq Ca1232_2.paired.fq Ca1232_2.paired.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36

fastp -i ${1}_1.fq -o ${1}_1.clean.fastq -I ${1}_2.fq -O ${1}_2.clean.fastq \
      -A -x --cut_mean_quality 20 --cut_front_window_size 4 --cut_front_mean_quality 20 \
      --cut_tail_window_size 4 --cut_tail_mean_quality 20 -n 5 --cut_front 20 \
      -f 15 -t 5 -F 15 -T 5

#fastuniq -i list.txt -t q -o MF04_1.tmp.uniq.fastq -p MF04_2.tmp.uniq.fastq

#mapping
bwa mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta ${1}_1.clean.fastq ${1}_2.clean.fastq > $1.sam

#samtools
samtools view -bS $1.sam > $1.bam
samtools sort $1.bam -o $1.sort.bam
samtools index $1.sort.bam

#mark duplicate
gatk MarkDuplicates -I $1.sort.bam -M $1.dupMetrix.txt -O $1.rmdup.bam --REMOVE_DUPLICATES true

#addReadGroup
gatk AddOrReplaceReadGroups -I $1.rmdup.bam -O $1.rmdup.rg.bam -LB Solexa-272222 -PL illumina -PU H0164ALXX140820.2 -SM $1

#CreateSequenceDict
#gatk CreateSequenceDictionary -R hg19_v0_Homo_sapiens_assembly19.fasta -O hg19_v0_Homo_sapiens_assembly19.dict

#BSQR
gatk BaseRecalibrator \
  -I $1.rmdup.rg.bam \
  -R /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  --known-sites /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_dbsnp_138.b37.vcf.gz  \
  --known-sites /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  -O $1.bqsr.grp 


gatk ApplyBQSR \
  -R /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -I $1.rmdup.rg.bam \
  --bqsr-recal-file $1.bqsr.grp \
  -O $1.bqsr.bam


#HaplotypeCaller
gatk HaplotypeCaller \
  -R /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
  -L exon.bed \
  -I $1.bqsr.bam \
  -O $1.vcf.gz
  #  -ERC GVCF \

#
gunzip $1.vcf.gz


#VariantFilter
#if you use genotype filter you will get FT values
#if you use filter, there is no FT values
gatk VariantFiltration \
   -R /home/ubuntu/guoping/07-test/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
   -V $1.vcf \
   -O $1.final.noFT.vcf \
   -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name 'FAIL'

#Get ONLY PASS FILTER
awk -F '\t' '{if($0 ~ /#/) print; else if($7 == "PASS") print}' $1.final.noFT.vcf > $1.final.PASSnoFT.vcf

#liftover
CrossMap.py vcf /home/ubuntu/software/annovar/humandb/GRCh37_to_GRCh38.chain.gz $1.final.PASSnoFT.vcf \
                /home/ubuntu/software/annovar/humandb/Homo_sapiens.GRCh38.dna.toplevel.fa $1.final.PASSnoFThg38.vcf


#annovar
convert2annovar.pl -format vcf4 $1.final.PASSnoFThg38.vcf -outfile ${1}PASSnoFTihg38.avinput --withfreq -includeinfo

#table_annovar.pl Ca1232PASSnoFT.avinput ~/software/annovar/humandb/ -buildver hg19 -out Ca1232PASSnoFT \
#                         -remove --protocol refGene,exac03,clinvar_20191231,cytoBand,esp6500siv2_all,genomicSuperDups,gwasCatalog,snp142 \
#                         -operation g,f,f,r,f,r,r,f \
#                         -nastring . \
#                         --polish \
#                         -xreffile ~/software/annovar/humandb/hg19_refGene.txt \
#                         --otherinfo \
#                         --thread 32


table_annovar.pl ${1}PASSnoFTihg38.avinput ~/software/annovar/humandb/ -buildver hg38 -out ${1}PASSnoFThg38 \
                         -remove --protocol refGene,ensGene,exac03,clinvar_20200113,esp6500siv2_all,avsnp150,dbnsfp35c,1000g2015aug,gnomad211_exome \
                         -operation g,g,f,f,f,f,f,f,f \
                         -nastring . \
                         --polish \
                         -xreffile ~/software/annovar/humandb/hg38_refGene.txt 
                         --otherinfo \
                         --thread 32


#reformat
#reformat.py /home/ubuntu/guoping/03-annotation/3.5-result/hg19_self_results/reformat/ /home/ubuntu/guoping/03-annotation/3.5-result/hg19_self_results/results/ hg19 20191231 Ca1232.hg19_multianno.txt 
