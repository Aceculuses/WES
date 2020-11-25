#!/bin/bash
echo Variant Calling starting at:
date
inputPath=$1
prefix=$2
outputPath=$3

fastp -i $inputPath${prefix}*R1*.fastq.gz -o $outputPath${prefix}_1.clean.fastq -I $inputPath${prefix}*R2*.fastq.gz -O $outputPath${prefix}_2.clean.fastq \
      -A -x --cut_mean_quality 20 --cut_front_window_size 4 --cut_front_mean_quality 20 \
      --cut_tail_window_size 4 --cut_tail_mean_quality 20 -n 5 --cut_front 20 \
      -f 15 -t 5 -F 15 -T 5 \
      --detect_adapter_for_pe

#fastuniq -i list.txt -t q -o MF04_1.tmp.uniq.fastq -p MF04_2.tmp.uniq.fastq

#bwa mapping
#bwa mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta ${1}_1.clean.fastq ${1}_2.clean.fastq > $1.sam

#bwa-mem2
#bwa-mem2 mem /data2/guoping/Databases/Broad_bundle_hg19/mem2/hg19_v0_Homo_sapiens_assembly19 $outputPath${prefix}_1.clean.fastq $outputPath${prefix}_2.clean.fastq -o $outputPath$prefix.sam -t 40

#bowtie2
bowtie2 -x /data2/guoping/Databases/Broad_bundle_hg19/bowtie2/hg19_v0_Homo_sapiens_assembly19 -1 $outputPath${prefix}_1.clean.fastq -2 $outputPath${prefix}_2.clean.fastq \
        -S $outputPath$prefix.sam --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 40
#samtools
samtools view -bS $outputPath$prefix.sam > $outputPath$prefix.bam
samtools sort $outputPath$prefix.bam -o $outputPath$prefix.sort.bam
samtools index $outputPath$prefix.sort.bam

echo exit status: $?

echo Variant Calling  ending at:
date
