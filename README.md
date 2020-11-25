# Best Practice for Germline Variants (SNP/InDel) discovery AND functional annotation workflow
To find short gene mutations using whole exome sequencing, it is recommend to follow the workflow of GATK4 and Annovar. [GATK4](https://gatk.broadinstitute.org/hc/en-us) is a famous widely-used tools developed by the Broad Instituite for detecting mutations, and it holds a indusctry standard. Annovar can help you to annotate a variety of databases onto your VCF files. 

**Download and install GATK4 and Annovar**
-------------------------------------------
Latest [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) Download:

http://download.openbioinformatics.org/annovar_download_form.php

Install GATK4 by conda

```
conda install GATK4
```

The latest GATK4 support VQSR variant correction based on CNN method, so it requires tensorflow as dependencies, which if we want to use CNN correction methods, we need to create a conda environment for GATK4 only follow by the special yaml file. The link explain how to do it indicates as follow. This GATK4 environment cannot share with other softwares otherwise it will collapse.

https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4

Create special enviroment for GATK4 after you downlaod and unzip the GATK file

```
conda env create -n gatk -f gatkcondaenv.yml
```

**Download Resource bundle**
-----------------------------

The resource bundle includes reference genomes hg19/hg38, and some well-known databases for storing SNPs/InDel variants information. It is recommend to download these resource, because we will use these databases for following analysis. We will treat these databases as valid and correct reference, and use them to correct errors occuring in our analysis. 

You could choose different resource bundle depending on your data. Typically, hg38 resource bundle is recommend. However, if your sequence library preparation kits are based on hg19 reference genome, you need to use same genome, as the coordinates are different which will affect PCR regions. This situation should be noticed in whole exome sequencing. 

hg38:

https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0

hg19:

https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37

https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0

**STEP1: QC of raw reads**
----------------------------------

For Illumina short reads, it is usually following workflow of Fastqc and raw reads cleaning steps. Fastqc allows you to visualized your reads quality. Typically, we want our reads quality to be higher than Q30. Then, we do reads trim. I use [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [fastp](https://github.com/OpenGene/fastp) to trim my raw reads into clean reads.

```
trimmomatic PE -phred33 -threads 6 -quiet \
               -validatePairs read_1.fq read_2.fq \
                              read_1.paired.fq read_1.unpaired.fq \
                              read_2.paired.fq read_2.paired.fq \
                              LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36
```

```
fastp -i read_1.fq -o read_1.clean.fastq -I read_2.fq -O read_2.clean.fastq \
      -A -x --cut_mean_quality 20 --cut_front_window_size 4 --cut_front_mean_quality 20 \
      --cut_tail_window_size 4 --cut_tail_mean_quality 20 -n 5 --cut_front 20 \
      -f 15 -t 5 -F 15 -T 5
```

**STEP2: Reads Mapping and Sort bam file**
----------------------------------
If you use [bwa](https://github.com/lh3/bwa) or [bwa-mem2](https://github.com/lh3/bwa), you need to index you reference file first. 

```
bwa index reference.fasta
```
OR
```
bwa-mem2 index reference.fasta
```
Then, you can map your clean reads onto the reference genome.
```
bwa mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta read_1.clean.fastq read_2.clean.fastq > sample.sam
```
OR
```
bwa-mem2 mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta read_1.clean.fastq read_2.clean.fastq > sample.sam
```
```
samtools view -bS sample.sam > sample.bam
samtools sort sample.bam -o sample.sort.bam
samtools index sample.sort.bam
```

**Germline Variant Call using Haplotype**
-----------------------------------------
There are several steps need to be done before we call the muations. As we need to eliminate PCR containmination, it will introduce PCR duplicates reads that will affect the calling model based on sequencing depth. For example, a false SNP cuased by sequencing error should has very shallow sequencing depth like one or two, but if this SNP was duplicateed multiple time during PCR process, the depth is increased, and this false SNP can be recogenised as a true SNP because of relatively high coverage. 


1. MarkDuplicates and Revome PCR effects
```
gatk MarkDuplicates -I sample.sort.bam -M sample.dupMetrix.txt -O sample.rmdup.bam --REMOVE_DUPLICATES true
```

2. Add Read Group

The second step is tricky, the GATK would not work if the bam file does not contain read group information. The read group is related with biological duplicates issues, by using read group tag, the reads from different biological duplicates can be identified. GATK expects all read groups appearing in the read data to be specified in the file header, and will fail with an error if it does not find that information (whether there is no read group information in the file, or a subset of reads do not have read groups).

```
gatk AddOrReplaceReadGroups -I sample.rmdup.bam -O sample.rmdup.rg.bam -LB Solexa-272222 -PL illumina -PU H0164ALXX140820.2 -SM $sample
```

FAQs from GATK official website

https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information

3. Create Sequence Dict file

The reference fasta dictionary file should be create. This file is a input file for the downstream input. This file just can be created only once.
```
gatk CreateSequenceDictionary -R hg19_v0_Homo_sapiens_assembly19.fasta -O hg19_v0_Homo_sapiens_assembly19.dict
```

4. Base Recalibration

Base quality score recalibration (BQSR) is a process in which is the application of machine learning to model these errors empirically and adjust the quality scores accordingly.For example we can identify that, for a given run, whenever we called two A nucleotides in a row, the next base we called had a 1% higher rate of error. So any base call that comes after AA in a read should have its quality score reduced by 1%. We do that over several different covariates (mainly sequence context and position in read, or cycle) in a way that is additive. So the same base may have its quality score increased for one reason and decreased for another.

https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-

```
gatk BaseRecalibrator \
  -I sample.rmdup.rg.bam \
  -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
  --known-sites /path/to/datasets/hg19_v0_dbsnp_138.b37.vcf.gz  \
  --known-sites /path/to/datasets/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  -O sample.bqsr.grp
```

The input BAM file should be the mark and remove the duplicates from the original bam file. This step is not apply BQSR on the bam file, and it just record the BQSR information in the file. The high confidence datasets are used to be as input for training information of the programme.

```
gatk ApplyBQSR \
  -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
  -I sample.rmdup.rg.bam \
  --bqsr-recal-file sample.bqsr.grp \
  -O sample.bqsr.bam
```

Apply BQSR on bam file, using the last step bqsr recall file as input to correct the error.

5. HaplotypeCaller

This is the main step calling variants from the bam file. The output is in vcf format (https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format) which record SNP/InDel.

```
gatk HaplotypeCaller \
  -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
  -L exon_hg19_v0.interval_list \
  -I sample.bqsr.bam \
  -O sample.vcf.gz
```
The -L option takes either BED file or Picard interval file. One thing should be noticed that the Interval file must be coordinated with reference FASTA file. The most common problem is that the 'chr' prefix before chrom field, e.g chr1 vs 1. The chrom field must be consistent in all the way.

6. Variant Filter

The latest GATK version4 combined CNN algorithm to score each mutations based on high quality reference databases such as 1000 Genome. This step will assign a quality score to SNP/InDel.
```
gatk CNNScoreVariants \
   -I sample.bqsr.bam \
   -V sample.vcf.gz \
   -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
   -O sample.CNN.vcf \
   --tensor-type read_tensor
```

Then, using CNN_2D to filter out the mutations based on assigned score in the previous step.
```
gatk FilterVariantTranches \
   -V sample.CNN.vcf \
   --resource /path/tp/reference/hg19_v0_1000G_omni2.5.b37.vcf.gz \
   --resource /path/to/reference/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf.gz \
   --resource /path/to/reference/hg19_v0_Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
   --info-key CNN_2D \
   --snp-tranche 99.95 \
   --indel-tranche 99.4 \
   --invalidate-previous-filters \
   -O $path$name.CNNfiltered.vcf
```

7. Calculating genotype Posteriors

Genotype is important in genetics analysis, especially in genetics consultation field. When analyze family trio samples, genotype probability should be correct.

```
gatk CalculateGenotypePosteriors \
    -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
    -V sample.CNNfiltered.vcf \
    -supporting /path/to/reference/hg19_v0_1000G_omni2.5.b37.vcf.gz \
    -supporting /path/to/reference/hg19_v0_hapmap_3.3.b37.vcf.gz \
    -O sample.genotypeCorrected.vcf
```

8. Hard filtering 

If you want to filter again, another choice is to do hard filtering. Just like the old version of GATK, hard filtering is popular.

```
gatk VariantFiltration \
   -R /path/to/reference/hg19_v0_Homo_sapiens_assembly19.fasta \
   -V sample.genotypeCorrected.vcf \
   -O sample.prefiltered.vcf \
   -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name 'FAIL' \
   --genotype-filter-expression 'GQ <20' \
   --genotype-filter-name 'FAIL'
```

**Annotation**
---------------------
There are some popular annotaiton scripts such as ANNOVAR, Snpneff and so forth. In my pipepline, I use ANNOVAR (Wang et al., 2010) 

**Other Issues in mutations detection**
--------------------------------------

[Haplotype phasing](http://data-science-sequencing.github.io/Win2018/lectures/lecture10/#:~:text=Haplotype%20phasing%20is%20the%20problem,problem%2C%20there%20are%20many%20methods.)

Sequencing depth

[Variant Interpretation](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)


**Citations**
-----------------------------------------

Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010
