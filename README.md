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

**STEP2: Reads Mapping**
----------------------------------
If you use [bwa](https://github.com/lh3/bwa) or [bwa-mem2](https://github.com/lh3/bwa), you need to index you reference file first. 

```
bwa index reference.fasta
```
OR
```
bwa-mem2 index reference.fasta
```
Then, you can map your clean reads onto the reference genome
```
bwa mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta read_1.clean.fastq read_2.clean.fastq > sample.sam
```
OR
```
bwa-mem2 mem Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta read_1.clean.fastq read_2.clean.fastq > sample.sam
```

**Other Issues in mutations detection**
--------------------------------------

[Haplotype phasing](http://data-science-sequencing.github.io/Win2018/lectures/lecture10/#:~:text=Haplotype%20phasing%20is%20the%20problem,problem%2C%20there%20are%20many%20methods.)



