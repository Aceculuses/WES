# Best Practice for Germline Variants (SNP/InDel) discovery AND functional annotation
To find short gene mutations using whole exome sequencing, it is recommend to follow the workflow of GATK4 and Annovar. [GATK4](https://gatk.broadinstitute.org/hc/en-us) is a famous widely-used tools developed by the Broad Instituite for detecting mutations, and it holds a indusctry standard. Annovar can help you to annotate a variety of databases onto your VCF files. 

**Download and install GATK4 and Annovar**

**Download Resource bundle**

The resource bundle includes reference genomes hg19/hg38, and some well-known databases for storing SNPs/InDel variants information. It is recommend to download these resource, because we will use these databases for following analysis. We will treat these databases as valid and correct reference, and use them to correct errors occur in our analysis. You could choose different resource bundle depending on your data. Typically, hg38 resource bundle is recommend. However, if your sequence library preparation kits are based on hg19 reference genome, you need to use same genome.

https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0

https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37
