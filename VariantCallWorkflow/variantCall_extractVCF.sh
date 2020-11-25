input=$1
path=$2

gatk IndexFeatureFile -I $path$1.final.vcf


gatk SelectVariants \
     -R /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.fasta \
     -V $input.final.vcf \
     -L /home/guoping/bed_file/PAH.interval_list \
     -O $path$input.PAH.vcf
