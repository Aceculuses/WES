gatk BedToIntervalList 
    -I exon_sorted_nochr_prefix.bed \
    -O exon_hg19_v0.interval_list \
    -SD /data2/guoping/Databases/Broad_bundle_hg19/hg19_v0_Homo_sapiens_assembly19.dict
