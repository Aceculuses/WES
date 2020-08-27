FQDIR="/data2/Data/200821_A00249_0074_AHMMHNDRXX/tumor_wes/"
SAMPLES,=glob_wildcards(FQDIR+"{sample}_R1_001.fastq.gz")

GENOME="/data2/database/DNA_Panel/biodata/hg1kv37/human_g1k_v37.fasta"
BED="/data2/Data/exon_sorted_nochr_prefix.bed"

rule all:
    input:
        expand('anno/{sample}.anno.filter.csv', sample=SAMPLES),
        expand("tmb/{sample}.tmb", sample=SAMPLES)

rule fastp_pe:
    input:
        #sample=expand("/data2/Data/tumor/{sample}_R{rep}_001.fastq.gz", sample=SAMPLES, rep=[1,2])
        fq=[FQDIR+"{sample}_R1_001.fastq.gz",FQDIR+"{sample}_R2_001.fastq.gz"]
    output:
        trimmed=["trimmed/pe/{sample}.1.fastq", "trimmed/pe/{sample}.2.fastq"],
        html="report/{sample}.html",
        json="report/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra=""
    threads:10
    shell:
        "/data2/software/fastp -i {input.fq[0]} -o {output.trimmed[0]} -I {input.fq[1]} -O {output.trimmed[1]} -j {output.json} -h {output.html} -w {threads} 1>{log} 2>&1"
rule bwa_mem:
    input:
        reads=["trimmed/pe/{sample}.1.fastq", "trimmed/pe/{sample}.2.fastq"]
    output:
        temp("mapped/{sample}.sam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=GENOME,
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads:10
    shell:
        "/data2/software/bwa-0.7.12/bwa mem -M {params.extra} -t {threads} -K 10000000 -P {params.index} {input} > {output} 2>{log}"

rule samtools_sort:
    input:
        "mapped/{sample}.sam"
    output:
        "mapped/{sample}.sort.bam"
    log:
        "logs/{sample}.sort.log"
    threads:10
    shell:
        "/data2/software/samtools-1.2/samtools view -bSh {input} -@ {threads} | /data2/software/samtools-1.2/samtools sort  -o {output} -@ {threads} 2>{log}"

rule mark_duplicates:
    input:
        "mapped/{sample}.sort.bam"
    output:
        bam="dedup/{sample}.dedup.bam",
        metrics="dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        java_opt="-Xmx30g",
        extra="REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"
    shell:
        "java {params.java_opt} -jar /data2/software/picard-tools-1.119/MarkDuplicates.jar I={input} O={output.bam}  {params.extra} METRICS_FILE={output.metrics} 2>{log}"

rule gatk_baserecalibrator:
    input:
        bam="dedup/{sample}.dedup.bam",
        ref=GENOME,
        known="/data2/software/muTect-1.1.7/dbsnp_132_b37.leftAligned.vcf"  # optional known sites
    output:
        recal_table="recal/{sample}.grp"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    params:
        extra="--known-sites /data2/database/DNA_Panel/biodata/bundle_b37/Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites /data2/database/DNA_Panel/biodata/bundle_b37/1000G_phase1.indels.b37.vcf",  # optional
        java_opts="--java-options '-Xmx30g'", # optional
    shell:
        "/data2/software/gatk-4.1.8.0/gatk  {params.java_opts} BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known} {params.extra} -O {output.recal_table} 2>{log}"

rule gatk_applybqsr:
    input:
        bam="dedup/{sample}.dedup.bam",
        ref = GENOME,
        recal_table="recal/{sample}.grp"
    output:
        bam=protected("recal/{sample}.bqsr.bam")
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="--java-options '-Xmx30g'", # optional
    shell:
        "/data2/software/gatk-4.1.8.0/gatk {params.java_opts} ApplyBQSR -R {input.ref}  -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam} 2>{log}"

rule mutect2:
    input:
        fasta = GENOME,
        map = "recal/{sample}.bqsr.bam",
        gnomad = "/data2/database/af-only-gnomad.raw.sites.b37.vcf.gz",
        pon = "/data2/zhangsy/PON/pon.vcf.gz"
    output:
        vcf = "variant/{sample}.vcf.gz"
    message:
        "Testing Mutect2 with {wildcards.sample}"
    log:
        "logs/mutect_{sample}.log"
    params:
        extra=BED,  # optional
        java_opts="--java-options '-Xmx30g'", # optional
    shell:
         "/data2/software/gatk-4.1.8.0/gatk {params.java_opts} Mutect2 -R {input.fasta} -I {input.map} -L {params.extra} --germline-resource {input.gnomad}  --panel-of-normals {input.pon} -O {output.vcf} 2>{log}"

rule getPileupSummay:
    input:
        "recal/{sample}.bqsr.bam"
    output:
        "contamination/{sample}_getpileupsummaries.table"
    params:
        exac_vcf="/data2/database/small_exac_common_3_b37.vcf.gz",
        interval=BED,
        java_opt='''--java-options "-Xmx30g" '''
    log:
        "logs/contamination/{sample}.getPileup.log"
    shell:
        "/data2/software/gatk-4.1.8.0/gatk {params.java_opt} GetPileupSummaries -I {input} -V {params.exac_vcf}  -L {params.interval}  -O {output} 2>{log}"

rule calContamination:
    input:
        "contamination/{sample}_getpileupsummaries.table"
    output:
        "contamination/{sample}_tumor_calculatecontamination.table"
    log:
        "logs/contamination/{sample}.calcontam.log"
    shell:
        "/data2/software/gatk-4.1.8.0/gatk CalculateContamination  -I {input} -O {output} 2>{log}"

rule filterMutectCall:
    input:
        vcf="variant/{sample}.vcf.gz",
        contamination_tab="contamination/{sample}_tumor_calculatecontamination.table"
    output:
        "variant/{sample}_somatic_oncefiltered.vcf.gz"
    params:
        ref= GENOME
    log:
        "logs/{sample}.filterMutectCall.log"
    shell:
        "/data2/software/gatk-4.1.8.0/gatk FilterMutectCalls  -R {params.ref} -V {input.vcf} --contamination-table {input.contamination_tab}  -O {output} 2>{log}"

rule filterPassVar:
    input:
        "variant/{sample}_somatic_oncefiltered.vcf.gz"
    output:
        "variant/{sample}.pass.vcf"
    run:
        import gzip
        with open(output[0], "w") as out_f:
            with gzip.open(input[0], "rb") as in_f:
                for line in in_f.readlines():
                    if line.startswith("#".encode()):
                        out_f.write(line.decode())
                    else:
                        if 'PASS'.encode() in line:
                            info = line.decode().split('\t')[-1]
                            af = info.split(':')[2]
                            #print(float(af))
                            if float(af) > 0.05:
                                out_f.write(line.decode())

rule annovar:
    input:
        vcf="variant/{sample}.pass.vcf",
        humandb="/data2/database/DNA_Panel/biodata/humandb"
    output:
        "anno/{sample}.anno.hg19_multianno.txt"
    log:
        "logs/annovar/{sample}.log"
    params:
        genome_ver="hg19",
        protocol="-protocol refGene,1000g2015aug_all,dbnsfp30a,cosmic74,clinvar_20200316,intervar_20180118,avsnp142 -operation g,f,f,f,f,f,f"
    shell:
        "perl /data2/software/annovar/table_annovar.pl {input.vcf} {input.humandb} -buildver {params.genome_ver} -out anno/{wildcards.sample}.anno -remove {params.protocol} -nastring . -vcfinput -polish 2>{log}"

rule filterAnno:
    input:
        "anno/{sample}.anno.hg19_multianno.txt"
    output:
        "anno/{sample}.anno.filter.txt"
    run:
        with open(output[0], "w") as out_f:
            with open(input[0], "r") as in_f:
                for line in in_f.readlines():
                    item_l = line.split('\t')
                    if line.startswith('Chr'):
                        out_f.write(line)
                        func_idx = item_l.index('Func.refGene')
                        exonic_idx = item_l.index('ExonicFunc.refGene')
                    else:
                        func_l = ['exonic', 'splicing', 'exonic;splicing']
                        if item_l[func_idx] in func_l and item_l[exonic_idx] != 'synonymous SNV':
                            out_f.write(line)

rule txt2csv:
    input:
        "anno/{sample}.anno.filter.txt"
    output:
        "anno/{sample}.anno.filter.csv"
    run:
        import csv
        import sys
        import os

        try:
            with open(input[0], 'r') as reader:
                list_data = reader.readlines()
                columns = list_data[0].rstrip().split('\t')
                list = []
                for i in list_data [1:]:
                    list.append(i.rstrip().split('\t'))
        except IOError as e:
            print(e)
            sys.exit(1)

        with open(output[0],"w") as csvfile:
            writer = csv.writer(csvfile)
            #先写入columns_name
            writer.writerow(columns)
            #写入多行用writerows
            writer.writerows(list)
rule tmb:
    input:
        "anno/{sample}.anno.filter.txt"
    output:
        "tmb/{sample}.tmb"
    run:
        n = 0
        with open(input[0], "r") as in_f:
            for line in in_f.readlines():
                item_l = line.split('\t')
                if line.startswith("Chr"):
#                    out_f.write(line)
                    func_idx = item_l.index('Func.refGene')
                    exonic_idx = item_l.index('ExonicFunc.refGene')
                else:
                    filter_l = ['nonframeshift', 'unknown']
                    af = item_l[-1].split(':')[2]
                    if item_l[exonic_idx] not in filter_l and float(af) >= 0.05:
                        n += 1

        with open(output[0], "w") as out_f:
            tmb = float(n/30)
            out_f.write(str(tmb))
