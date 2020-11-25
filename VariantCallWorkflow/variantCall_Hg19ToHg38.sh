echo VariantCall_Hg19toHg38 starting at:
date

inputPath=$1
inputVCFprefix=$2

CrossMap.py vcf /data2/guoping/software/annovar/humandb/GRCh37_to_GRCh38.chain.gz $inputPath/$inputVCFprefix.final.vcf \
                /data2/guoping/software/annovar/humandb/Homo_sapiens.GRCh38.dna.toplevel.fa $inputPath/$inputVCFprefix.hg38.vcf

convert2annovar.pl -format vcf4 $inputPath/$inputVCFprefix.hg38.vcf -outfile $inputPath/$inputVCFprefix.hg38.avinput --withfreq -includeinfo

table_annovar.pl $inputPath/$inputVCFprefix.hg38.avinput /data2/guoping/software/annovar/humandb -buildver hg38 -out $inputPath/$inputVCFprefix \
                         -remove --protocol refGene,ensGene,exac03,clinvar_20200113,esp6500siv2_all,avsnp150,dbnsfp35c,1000g2015aug,gnomad211_exome \
                         -operation g,g,f,f,f,f,f,f,f \
                         -nastring . \
                         --polish \
                         -xreffile ~/software/annovar/humandb/hg38_refGene.txt \
                         --otherinfo \
                         --thread 32

python /home/guoping/scripts/wes/variantCall_Hg38Reformat.py $inputPath/$inputVCFprefix.hg38_multianno.txt

echo exit status: $?

echo VariantCall_Hg19toHg38 ending at:
date
