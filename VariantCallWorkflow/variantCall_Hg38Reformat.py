#!/data2/guoping/software/miniconda3/envs/ngs/bin/python3
import sys
import subprocess
import time

def getClinvarDate(file):
    cmd="head -n1 " +file+" | tr '\t' '\n' | grep clinvar"
    clinvar=subprocess.check_output(cmd,shell=True)
    clinvar_title=clinvar.decode('utf-8').rstrip('\n')
    clinvar_date="".join(filter(str.isdigit, clinvar_title))
#    print(clinvar_date)
    return clinvar_date

def getFileHeaderIndex(file):
    cmd="head -n1 " +file
    head=subprocess.check_output(cmd,shell=True)
    FileHead=head.decode('utf-8').rstrip('\n').split('\t')
    IndexDict={}
    for i in FileHead:
        IndexDict[i]=int(FileHead.index(i))+1
    return IndexDict
    
def Clinsig(clinvar_date):
    ClininfoPath='/data2/guoping/software/annovar/humandb/clinvarinfo_db/'
    infoFile='hg38_clinvarinfo_'+clinvar_date
    file=open('/data2/guoping/software/annovar/humandb/clinvarinfo_db/'+infoFile,'r').readlines()
    Clinvar_Dict={}
    for i in file:
        x=i.rstrip('\n').split('\t')
        Clinvar_Dict[x[0]]=x[-1]
    return Clinvar_Dict

def getClinvar(clinvar_date,IndexDict,Clinvar_Dict,file):
    clinvar_index=IndexDict['clinvar_'+clinvar_date]
    cmd='cut -f'+str(clinvar_index)+' '+str(file)
    clinsig=subprocess.check_output(cmd,shell=True)
    clinsig=clinsig.decode('utf-8').rstrip('\n').split('\n')
    output=open('clinsig.hg38.o','w')
    output.write('clinsig'+'\n')
    for i in clinsig:
        if i == '.':
            output.write(i+'\n')
        elif i =='clinvar_'+clinvar_date:
            pass
        else:
            output.write(Clinvar_Dict[i]+'\n')
        
def getOtherinfo13(IndexDict,file):
    otherinfo13Index=IndexDict['Otherinfo13']
    cmd='cut -f'+str(otherinfo13Index)+' '+str(file)
    otherinfo13=subprocess.check_output(cmd,shell=True)
    otherinfo13=otherinfo13.decode('utf-8').rstrip('\n').split('\n')
    output=open('Otherinfo13.hg38.o','w')
    output.write('GT'+'\t'+'AD'+'\t'+'DP'+'\t'+'FT'+'\t'+'GQ'+'\t'+'PL'+'\t'+'PP'+'\n')
    for i in otherinfo13:
        x=i.split(':')
#        print(x)
        if len(x) == 7:
            output.write(x[0]+'\t'+x[1]+'\t'+x[2]+'\t'+x[3]+'\t'+x[4]+'\t'+x[5]+'\t'+x[6]+'\n')
        elif len(x) ==6:
            output.write(x[0]+'\t'+x[1]+'\t'+x[2]+'\t'+x[3]+'\t'+'PASS'+'\t'+x[4]+'\t'+x[5]+'\n')
        else:
            pass

def getRest(Rest,file,IndexDict):
    for i in Rest:
        cmd='cut -f'+str(IndexDict[i])+' '+str(file)+' '+'>'+' '+i+'.hg38.o'
        subprocess.check_output(cmd,shell=True)
        time.sleep(1)

def getNA(NA,file):
    cmd='wc -l'+' '+file
    line=subprocess.check_output(cmd,shell=True)
    lineNum=line.decode('utf-8').rstrip('\n').split(' ')
    for i in NA:
        output=open(i+'.hg38.o','w')
        output.write(i+'\n')
        for i in range(1,int(lineNum[0])):
            output.write('.'+'\n')
        output.close() 

def concate(template,file):
    x=''
    for i in template:
        x=x+i+'.hg38.o '
    Name=str(file).split('.')[0]
    cmd='paste '+x+' > '+Name+'.hg38_snp_indel.txt'
    subprocess.check_output(cmd,shell=True)
        
if __name__=="__main__":
    file=sys.argv[1]
    template=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene',
            'ensembl_ID','avsnp150','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all',
            'Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL','gnomAD_exome_EAS',
            'Otherinfo13','clinsig','Func.refGene']
    NA=['ensembl_ID','Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL','gnomAD_exome_EAS']

    Rest=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene','avsnp150','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all','Func.refGene']

    clinvar_date=getClinvarDate(file)
    Clinvar_Dict=Clinsig(clinvar_date)
    IndexDict=getFileHeaderIndex(file)

    getClinvar(clinvar_date,IndexDict,Clinvar_Dict,file)
    getOtherinfo13(IndexDict,file)
    getRest(Rest,file,IndexDict)
    getNA(NA,file)
    concate(template,file)
