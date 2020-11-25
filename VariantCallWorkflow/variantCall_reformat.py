#!/data2/guoping/software/miniconda3/envs/gatk/bin/python3
import subprocess
import sys
import time

def readInfoTable(inputInfoTable):
    dict={}
    file=open(inputInfoTable,'r').readlines()
    for i in file:
        line=i.rstrip('\n').split(':')
        dict[line[0]]=line[1]
    return dict
    
def exportIndexList(inputFile):
    fi=[]
    file=open(inputFile,'r').readlines()
    for i in file:
        NR=i.rstrip('\n').split('\t')
        fi.append(NR)
    index={}
    header=fi[0]
    for NF in header:
        index[NF]=header.index(NF)+1
    return index

def splitInputFileByHeader(inputFile,indexDict,header,inputPath):
    for i in header:
        cmd='cat '+inputFile+' | '+'cut -f'+str(indexDict[i])+' > '+inputPath+i+'.o'
        subprocess.Popen(cmd,shell=True)
        time.sleep(0.5)
        print(cmd)

def turnClinvarNumberIntoClinsig(clinvarNumber,clinvarDate,inputPath):
    clinvarNumberList=[]
    file=open(clinvarNumber,'r').readlines()
    outputFile=open(inputPath+'clinsig_'+clinvarDate+'.o','w')
    clinDict={}
    for i in file:
        number=i.rstrip('\n').split(' ')
        clinvarNumberList.append(number[0])

    clinInfo='/data2/guoping/software/annovar/humandb/clinvarinfo_db/hg19_clinvarinfo_'+str(clinvarDate)
    
    file2=open(clinInfo,'r',encoding='utf-8')
    for i in file2.readlines():
        i=i.encode('utf-8')
        x=str(i).split('\\t')
        number="".join(filter(str.isdigit,x[0]))
#        print(number)
        clinsig=x[-1].replace("\\n'",'')
        clinDict[int(number)]=clinsig
#    print(clinDict) 
    for i in clinvarNumberList:
        if i.isdigit() == False:
            outputFile.write(i+'\n')
        else:
            outputFile.write(clinDict[int(i)]+'\n')

def turnOtherinfo13(inputPath):
    allTheLine=[]
    file=open(inputPath+'Otherinfo13.o','r').readlines()
    fileTurned=open(inputPath+'Otherinfo13.o','w')
    for i in file:
        line=i.rstrip('\n').split(':')
        if len(line)==6:
            line.insert(3,'PASS')
            allTheLine.append('\t'.join(line))
        else:
            allTheLine.append('\t'.join(line))
    allTheLine[0]='GT\tAD\tDP\tFT\tGQ\tPL\tPP'
#    print(allTheLine)
    for i in allTheLine:
        fileTurned.write(i+'\n')
    return(len(allTheLine))
        
def createNewTag(NumberOfallTheLine,noThisField,inputPath):
    allTheMissFile=[]
    for i in noThisField:
        missField=['.']*NumberOfallTheLine
        missField[0]=i
        f=open(inputPath+i+'.o','w')
        for t in missField:
            f.write(t+'\n')
        missField=None

def concatateFile(template,inputPath,inputFile):
    newTemplate=[]
    for i in template:
        i=inputPath+i+'.o'
        newTemplate.append(i)
    string=' '.join(newTemplate)
    cmd='paste '+string+' > '+inputPath+inputFile.replace('.hg19_multianno.txt','')+'_snp_indel.txt'
    print(cmd)
    subprocess.Popen(cmd,shell=True) 
                                 
if __name__=='__main__':
    inputInfoTable=sys.argv[1]
    infoTableDict=readInfoTable(inputInfoTable)
    inputFile=infoTableDict['file']
    inputPath=infoTableDict['inputPath']
    template=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene',
            'Ensembl_ID','snp142','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all',
            'Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
            'gnomAD_exome_EAS','Otherinfo13','clinsig_'+infoTableDict['clinvar_date'],'Func.refGene']
    noThisField=['Ensembl_ID','Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
          'gnomAD_exome_EAS']
     
#    print(infoTableDict)
    header=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene',
            'snp142','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all','Otherinfo13','clinvar_'+infoTableDict['clinvar_date'],'Func.refGene']

    indexDict=exportIndexList(inputPath+inputFile)
    splitInputFileByHeader(inputPath+inputFile,indexDict,header,inputPath)
    
    clinvarNumber=inputPath+'clinvar_'+str(infoTableDict['clinvar_date'])+'.o'

    turnClinvarNumberIntoClinsig(clinvarNumber,infoTableDict['clinvar_date'],inputPath)
    
    NumberOfallTheLine=turnOtherinfo13(inputPath)

    createNewTag(NumberOfallTheLine,noThisField,inputPath)
 
    concatateFile(template,inputPath,inputFile)
