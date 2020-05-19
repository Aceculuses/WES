#!/usr/bin/python
import subprocess
import sys
    
def indexing(txt):
    container=[]
    index={}
    header=open(txt,'r').readlines()[0].rstrip('\n')
    headerL=list(header.split('\t'))
    container.append(headerL)
    for i in headerL:
        index[i]=headerL.index(i)+1
    container.append(index)
    return container
    
def subfileIndex(container,header):
    subfileindex={}
    overlape=set(container[0]).intersection(set(header))
    for i in overlape:
        subfileindex[i]=container[1][i]
    return subfileindex

def output(subfileindex,cmd_Output,txt):
    pidd=set()
    cmd=''.join(cmd_Output)
    for i in subfileindex.keys():
        command=cmd.format(txt,subfileindex[i],i)
#        print('command :',command)
        p=subprocess.Popen(command,shell=True)
        while p.poll() ==None:
            pidd.add(p.pid)
    
def clinsig(index,clinvar,clinvarcol,cmd_Output,txt,outputPath):
    clinSig=open(outputPath+'clinsig.o','w')
    clindict={}
    clinfo=open(''.join(clinvar),'r').readlines()
    for i in clinfo:
        clinline=i.rstrip('\n').split('\t')
        clindict[clinline[0]]=clinline[-1]

    clindate=''.join(clinvarcol)
    cmd=''.join(cmd_Output).format(txt,index[clindate],clindate)
#    print(cmd)
    p=subprocess.Popen(cmd,shell=True)
    while p.poll() ==None:
        p.pid
    
    clinlist=[]
    clinOpen=open(outputPath+clindate+'.o','r').readlines()
    for i in clinOpen:
        x=i.rstrip('\n').split('\t')
        clinlist.append(x[0])
    
    for i in clinlist:
        if i=='.':
            clinSig.write('.'+'\n')
        elif i == clindate:
            clinSig.write('clinsig'+'\n')
        else:
            clinSig.write(clindict[i]+'\n')

    numOfL=len(clinlist)
    return numOfL

def meiyou(na,numOfL,outputPath):
    for i in na:
        file=open(outputPath+i+'.o','w')
        for x in range(0,numOfL):
            if x ==0:
                file.write(i+'\n')
            else:
                file.write('.'+'\n')

def GT(txt,cmd_Output,index,cmd_sed,outputPath):
    cmd1=''.join(cmd_Output).format(txt,index['Otherinfo13'],'Otherinfo13.tmp')
    cmd2=''.join(cmd_sed).format(outputPath+'Otherinfo13.tmp.o',index['Otherinfo13'],'Otherinfo13')
    p1=subprocess.Popen(cmd1,shell=True)
    while p1.poll()==None:
        p1.pid
    
    p2=subprocess.Popen(cmd2,shell=True)
    while p2.poll()==None:
        p2.pid
    

def merge(header,prefix,inputPath,outputPath):
    for i in range(len(header)):
        header[i]=outputPath+header[i]+'.o'
    filename=' '.join(header)
    cmd1='paste '+filename+' > ' +outputPath+prefix+'_snp_indel.txt'
#    print(cmd1)
    p1=subprocess.Popen(cmd1,shell=True)
    while p1.poll()==None:
        p1.pid
   
    cmd2='rm '+outputPath+'*.o'
    p2=subprocess.Popen(cmd2,shell=True)
    while p2.poll()==None:
        p2.pid

def transfer(outputPath,pathFortxt,pathForcsv,cmd_mv,prefix,cmd_csv,cmd_nospace):
    cmd1=''.join(cmd_csv).format(outputPath,prefix,pathForcsv,prefix)
    cmd2=''.join(cmd_mv).format(outputPath,prefix,pathFortxt)
    cmd3=''.join(cmd_nospace).format(pathFortxt,prefix,pathFortxt,prefix)
#    print(cmd2)
    p1=subprocess.Popen(cmd1,shell=True)
    while p1.poll()==None:
        p1.pid
    p2=subprocess.Popen(cmd2,shell=True)
    while p2.poll()==None:
        p2.pid
    p3=subprocess.Popen(cmd3,shell=True)
    while p3.poll()==None:
        p3.pid
   

if __name__=='__main__':
#   pathForReformat
    inputPath=sys.argv[1]
#   pathForOutput
    outputPath=sys.argv[2]
#   pathFortxt
    pathFortxt=outputPath+'txt/'
#   pathForcsv
    pathForcsv=outputPath+'csv/'

    genomeVersion=sys.argv[3]
    clinvar_date=sys.argv[4]
    file=sys.argv[5]
#    print(file)
#    txt='CS14.hg19_multianno.txt'
    prefix=file.split('.')[0]
    print('The file prefix is: '+prefix)
    txt=inputPath+file
#    print(txt)

    if genomeVersion == 'hg19':
        header=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene',
                'Ensembl_ID','snp142','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all',
                'Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
                'gnomAD_exome_EAS','Otherinfo13','clinsig','Func.refGene']
        na=['Ensembl_ID','Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
          'gnomAD_exome_EAS']
    else:
        header=['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene',
                'ensembl_ID','avsnp150','GeneDetail.refGene','ExAC_ALL','esp6500siv2_all',
                'Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
                'gnomAD_exome_EAS','Otherinfo13','clinsig','Func.refGene']
        na=['ensembl_ID','Freq_1000g2015aug_all','Freq_1000g2015aug_EAS','gnomAD_exome_ALL',
          'gnomAD_exome_EAS']
    
    clinvar=['/home/ubuntu/software/annovar/humandb/clinvarinfo_db/',genomeVersion,'_clinvarinfo_',clinvar_date]
    
    clinvarcol=['clinvar_',clinvar_date]

    cmd_Output=['cat ','{} ','| ','cut -f','{} ','> ',outputPath,'{}','.o']
    cmd_sed=['cat ','{} ','| ','cut -f','{} ',"| sed 's#:#\t#g' ","| sed 's#^# #g' ","| sed 's# Otherinfo13#GT\tAD\tDP\tGQ\tPL#g' ",'> ',outputPath,'{}.o']
    cmd_mv=['mv ','{}',"{}* ",'{}']
    cmd_csv=['cat ','{}{}* ','| ',"sed 's#,# #g' ","| ","sed 's#\t#,#g' ",'> ','{}{}_snp_indel.csv']    
    cmd_nospace=['cat ','{}{}* ','| ',"sed 's# ##g' ",'> ','{}{}_snp_indel.final.txt']

    container=indexing(txt)

    output(subfileIndex(container,header),cmd_Output,txt)

    numOfL=clinsig(container[1],clinvar,clinvarcol,cmd_Output,txt,outputPath) 

    meiyou(na,numOfL,outputPath)
    
    GT(txt,cmd_Output,container[1],cmd_sed,outputPath)

    merge(header,prefix,inputPath,outputPath) 
    
    transfer(outputPath,pathFortxt,pathForcsv,cmd_mv,prefix,cmd_csv,cmd_nospace)
