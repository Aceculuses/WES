#!/usr/bin/python
import sys
import subprocess
import os
import time
import getopt
import numpy as np
import datetime

#--------------------------------------------------------------------------------------
#    ***The current version is Version 4***
#    ***The control flow has been added to manage multi-subprocess***
#    ***Implement Object Orientation Programming
#    ***Author: Guoping Gu
#    ***Date: 5/14/2020
#--------------------------------------------------------------------------------------

#create folders
def mkdir(int,list,str1):
    index=int
    pathList=list
    genomeVersion=str1
    if os.path.exists(pathList[index]):
        print(pathForResults+'\033[1;35m'+' Exists!'+"\033[0m")
        print('Use another folder name next time. PROGRAMME\033[1;31m EXITS\033[0m')
        sys.exit()

    elif genomeVersion=='hg19':
        for dir in pathList[0:8]:
            os.makedirs(dir)
        print('Generate analysis folders.......\033[1;32mDone\033[0m')

    elif genomeVersion=='hg38':
        for dir in pathList:
            os.makedirs(dir)
#    else:
#        for dir in pathList:
#            if genomeVersion == 'hg19':
#                os.makedirs()
        print('Generate analysis folders.......\033[1;32mDone\033[0m')

#copy files
def copyFile(list):
    command=list
    for cmd in command:
        subprocess.Popen(cmd,shell=True)
    print('Copy files.......\033[1;32mDone\033[0m')
    time.sleep(1)

#
def returnPath(str1,list1):
    pathList=list1
    genomeVersion=str1
    if genomeVersion == 'hg19':
        return pathList[2]
    elif genomeVersion == 'hg38':
        return pathList[8]
    else:
        print('\033[1;31mError: \033[0mYou need to input correct genome version hg19/hg38')  
    
#concat files
def concat(str1,str2,list1,list2,list3,str3):
    path=str1
    genomeVersion=str2
    cmd_bgzip=list1
    cmd_tabix=list2
    cmd_bcf=list3
    bcfPath=str3
    file=os.listdir(path)
    filenum=len(file)
    prefix=set()
    #Compress files
    pidd=set()
    for filename in file:
        prefix.add(filename.split('.')[0])
        #It is better to use string as subprocess input
        cmd=' '.join(cmd_bgzip[0:2])+filename
        print('\033[1;34mNOTICE\033[0m'+': <'+cmd+'>')
        p1=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p1.poll() == None:
            pidd.add(p1.pid)        
    print('Compressing files..............\033[1;32mDone\033[0m')
    time.sleep(1)    

    #Indexing 
    gz=os.listdir(path)
    for bgz in gz:
        cmd=' '.join(cmd_tabix[0:3])+bgz
        print('\033[1;34mNOTICE\033[0m'+': <'+cmd+'>')
        p2=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p2.poll() == None:
            pidd.add(p2.pid)
    print('Indexing files..............\033[1;32mDone\033[0m')
    time.sleep(1)   
    
    #Concat files
    for id in prefix:
        cmd1=''.join(cmd_bcf[0:12])
        cmd=cmd1.format(id,id,bcfPath,id)
        print('\033[1;34mNOTICE\033[0m'+': <'+cmd+'>')
        p3=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p3.poll() ==None:
            pidd.add(p3.pid)
    print('Concate snp/indel files.......\033[1;32mDone\033[0m')
    time.sleep(1)
 
#LiftOver hg19Thg38
def liftover(str1,str2,str3,list1):
    pidd=set()
    genomeVersion=str1
    cmd_liftover=list1
    inpath=str2
    outpath=str3
    if genomeVersion == 'hg19':
        pass
    elif genomeVersion == 'hg38':
        file=os.listdir(inpath)
        for filename in file:
            cmd1=''.join(cmd_liftover)
            cmd=cmd1.format(filename,filename)
            p=subprocess.Popen(cmd,shell =True,stdin=subprocess.PIPE,close_fds=True)
            while p.poll()==None:
                pidd.add(p.pid)
        print('LiftOver hg19 To hg38.......\033[1;32mDone\033[0m')
    else:
        print('Please enter correct genome version hg19/hg38')
        sys.exit()
     
    if genomeVersion == 'hg38':
        p2=subprocess.Popen('rm '+outpath+'*unmap',shell =True,stdin=subprocess.PIPE,close_fds=True)
        while p2.poll()==None:
            pidd.add(p2.pid)
    else:
        pass
    time.sleep(1)

#Convert Format
def convert(str,list):
    pidd=set()
    path=str
    cmd_convert2annovar=list
    file=os.listdir(path)
    
    #Format convert
    for filename in file:
        cmd1=''.join(cmd_convert2annovar[0:11])
        cmd=cmd1.format(filename,filename.split('_')[0])
        print(cmd)
        p1=subprocess.Popen(cmd,shell =True,stdin=subprocess.PIPE,close_fds=True)
        while p1.poll() ==None:
            pidd.add(p1.pid)
    print('Convert to annovar.......\033[1;32mDone\033[0m')
    time.sleep(1)

#Annotate files---------------------------------------------------------------------------------
class people:
    
    def __init__(self,name,cmdList):
        self.name=name
        self.cmd=cmdList

    def annotate(self):
        prefix=self.name.split('.')[0]
        command=''.join(self.cmd).format(prefix,prefix)
        p=subprocess.Popen(command,shell=True)
        return p

def prepare_sample(genomeVersion,inpath,cmd_hg19_dbannotate,cmd_hg38_dbannotate):
    file=os.listdir(inpath)
    container=[]
    if genomeVersion == 'hg19':
        cmd=cmd_hg19_dbannotate
    else:
        cmd=cmd_hg38_dbannotate
    for filename in file:
        sample=people(filename,cmd)
        container.append(sample)
    return container

def dbannotate(container,numOfSubprocess):
    pidd=[]
    process=[]
    for i in range(0,len(container),numOfSubprocess):
        subContainer=container[i:i+numOfSubprocess]
        for sample in subContainer:
            p=sample.annotate()
            pidd.append(p.poll())
            process.append(p)
        while None in pidd:
            pidd=update(process)
        process=[]
    print('DataBases annotation.......\033[1;32mDone\033[0m')

#update the pid list
def update(process):
    update_pidd=[]
    for i in process:
        update_pidd.append(i.poll())
    return update_pidd
#------------------------------------------------------------------------------------------


#Reformat
def reformat(str,list):
    pidd=set()
    inpath=str
    cmd_reformat=list
    file=os.listdir(inpath)
    for filename in file:
        cmd1=' '.join(cmd_reformat)
        cmd=cmd1.format(filename)
        print(cmd)
        p1=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p1.poll() == None:
            pidd.add(p1.pid)
    time.sleep(1)

#
def transfer(str1,list1,list2):
    pidd=set()
    path=str1
    cmd_transfer=list1
    cmd_rm=list2
    file=os.listdir(path)
    for filename in file:
        cmd1=''.join(cmd_transfer)
        cmd2=''.join(cmd_rm)
        cmd=cmd1.format(filename,filename)
        p1=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p1.poll() == None:
            pidd.add(p1.pid)

        cmdd=cmd2.format(filename)
        p2=subprocess.Popen(cmdd,shell=True,stdin=subprocess.PIPE,close_fds=True)
        while p2.poll() == None:
            pidd.add(p2.pid)

#Check available clinvar edition
def checkClinAvail(str1,str2):
    genomeVersion=str1
    clinDate=str2
    template=genomeVersion+'_clinvarinfo_'+clinDate
    hg=[]
    file=os.listdir('/home/ubuntu/software/annovar/humandb/clinvarinfo_db')

    for item in file:
        if item.split('_')[0] == 'hg19' or 'hg38':
            hg.append([item.split('_')[0],'------',item.split('_')[2]])
        else:
            print('\033[1;31mError: \033[0m You should enter correct genom version: hg19 or hg38')
    
    if template in file:
        pass
    else:
        hgarraystr=str(np.array(hg))
        new=hgarraystr.replace('[','').replace(']','').replace("'",'').replace(' ','')
        print('\033[1;31mError: \033[0mCannot find the Clinvar database\n')
        print('Clinvar Available:')
        print('++++++++++++++++++\n')
        print(new+'\n')
        print('++++++++++++++++++\n')
        sys.exit()

#I/O 
def main():
    try:
        if sys.argv[1:]!=[]:
            opts,args=getopt.getopt(sys.argv[1:],'hg:i:o:d:p:')
            optslist=[]
            for i in opts:
                optslist.append(i[0])
            if len(args) != 0:
                print('\033[1;31mError: \033[0mYou need to specify input options and parameters. Check -h to see how to use this programme')
                sys.exit()
            elif len(opts) == 1:
                if opts[0][0] == '-h':
                    pass
                else:
                    print('\033[1;31mError: \033[0mSome options are missing. Check -h to see how to use this programme')
                    sys.exit()
            elif '-g' and '-n' and '-f' and '-d' not in optslist:
                print('\033[1;31mError: \033[0mSome options are missing. Check -h to see how to use this programme')
        else:
            print('\033[1;31mError: \033[0mNo input options')
            sys.exit()

    except getopt.GetoptError as err:
        print(err)
        print('Please check your options')
        print("Please type in: python annotate.py -h")
        sys.exit(2)
    
    values_dict={}
    for op,value in opts:
        if op == '-h':
            print("Usage: annotation.py [-h] [-g genome -o name -i NWZfolder -d date_clinvar] [-p number]")
            print("Options:")
            print('  -g    STR      The genome version hg19/hg38')
            print('  -o    STR      The name of the output folder')
            print('  -i    STR      The NWZ directory containing VCF files')
            print('  -d    STR      The date of clinvar (eg.20191231)')
            print('  -p    INT      The number of external subprocessess (default: 3)')
            print('NOTICE: -g genome -n name -f folder -d date  must be provided together')
            sys.exit()
        else:
            values_dict[op]=value
    return(values_dict)

#
if __name__=='__main__':
    args=main()
    #Prepare input
    genomeVersion=args['-g']
    dirName=args['-o']
    projectID=args['-i']
    clinDate=args['-d']
    if args['-p'].isdigit() == True:
            numOfSubprocess=int(args['-p'])
    else:
        numOfSubprocess=3
        
    dirPrefix=genomeVersion+'_'+dirName
#Thie project ID should be the name of NWZ platform project

#Prepare path
    myPath='/home/ubuntu/guoping/03-annotation/3.5-result/'
#Main folder
    pathForResults=myPath+dirPrefix+'_results/'
#sub foldfers
    pathForConcat=pathForResults+'concat/'
    pathForLiftOver=pathForResults+'liftover/'
    pathForConvert=pathForResults+'formatcovert/'
    pathForAnnotate=pathForResults+'dbannotate/'
    pathForReformat=pathForResults+'reformat/'
    pathForOutput=pathForResults+'results/'
    pathFortxt=pathForOutput+'txt/'
    pathForcsv=pathForOutput+'csv/'
#nwz folders
    nwzSnp='/mnt/Vazyme/cloud_platform/platform/vazyme/project/'+projectID+'/'+projectID+'/final_result/result/result_variation/snp/'
    nwzIndel='/mnt/Vazyme/cloud_platform/platform/vazyme/project/'+projectID+'/'+projectID+'/final_result/result/result_variation/indel/'
#pathList
    pathList=[pathForResults,pathForConcat,pathForConvert,pathForAnnotate,pathForReformat,pathForOutput,pathFortxt,pathForcsv,pathForLiftOver]
#Command
    cmd_copyfile=['cp '+nwzSnp+'*vcf* '+pathList[1],'cp '+nwzIndel+'*vcf* '+pathList[1]]
    cmd_bgzip=['bgzip',pathList[1]]
    cmd_tabix=['tabix','-C',pathList[1]]
    cmd_bcftools=['bcftools ','concat',' -a ',pathList[1],'{}','.snp.vcf.xls.gz ',pathList[1],'{}','.indel.vcf.xls.gz ','-o ','{}','{}','_all.vcf']
    cmd_liftover=['CrossMap.py ','vcf ','/home/ubuntu/software/annovar/humandb/','GRCh37_to_GRCh38.chain.gz ',pathList[8],'{} ',
                  '/home/ubuntu/software/annovar/humandb/','Homo_sapiens.GRCh38.dna.toplevel.fa ',pathList[2],'{}']
    cmd_convert2annovar=['convert2annovar.pl',' -format',' vcf4 ',pathList[2],'{}',' -outfile ',pathList[3],'{}','.avinput',' --withfreq',' -includeinfo']

    cmd_hg19_dbannotate=['table_annovar.pl ',pathList[3],'{}','.avinput ','~/software/annovar/humandb/ ','-buildver hg19 ','-out ',pathList[4],'{}',
                         ' -remove ','--protocol refGene,exac03,clinvar_',clinDate,',cytoBand,esp6500siv2_all,genomicSuperDups,gwasCatalog,snp142 ',
                         '-operation g,f,f,r,f,r,r,f ',
                         '-nastring . ','--polish ', '-xreffile ~/software/annovar/humandb/hg19_refGene.txt ','--otherinfo ','--thread 6']

    cmd_hg38_dbannotate=['table_annovar.pl ',pathList[3],'{}','.avinput',' ~/software/annovar/humandb/ ','-buildver hg38 ','-out ',pathList[4],'{}',
                         ' -remove ','--protocol refGene,ensGene,exac03,clinvar_',clinDate,',esp6500siv2_all,avsnp150,dbnsfp35c,1000g2015aug,gnomad211_exome ',
                         '-operation g,g,f,f,f,f,f,f,f ',
                         '-nastring . ','--polish ','-xreffile ~/software/annovar/humandb/hg38_refGene.txt ','--otherinfo ','--thread 6']

    cmd_reformat=['reformat.py',pathList[4],pathList[5],genomeVersion,clinDate,' {}']
    cmd_transfer=['cat ',pathList[6],'{}',' | ',"sed 's# ##g' | sed 's#1/2#1/1#g' > ",pathList[6],'{}_snp_indel.txt']
    cmd_rm=['rm ',pathList[6],'{}']

#Start time
    start = time.asctime(time.localtime(time.time()))
#Function 
    checkClinAvail(genomeVersion,clinDate)
#Creating folders
#    mkdir(0,pathList,genomeVersion)
#Copy files
#    copyFile(cmd_copyfile)
#Return path for cmd_bcf
#    bcfPath=returnPath(genomeVersion,pathList)
#Concat files
#    concat(pathList[1],genomeVersion,cmd_bgzip,cmd_tabix,cmd_bcftools,bcfPath)
#LiftOver hg19Thg38
#    liftover(genomeVersion,pathList[8],pathList[2],cmd_liftover)
#Convert files
#    convert(pathList[2],cmd_convert2annovar)
#Annotate file
#    dbannotate(prepare_sample(genomeVersion,pathList[3],cmd_hg19_dbannotate,cmd_hg38_dbannotate),numOfSubprocess)
#Reformat
    reformat(pathList[4],cmd_reformat)
#transfer
    transfer(pathList[6],cmd_transfer,cmd_rm)
#End time
    end = time.asctime(time.localtime(time.time()))
#Finish
    print('\033[1;32mCongratulation! ALL FINISHED!\033[0m')
    print('START Time:',start)
    print('END Time: ',end)
