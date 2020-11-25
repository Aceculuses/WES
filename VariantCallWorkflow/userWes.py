#!/data2/guoping/software/miniconda3/envs/gatk/bin/python3
import subprocess
import time
import sys
def mapping():
#deal with input
    neg=['N','n','NO','No','no','nO']
    pos=['Y','y','Yes','YEs','YES','yes','yES','yES']
    pathlist=[]
    samplelist=[]
    outpathlist=[]

    path = input('Please enter the input directory path: ')
    if len(path) != 0:
        pathlist.append(path)
    else:
        while len(path) == 0:
            path = input('No input! Please Re-enter the directory path: ')
        pathlist.append(path)

    output=input('Please enter the output directory path: ')
    if len(output) != 0:
        outpathlist.append(output)
    else:
        while len(output)==0:
            output = input('No input! Please Re-enter the output directory path: ')
            outpathlist.append(output)

    samples = input('Please enter the samples name: ')
    samplelist.append(samples)
    while len(samples) == 0:
        samples = input('No input! Please Re-enter the samples name: ')
        samplelist.append(samples)

    condition=input('Are there all the samples to be analyzed?[Y/N]: ')
    while len(condition) == 0:
        condition=input('No input! Please Re-enter [Y/N]: ')

    while condition in neg:
        samples = input('Please enter the next sample: ')
        samplelist.append(samples)
#    print(samplelist)
        while len(samples) == 0:
            samples = input('No input! Please Re-enter the samples name: ')
            samplelist.append(samples)
        condition=input('Are there all the samples to be analyzed?[Y/N]: ')
        while len(condition) == 0:
           condition=input('No input! Please Re-enter [Y/N]: ')

#create wesinfo.txt and directory and run variantCall_mapping.sh
    calling=[]
    for i in samplelist:
        if i != '':
            file=open(outpathlist[0]+i+'_wesinfo.txt','w')
            file.write(pathlist[0]+' '+i+' '+outpathlist[0]+i+'/'+'\n')
            info=pathlist[0]+' '+i+' '+outpathlist[0]+i+'/'
            calling.append(info)
            cm=['mkdir -p ',outpathlist[0],i+'/']
            cmd=''.join(cm)
            subprocess.Popen(cmd,shell=True)
    time.sleep(0.5)
    print('Start running variantCall_mapping.sh. It will take a moment to complete.')
    for x in calling:
        name=x.split('/')[-2]
        cm='nohup variantCall_mapping.sh'+' '+x+' '+'>'+' '+x.split(' ')[2]+name+'_mapping.log'+' '+'2>&1 &'
        subprocess.Popen(cm,shell=True)
        print(cm)
        time.sleep(0.5)

def calling():
#deal with input
    neg=['N','n','NO','No','no','nO']
    pos=['Y','y','Yes','YEs','YES','yes','yES','yES']
    pathlist=[]
    samplelist=[]
    bedlist=[]
    default='20200602'
    glist=['hg19','hg38']
    genomelist=[]
    clinvarlist=[]

    path = input('Please enter the wesinfo.txt directory path: ')
    if len(path) != 0:
        pathlist.append(path)
    else:
        while len(path) == 0:
            path = input('No input! Please Re-enter the directory path: ')
        pathlist.append(path)   

    genome = input('Please enter the genome version hg19/hg38 (default: hg19): ')
    if len(genome) != 0:
        while genome not in glist:
            genome=input('Please enter correct genome version hg19 or hg38: ')
        genomelist.append(genome)
    else:
        while len(genome) == 0:
            genome = input('No input! Please Re-enter the genome hg19 or hg38: ')
        genomelist.append(genome)
        while genome not in glist:
            genome=input('Please enter correct genome version hg19 or hg38: ')
        genomelist.append(genome)
    version=genomelist[-1]
                
    bed = input('Please enter the interval file (bed/interval): ')
    if len(bed) != 0:
        bedlist.append(bed)
    else:
        while len(bed) == 0:
            bed = input('No input! Please Re-enter the interval file (bed/interval): ')
        bedlist.append(bed) 
    
    clinvar = input('Please enter the clinvar date (default : 20200602): ')
    if len(clinvar) != 0:
        clinvarlist.append(clinvar)
    else:
        clinvarlist.append(default)
        print('Using the deafult clinvar database')

    samples = input('Please enter the samples name: ')
    samplelist.append(samples)
    while len(samples) == 0:
        samples = input('No input! Please Re-enter the samples name: ')
        samplelist.append(samples)

    condition=input('Are there all the samples to be analyzed?[Y/N]: ')
    while len(condition) == 0:
        condition=input('No input! Please Re-enter [Y/N]: ')

    while condition in neg:
        samples = input('Please enter the next sample: ')
        samplelist.append(samples)
        while len(samples) == 0:
            samples = input('No input! Please Re-enter the samples name: ')
            samplelist.append(samples)
        condition=input('Are there all the samples to be analyzed?[Y/N]: ')
        while len(condition) == 0:
           condition=input('No input! Please Re-enter [Y/N]: ')

    for i in samplelist:
        if i != '':
            file=open(pathlist[0]+i+'_wesinfo.txt','r').readlines()
            for x in file:
                outpath=x.rstrip('\n').split(' ')[-1]
                infotable=open(outpath+i+'_infotable.txt','w')
                infotable.write('inputPath:'+outpath+'\n'+'outputPath:'+outpath+'\n'+'genomeVersion:'+version+'\n'+'clinvar_date:'+clinvarlist[-1]+'\n'+'file:'+i+'.hg19_multianno.txt')
                cmd=['nohup variantCall_combined.sh ',outpath,' ',i,' ',i+'_infotable.txt ',bedlist[-1],' > ',outpath,i,'_calling.log 2>&1 &']
                command=''.join(cmd)
                subprocess.Popen(command,shell=True)
                print(command)
                time.sleep(0.5)
#                print(outpath)
#    print(samplelist)
#    print(pathlist)
#    print(bedlist)
#    print(version)

if __name__=='__main__':
    tag=sys.argv[1]
    if tag == '1':
        mapping()
    else:
        calling()

