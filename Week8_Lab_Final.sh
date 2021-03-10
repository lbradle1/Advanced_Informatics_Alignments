#####################################
### NEED TO INDEX REF GENOME FIRST ##
#####################################

#Run the following in my AdvInf_Winter2021 (above ref) dir

# grab a non-head node
srun -A ecoevo283 --pty bash -i

# check that you have both the dmell fasta and gtf files in a ref driretory 

# loaded Modules that Tony said I'd need 
module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load java/1.8.0
module load hisat2/2.2.1

# this makes my .fasta file in ref directory be called 'ref'
# makes more files that apparently i need to make a SAM file
# https://informatics.fas.harvard.edu/short-introduction-to-bwa.html 
# makes a SAM file that we need to align stuff 

#index in each language that you'll be using so bwa/sam/java/hisat
ref="ref/dmel-all-chromosome-r6.13.fasta"
bwa index $ref
samtools faidx $ref
# now we need to do something like making a dictionary or something idk 
# need to replace .fasta.dict (auto-output) wuth just .dict to make the programs play nice
java -jar /opt/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=$ref O=ref/dmel-all-chromosome-r6.13.dict
#makes alignment program to ref possible 
hisat2-build $ref $ref # $ref $ref = input and output 

## THIS WILL MAKE LOTS OF NEW FILES IN REF YAY

#if trying to find picard-tools 
module show picard-tools

#make the relevant directories for bam files where alignments will be output 
mkdir DNAseq/bam
mkdir RNAseq/bam
mkdir ATACseq/bam

######################################################################

##################################
### gDNA to Start ################
##################################

#running this also in AdvInf_Winter2021 dir - not sure that's right
## prefixes.txt needs to be in same DIR as myDNAjob.sh and the DIR where you RUN that with sbatch 

## now do some command line stuff to make a file of "prefixes" you can step through to do the alignment

# this line lists all files ending in _1.fq.gz in DNA/rawdata and filter/edits them (sed)
# 's/regexp/replacement/' - replaces _1.fq.gz with nothing i.e. leaving only the prefixes in a new text file   
ls DNAseq/rawdata/*_1.fq.gz | sed 's/_1.fq.gz//' >prefixes.txt


####################################
##### DNA Alignments ###############
##### Submitted myDNAjob.sh ########
####################################

# add myDNAjob.sh to DIR containing prefixes.txt (for me this was DNAseq) 

# make sure file is in unix format 
dos2unix myDNAjob.sh

#submit job - job ID 3566845
sbatch myDNAjob.sh

#check that it's running 
squeue -u lbradle1
# 3566845 = job ID


####################
## myDNAjob.sh #####


#!/bin/bash
#SBATCH --job-name=DNA     ## Name of job
#SBATCH -A ecoevo283       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-12       ## number of tasts to launch (prefixes has 12 lines so 12)
#SBATCH --cpus-per-task=2  ## number of cores needed 

module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load java/1.8.0
module load picard-tools/1.87

# or pass the file name to the shell script, how would I do this?
file="prefixes.txt"
# is the file indexed for bwa?
ref="ref/dmel-all-chromosome-r6.13.fasta"

# here is a hint if you had a tab delimited input file
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1`
# whats this do - takes first n lines of file $SLURM ARRAY (lines 1-12) pipes that into tail -n 1
# prefix = nth line - allows us to process things in parallel 
### IMPORTANT CAN RUN THINGS IN PARALLEL ###

#changes the path - writes it out to rawdata/bam
samplename=`echo $prefix | sed -e "s/rawdata/bam/"`

#pipe prefix to cut, take 3rd field, pass to 2nd cut command at _ gives sample id's without technical replicate number 
idname=`echo $prefix | cut -d "/" -f 3 | cut -d "_" -f 1`


# alignments
# -t 2 is # or cores

bwa mem -t 2 -M $ref ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -bS - > $samplename.bam
samtools sort $samplename.bam -o $samplename.sort.bam
# GATK likes readgroups
java -jar  /opt/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=$samplename.sort.bam O=$samplename.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$idname RGSM=$idname VALIDATION_STRINGENCY=LENIENT
samtools index $samplename.RG.bam

### end myDNAjob.sh ###
#######################


#############################################################################


#######################
## ATAC seq ###########
#######################

#ran this inside the ATACsec/ DIR - thus .txt is in ATACseq DIR
ls rawdata/*_F.fq.gz | sed 's/_F.fq.gz//' >ATACprefixes.txt

# add ATACjob.sh to ATACseq DIR

# make sure file is in unix format 
dos2unix ATACseq.sh

#submit job - job ID 3566845
sbatch ATACseq.sh

#check that it's running 
squeue -u lbradle1
# job ID 3569526 

###################################################################################

####################
## RNA Stuff #######
####################

# MAKE A SUBSET OF RNAseq DATA 
###  I will focus on a small subset of the dataset for speed

R
mytab = read.table("/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384_SampleCoding.txt",header=TRUE)
# drop some columns
mytab2 = mytab[,9:12]
# select a subset of samples
mytab3 = mytab2[c(1:10,93:102),]
# rearrange and drop another column
mytab4 = mytab3[,c(4,1,2)]
if (file.exists("shortRNAseq.names.txt")) {file.remove("shortRNAseq.names.txt")}
for(i in 1:nrow(mytab4)){
	cat("RNAseq/bam/",mytab4$FullSampleName[i],".bam\n",file="shortRNAseq.names.txt",append=TRUE,sep='')
	}
write.table(mytab4,"shortRNAseq.txt")


write.table(mytab4[,1],"shortRNAseq.prefix.txt",row.names=FALSE, col.names=FALSE)
# currently getting error about quotes and colnames for the line above 


### DO STUFF WITH THAT SUBSET 
file="shortRNAseq.prefix.txt"
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1`

####################################
##### RNA Alignments ###############
##### Submitted RNAjob.sh ########
###################################

srun -A ecoevo283 --pty bash

sbatch RNAseq.sh
# job ID 3637406

#check that something is running
squeue -u lbradle1

