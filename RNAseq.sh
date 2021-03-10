#!/bin/bash
#SBATCH --job-name=RNA    ## Name of job
#SBATCH -A ecoevo283       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-20       ## number of tasts to launch (prefixes has 24 lines so 24)
#SBATCH --cpus-per-task=2  ## number of cores needed 

module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load java/1.8.0
module load picard-tools/1.87

# or pass the file name to the shell script, how would I do this?
file="shortRNAseq.prefix.txt"
# is the file indexed for bwa?
ref="/data/class/ecoevo283/lbradle1/AdvInf_Winter2021/ref/dmel-all-chromosome-r6.13.fasta"
dir="/data/class/ecoevo283/lbradle1/AdvInf_Winter2021/RNAseq"
# here is a hint if you had a tab delimited input file
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1`
samplename=`echo $prefix | sed -e "s/rawdata/bam/"`


## loop through files as a SLURM job (you have to write).... as an array job
## see how I define an output directory
hisat2 -p 2 -x $ref -1 $prefix_1.fq.gz -2 $prefix_2.fq.gz -S $dir/$samplename.sam
samtools view -bS $dir/$samplename.sam > $dir/$samplename.bam
samtools sort $dir/$samplename.bam $dir/$samplename.sorted
samtools index $dir/$samplename.sorted.bam