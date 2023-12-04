#!/bin/bash
# source /app/Lmod/lmod/lmod/init/bash
#Pavitra Roychoudhury, Sep 2017

#Load required tools before running this script
# module load bowtie2/2.2.5
# module load SAMtools/1.8-foss-2016b
# module load FastQC/0.11.8-Java-1.8
# module load R/3.6.2-foss-2016b-fh1
# module load BBMap/38.44-foss-2016b

PATH=$PATH:$HOME/.local/bin:

echo $SLURM_CPUS_PER_TASK

while getopts ":1:2:s:r:g:" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		s) in_fastq="$OPTARG"
		paired="false"
		;;
		r) ref="$OPTARG"
		;;
        g) tgt="$OPTARG"
        ;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

if [[ -z "$ref" ]] 
then 
echo "Reference missing"
exit 1
fi

##  PAIRED-END  ##
if [[ $paired == "true" ]]
then
if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq_r1%%_R1.fastq*})

#FastQC report on raw reads
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2  out1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz' ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'  out1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' './trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' out1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' out2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20


#FastQC report on processed reads
mkdir -p ./fastqc_reports_trimmed
fastqc './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_unpaired_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_unpaired_r2.fastq.gz' -o ./fastqc_reports_trimmed

#Map reads to reference
mkdir -p ./mapped_reads
bowtie2 -x $ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'.sam'


##  SINGLE-END  ##
else 
if [[ $paired == "false" ]]
then
if [ -z $in_fastq ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq%%.fastq*})

#FastQC report on raw reads
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq -o ./fastqc_reports_raw

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in=$in_fastq out='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'  out='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'

#Quality trimming
mkdir -p ./preprocessed_fastq
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' out='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#FastQC report on processed reads
mkdir -p ./fastqc_reports_trimmed
fastqc './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o ./fastqc_reports_trimmed

#Map reads to reference
mkdir -p ./mapped_reads
bowtie2 -x $ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'.sam'
 
fi
fi

#Make a sorted bam file
samtools view -bh -q 10 -o './mapped_reads/'$sampname'.bam' -T $ref'.fasta' './mapped_reads/'$sampname'.sam'  
rm './mapped_reads/'$sampname'.sam'
samtools sort -o './mapped_reads/'$sampname'.sorted.bam' './mapped_reads/'$sampname'.bam' 
rm './mapped_reads/'$sampname'.bam'

reffasta=$ref'.fasta';


#Call R script to count variant reads
mkdir -p ./results
if [[ $paired == "true" ]]
then
# Rscript --vanilla hsv_variant_analysis.R tgt_region=\"$tgt\" paired=TRUE s1=\"$(basename $in_fastq_r1)\" s2=\"$(basename $in_fastq_r2)\" refseq_fname=\"$reffasta\"
Rscript --vanilla hsv_variant_analysis.R tgt_region=\"$tgt\" paired=TRUE s1=\"$in_fastq_r1\" s2=\"$in_fastq_r2\" refseq_fname=\"$reffasta\"
else
if [[ $paired == "false" ]]
then
Rscript --vanilla hsv_variant_analysis.R tgt_region=\"$tgt\" paired=FALSE s1=\"$(basename $in_fastq)\" refseq_fname=\"$reffasta\"
fi
fi



