#!/bin/bash 

#SBATCH --partition=amilan
#SBATCH --job-name=micc_pipeline
#SBATCH --output=test_micro-c_pipeline.%j.out
#SBATCH --time=168:00:00
#SBATCH --qos=long
#SBATCH --nodes=18
#SBATCH --ntasks=18

module purge
module load anaconda
conda activate MicC # you will need an active conda environment with the bowtie2, samtools, pairtools, and cooler packages 

fastq1=$1 
fastq2=$2

gunzip $fastq1 #if your fastq files are already unzipped feel free to # these out
gunzip $fastq2

file_path_to_genome="/pl/active/swygertlab/jasonher/Saccer3/Saccer3" #you need to change this to the location of your genome data
file_path_to_chrom_sizes="pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes" # you need to change this to the location of your chromomsomes' sizes file
distance_graphed=201 #change to distance you want to get graphed
sample3=$3 # include how you want the files to be named 

sample1=$(echo ${fastq1} | sed 's/\..*$//')
sample2=$(echo ${fastq2} | sed 's/\..*$//')

bzip2 -d ${fastq1}
bzip2 -d ${fastq2}

#bowtie2 will align the sample reads to the reference genome
bowtie2 --very-sensitive -p 16 --reorder -x $file_path_to_genome -1 ${sample1}.fastq -2 ${sample2}.fastq -S ${sample3}.sam

#turning the bowtie2 .sam file output into a .bam file and using the .bam file to make a .pairs file
samtools view -S -b ${sample3}.sam > ${sample3}.bam
samtools view -h ${sample3}.bam | pairtools parse -c $file_path_to_chrom_sizes -o ${sample3}_parsed.pairs.gz

#pairtools sort puts the reads in base pair sequential order
#dedup removes duplicates
#select removes based on some criteria here we only want unrescued reads for more info check out the pairtools select documentation
pairtools sort --nproc 8 --tmpdir=./ -o ${sample3}_sorted.pairs.gz ${sample3}_parsed.pairs.gz
pairtools dedup --mark-dups -o ${sample3}_deduped.pairs.gz ${sample3}_sorted.pairs.gz
pairtools select '(pair_type == "UU")' -o ${sample3}_filtered.pairs.gz ${sample3}_deduped.pairs.gz 
pairtools split --output-pairs ${sample3}_output.pairs.gz ${sample3}_filtered.pairs.gz

#in order to access files in python script unzipping them
#the filter_orientation_heading.py python script will generate IN, OUT, SAME, and NoFilter .pairs files the python scripts must be in working directory with command as written
gunzip ${sample3}_output.pairs.gz
python filter_orientations_heading.py ${sample3}_output.pairs

java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample3}_output_IN_reads.pairs ${sample3}'_IN_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample3}_output_OUT_reads.pairs ${sample3}'_OUT_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample3}_output_SAME_reads.pairs ${sample3}'_SAME_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample3}_output_noIN.pairs ${sample3}'_noIN.hic' sacCer3

#rezipping files
in=$(wc -l ${sample3}_output_IN_reads.pairs)
out=$(wc -l ${sample3}_output_OUT_reads.pairs)
same=$(wc -l ${sample3}_output_SAME_reads.pairs)
noIN=$(wc -l ${sample3}_output_noIN.pairs)
in_reads=$(echo ${in} | cut -d ' ' -f 1)
out_reads=$(echo ${out} | cut -d ' ' -f 1)
same_reads=$(echo ${same} | cut -d ' ' -f 1)
noIN_reads=$(echo ${noIN} | cut -d ' ' -f 1)
sum=$((${in_reads}+${out_reads}+${same_reads}))
#python distance_decay.py script generates short distance decay plots from 0 to 2000 bp
python distance_decay.py ${sample3}_output_IN_reads.pairs ${in_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample3}_output_OUT_reads.pairs ${out_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample3}_output_SAME_reads.pairs ${same_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample3}_output_noIN.pairs ${noIN_reads} ${sum} $distance_graphed "False"

#rezipping files
bgzip ${sample3}_output_IN_reads.pairs
bgzip ${sample3}_output_OUT_reads.pairs
bgzip ${sample3}_output_SAME_reads.pairs
bgzip ${sample3}_output_noIN.pairs

#creating cooler files
pairix -f ${sample3}_output_IN_reads.pairs.gz
pairix -f ${sample3}_output_OUT_reads.pairs.gz
pairix -f ${sample3}_output_SAME_reads.pairs.gz
pairix -f ${sample3}_output_noIN.pairs.gz

cooler cload pairix $file_path_to_chrom_sizes:150 ${sample3}_output_IN_reads.pairs.gz ${sample3}_output_IN_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample3}_output_OUT_reads.pairs.gz ${sample3}_output_OUT_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample3}_output_SAME_reads.pairs.gz ${sample3}_output_SAME_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample3}_output_noIN.pairs.gz ${sample3}_output_noIN.cool

cooler balance ${sample3}_output_IN_reads.cool
cooler zoomify ${sample3}_output_IN_reads.cool
cooler balance ${sample3}_output_OUT_reads.cool
cooler zoomify ${sample3}_output_OUT_reads.cool
cooler balance ${sample3}_output_SAME_reads.cool
cooler zoomify ${sample3}_output_SAME_reads.cool
cooler balance ${sample3}_output_noIN.cool
cooler zoomify ${sample3}_output_noIN.cool

gzip $sample1.fastq
gzip $sample2.fastq
gzip $sample3.sam
