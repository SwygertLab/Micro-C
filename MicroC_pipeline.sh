#!/bin/bash 

#SBATCH --partition=amem
#SBATCH --job-name=micro-c
#SBATCH --output=micro-c.%j.out
#SBATCH --time=48:00:00
#SBATCH --qos=mem
#SBATCH --nodes=2
#SBATCH --ntasks=36

module purge
module load anaconda
conda activate MicC # you will need an active conda environment with the bowtie2, samtools, pairtools, and cooler packages 

mkdir ./data
mkdir ./results
mkdir ./log

fastq1=$1 
fastq2=$2
file_path_to_genome="/pl/active/swygertlab/jasonher/Saccer3/Saccer3" #you need to change this to the location of your genome data
file_path_to_chrom_sizes="/pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes" # you need to change this to the location of your chromomsomes' sizes file

distance_graphed=201 #change to distance you want to get graphed
sample_name=$3 # include how you want the files to be named 
fastq1_name=$(echo ${fastq1} | sed 's/\..*$//')
fastq2_name=$(echo ${fastq2} | sed 's/\..*$//')

gunzip $fastq1 #if your fastq files are already unzipped feel free to # these out
gunzip $fastq2



#bowtie2 will align the sample reads to the reference genome
bowtie2 --very-sensitive -k 1 -p $SLURM_NTASKS --reorder -x $file_path_to_genome -1 ${fastq1_name}.fastq -2 ${fastq2_name}.fastq -S ${sample_name}.sam

gzip ${sample1}.fastq
gzip ${sample2}.fastq
mv ${fastq1} ./data/
mv ${fastq2} ./data/

turning the bowtie2 .sam file output into a .bam file and using the .bam file to make a .pairs file
samtools view -S -b ${sample_name}.sam > ${sample_name}.bam
samtools view -h ${sample_name}.bam | pairtools parse -c $file_path_to_chrom_sizes -o ${sample_name}_parsed.pairs.gz

#pairtools sort puts the reads in base pair sequential order
#dedup removes duplicates
#select removes based on some criteria here we only want unrescued reads for more info check out the pairtools select documentation
pairtools sort --nproc 8 --tmpdir=./ -o ${sample_name}_sorted.pairs.gz ${sample_name}_parsed.pairs.gz
pairtools dedup --mark-dups -o ${sample_name}_deduped.pairs.gz ${sample_name}_sorted.pairs.gz
pairtools select '(pair_type == "UU")' -o ${sample_name}_filtered.pairs.gz ${sample_name}_deduped.pairs.gz 
pairtools split --output-pairs ${sample_name}_output.pairs.gz ${sample_name}_filtered.pairs.gz

#in order to access files in python script unzipping them
#the filter_orientation_heading.py python script will generate IN, OUT, SAME, and NoFilter .pairs files the python scripts must be in working directory with command as written
gunzip ${sample_name}_output.pairs.gz
python filter_orientations_heading.py ${sample_name}_output.pairs

java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample_name}_output_IN_reads.pairs ./results/${sample_name}'_IN_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample_name}_output_OUT_reads.pairs ./results/${sample_name}'_OUT_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample_name}_output_SAME_reads.pairs ./results/${sample_name}'_SAME_reads.hic' sacCer3
java -Xmx22g -jar juicer_tools_1.22.01.jar -r 10,50,100,150,200,500,1000,3200,5000 pre ${sample_name}_output_noIN.pairs ./results/${sample_name}'_noIN.hic' sacCer3

in=$(wc -l ${sample_name}_output_IN_reads.pairs)
out=$(wc -l ${sample_name}_output_OUT_reads.pairs)
same=$(wc -l ${sample_name}_output_SAME_reads.pairs)
noIN=$(wc -l ${sample_name}_output_noIN.pairs)
in_reads=$(echo ${in} | cut -d ' ' -f 1)
out_reads=$(echo ${out} | cut -d ' ' -f 1)
same_reads=$(echo ${same} | cut -d ' ' -f 1)
noIN_reads=$(echo ${noIN} | cut -d ' ' -f 1)
sum=$((${in_reads}+${out_reads}+${same_reads}))
#python distance_decay.py script generates short distance decay plots from 0 to 2000 bp
python distance_decay.py ${sample_name}_output_IN_reads.pairs ${in_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample_name}_output_OUT_reads.pairs ${out_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample_name}_output_SAME_reads.pairs ${same_reads} ${sum} $distance_graphed "False"
python distance_decay.py ${sample_name}_output_noIN.pairs ${noIN_reads} ${sum} $distance_graphed "False"

#creating cooler files
pairix -f ${sample_name}_output_IN_reads.pairs.gz
pairix -f ${sample_name}_output_OUT_reads.pairs.gz
pairix -f ${sample_name}_output_SAME_reads.pairs.gz
pairix -f ${sample_name}_output_noIN.pairs.gz

cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_IN_reads.pairs.gz ${sample_name}_output_IN_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_OUT_reads.pairs.gz ${sample_name}_output_OUT_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_SAME_reads.pairs.gz ${sample_name}_output_SAME_reads.cool
cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_noIN.pairs.gz ${sample_name}_output_noIN.cool

cooler balance ${sample_name}_output_IN_reads.cool
cooler zoomify ${sample_name}_output_IN_reads.cool
cooler balance ${sample_name}_output_OUT_reads.cool
cooler zoomify ${sample_name}_output_OUT_reads.cool
cooler balance ${sample_name}_output_SAME_reads.cool
cooler zoomify ${sample_name}_output_SAME_reads.cool
cooler balance ${sample_name}_output_noIN.cool
cooler zoomify ${sample_name}_output_noIN.cool

#rezipping files
bgzip ${sample_name}_output_IN_reads.pairs
bgzip ${sample_name}_output_OUT_reads.pairs
bgzip ${sample_name}_output_SAME_reads.pairs
bgzip ${sample_name}_output_noIN.pairs
gzip ${sample_name}.sam

#removing unneeded files
rm ${sample_name}_parsed.pairs.gz
rm ${sample_name}_sorted.pairs.gz
rm ${sample_name}_deduped.pairs.gz
rm ${sample_name}_filtered.pairs.gz
