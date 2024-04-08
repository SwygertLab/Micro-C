#!/bin/bash 

#SBATCH --partition=amem
#SBATCH --job-name=micro-c
#SBATCH --output=micro-c.%j.out
#SBATCH --time=24:00:00
#SBATCH --qos=mem
#SBATCH --nodes=2
#SBATCH --ntasks=36

################################################
# PROGRAM:
# MicroC.sh
#
# DESCRIPTION:
#   This is a pipeline for taking micro-c fastq data and
#   converting them to hic and cooler files for contact 
#   heatmap visualization. It also generates some distance
#   decay plots and logs the processes. 
#
# AUTHOR:
# Jason Hernandez
#
# START DATE:
# May 18, 2021
#
# DEPENDENCIES:
# 	Requires the installation of the follwing software: 
#       bowtie2
#       pairtools version 0.3.0+
#       cooltools
#       numpy
#       matplotlib
#       python version 3.10+
#
# REQUIRES:
#    INPUT: .fastq files:    For each sample, paired forward and reverse sequencing files
#								are required. These should be placed in an input
#								directory.
#
#    INPUT: metadata.txt file: A metadata file with two columns. The first two columns
#								are fastq file names. The third column is a "nickname"
#								of each sample. Later columns can be included with other
#								metadata information. Metadata file should be placed
#								within the inputdir directory.
#
#
#
#    GENOME SEQUENCE: .fa  or .tar.gz file for the genome. This is the sequence of the 
#                                genome.
#
# USAGE:
# $ sbatch MicroC.sh <metadata.txt> <number of threads>
#
# OUTPUT: .sam file:
#              A file indicating the location of reads in reference to the reference genome.
#
# OUTPUT: _orientation_output.pairs files:
#              A series of files witht the different read orientations and their associated
#              read information in the .pairs file format. 
#
# OUTPUT: _distance_decay_orientation.png:
#              A series of graphical images containing the binned in 10 base pairs bins and their 
#              relationships between relative to that orientations total amount of reads over distance. 
#
# OUTPUT: contacts.txt:
#              A file containing the exact location pair of contacts.
#
# OUTPUT: microc_heatmap.hic / microc_heatmap.cool:
#              Micro-C heatmaps in both .hic format and .cool format are provided. .hic files are 
#              opened with JuiceBox, and .cool files are visualized using Hi-C glass. 
#
# KNOWN BUGS:
#
# THINGS TO IMPROVE:
#
################################################
module purge
module load anaconda
####### MODIFY THIS SECTION #############
conda activate MicC # you will need an active conda environment with the bowtie2, samtools, pairtools, and cooler packages 

#Metadata file. This pulls the metadata path and file from the command line
metadata=$1

#This is where the genome sequence lives:
file_path_to_genome="/pl/active/swygertlab/jasonher/Saccer3/Saccer3" #you need to change this to the location of your genome data

#This is where the chromosome sizes live:
file_path_to_chrom_sizes="/pl/active/swygertlab/jasonher/micro-c/sacCer3.chrSizes" # you need to change this to the location of your chromomsomes' sizes file

file_path_to_juicer_jar="/pl/active/swygertlab/jasonher/juicer_jar/juicer_tools_1.22.01.jar"

DATE=$(date +%Y-%m-%d)
#OR
#DATE='2022-12-03'
inputdir="/pl/active/swygertlab/jasonher/Micro-C/SCeres_logWT/01_input"
scriptsdir="/pl/active/swygertlab/jasonher/Micro-C/SCeres_logWT/02_scripts"
outputdir="/pl/active/swygertlab/jasonher/Micro-C/SCeres_logWT/03_results"
#note that these do not have the foward slash "/" on purpose

distance_graphed=201 #change to distance you want to get graphed by distance_decay.py

mapq_filter=2 #the filter by which you don't want reads to be included if they are lower quality, higher mapq, than

########## DONE MODIFYING ###############

########## BEGIN CODE ###############
echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"

samples=$(wc -l $metadata | awk '{ print $1 }')

fastq1=$(awk 'NR==1 { print $1 }' $metadata)
fastq2=$(awk 'NR==1 { print $2 }' $metadata)
sample_name=$(awk 'NR==1 { print $3 }' $metadata | sed 's/\..*$//') #this will the prefix for all the new files made

mkdir $outputdir/$DATE"_output"
mkdir $outputdir/pairs_files
mkdir $outputdir/cool_files
mkdir $outputdir/hic_files
mkdir $outputdir/distance_decay_plots

#these lines are only necessary if you did not make your own input and scripts directory
#mkdir $inputdir
#mkdir $scriptsdir
#mv fastq1 $inputdir
#mv fastq2 $inputdir
#mv separate_by_orientation.py $scriptsdir
#mv distance_decay.py $scriptsdir
#cd $inputdir

#bowtie2 will align the sample reads to the reference genome
cmd1="bowtie2 --very-sensitive -p $pthread --reorder -x $file_path_to_genome -1 ${fastq1} -2 ${fastq2} -S ${sample_name}.sam 2> ${sample_name}.bowtie2.log"
echo "Completing bowtie2 alignment with reference genome, here is the command used: "
echo $cmd1
time eval $cmd1

grep '@' ${sample_name}.sam > ${sample_name}'_'${mapq_filter}_filtered.sam
grep -v '@' ${sample_name}.sam | awk '{ if($5 >= '$mapq_filter') print $0;}' >> ${sample_name}'_'${mapq_filter}_filtered.sam

#turning the bowtie2 .sam file output into a .bam file and using the .bam file to make a .pairs file
samtools view -S -b ${sample_name}'_'${mapq_filter}_filtered.sam > ${sample_name}'_'${mapq_filter}_filtered.bam
samtools view -h ${sample_name}'_'${mapq_filter}_filtered.bam | pairtools parse -c $file_path_to_chrom_sizes -o ${sample_name}_parsed.pairs.gz

#pairtools sort puts the reads in base pair sequential order
#dedup removes duplicates
#select removes based on some criteria here we only want unrescued reads for more info check out the pairtools select documentation
echo "Beginning conversion of sam file to pairs file and cleaning up data"
cmd2="pairtools sort --nproc 8 --tmpdir=./ -o ${sample_name}_sorted.pairs.gz ${sample_name}_parsed.pairs.gz"
cmd3="pairtools dedup --mark-dups -o ${sample_name}_deduped.pairs.gz ${sample_name}_sorted.pairs.gz"
cmd4="pairtools select '(pair_type == \"UU\")' -o ${sample_name}_filtered.pairs.gz ${sample_name}_deduped.pairs.gz"
cmd5="pairtools split --output-pairs $outputdir/pairs_files/${sample_name}_output.pairs.gz ${sample_name}_filtered.pairs.gz"
echo $cmd2
time eval $cmd2
echo $cmd3
time eval $cmd3
echo $cmd4
time eval $cmd4
echo $cmd5
time eval $cmd5

#in order to access files in python script unzipping them
#the filter_orientation_heading.py python script will generate IN, OUT, SAME, and NoFilter .pairs files the python scripts must be in working directory with command as written
echo "Separating the one _output.pairs file into 4 different pairs files separated by orientation"
cd $outputdir/pairs_files/
gunzip ${sample_name}_output.pairs.gz
cmd6="python $scriptsdir/separate_by_orientation.py $outputdir/pairs_files/${sample_name}_output.pairs"
echo $cmd6
time eval $cmd6

echo "Creating hic files containing Micro-C heatmaps for each orientation"
cmd7="java -Xmx22g -jar $file_path_to_juicer_jar -r 10,50,100,150,200,500,1000,3200,5000 pre $outputdir/pairs_files/${sample_name}_output_IN_reads.pairs $outputdir/hic_files/${sample_name}'_IN_reads.hic' sacCer3 > $outputdir/hic_files/log_hicfile_generation.txt"
cmd8="java -Xmx22g -jar $file_path_to_juicer_jar -r 10,50,100,150,200,500,1000,3200,5000 pre $outputdir/pairs_files/${sample_name}_output_OUT_reads.pairs $outputdir/hic_files/${sample_name}'_OUT_reads.hic' sacCer3 >> $outputdir/hic_files/log_hicfile_generation.txt"
cmd9="java -Xmx22g -jar $file_path_to_juicer_jar -r 10,50,100,150,200,500,1000,3200,5000 pre $outputdir/pairs_files/${sample_name}_output_SAME_reads.pairs $outputdir/hic_files/${sample_name}'_SAME_reads.hic' sacCer3 >> $outputdir/hic_files/log_hicfile_generation.txt"
cmd10="java -Xmx22g -jar $file_path_to_juicer_jar -r 10,50,100,150,200,500,1000,3200,5000 pre $outputdir/pairs_files/${sample_name}_output_noIN.pairs $outputdir/hic_files/${sample_name}'_noIN.hic' sacCer3 >> $outputdir/hic_files/log_hicfile_generation.txt"
echo $cmd7
time eval $cmd7
echo $cmd8
time eval $cmd8
echo $cmd9 
time eval $cmd9
echo $cmd10
time eval $cmd10

echo "Creating Distance Decay plots and contacts.txt"
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
cd $outputdir/distance_decay_plots/
cmd11="python $scriptsdir/distance_decay.py $outputdir/pairs_files/${sample_name}_output_IN_reads.pairs ${in_reads} $distance_graphed 'False'"
cmd12="python $scriptsdir/distance_decay.py $outputdir/pairs_files/${sample_name}_output_OUT_reads.pairs ${out_reads} $distance_graphed 'False'"
cmd13="python $scriptsdir/distance_decay.py $outputdir/pairs_files/${sample_name}_output_SAME_reads.pairs ${same_reads} $distance_graphed 'False'"
cmd14="python $scriptsdir/distance_decay.py $outputdir/pairs_files/${sample_name}_output_noIN.pairs ${noIN_reads} $distance_graphed 'False'"
echo $cmd11
time eval $cmd11
echo $cmd12
time eval $cmd12
echo $cmd13
time eval $cmd13
echo $cmd14
time eval $cmd14
mv $outputdir/pairs_file/${sample_name}_output_IN_readsori_decay.png $outputdir/distance_decay_plots
mv $outputdir/pairs_file/${sample_name}_output_OUT_readsori_decay.png $outputdir/distance_decay_plots
mv $outputdir/pairs_file/${sample_name}_output_SAME_readsori_decay.png $outputdir/distance_decay_plots
mv $outputdir/pairs_file/${sample_name}_output_noINori_decay.png $outputdir/distance_decay_plots

bgzip $outputdir/pairs_files/${sample_name}_output_IN_reads.pairs
bgzip $outputdir/pairs_files/${sample_name}_output_OUT_reads.pairs
bgzip $outputdir/pairs_files/${sample_name}_output_SAME_reads.pairs
bgzip $outputdir/pairs_files/${sample_name}_output_noIN.pairs

#creating cooler files
echo "Creating and balancing cooler files"
cd $outputdir/pairs_files/
cmd15="pairix -f ${sample_name}_output_IN_reads.pairs.gz"
cmd16="pairix -f ${sample_name}_output_OUT_reads.pairs.gz"
cmd17="pairix -f ${sample_name}_output_SAME_reads.pairs.gz"
cmd18="pairix -f ${sample_name}_output_noIN.pairs.gz"
cmd19="cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_IN_reads.pairs.gz $outputdir/cool_files/${sample_name}_output_IN_reads.cool"
cmd20="cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_OUT_reads.pairs.gz $outputdir/cool_files/${sample_name}_output_OUT_reads.cool"
cmd21="cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_SAME_reads.pairs.gz $outputdir/cool_files/${sample_name}_output_SAME_reads.cool"
cmd22="cooler cload pairix $file_path_to_chrom_sizes:150 ${sample_name}_output_noIN.pairs.gz $outputdir/cool_files/${sample_name}_output_noIN.cool"
cmd23="cooler zoomify $outputdir/cool_files/${sample_name}_output_IN_reads.cool"
cmd24="cooler balance $outputdir/cool_files/${sample_name}_output_IN_reads.cool"
cmd25="cooler zoomify $outputdir/cool_files/${sample_name}_output_OUT_reads.cool"
cmd26="cooler balance $outputdir/cool_files/${sample_name}_output_OUT_reads.cool"
cmd27="cooler zoomify $outputdir/cool_files/${sample_name}_output_SAME_reads.cool"
cmd28="cooler balance $outputdir/cool_files/${sample_name}_output_SAME_reads.cool"
cmd29="cooler zoomify $outputdir/cool_files/${sample_name}_output_noIN.cool"
cmd30="cooler balance $outputdir/cool_files/${sample_name}_output_noIN.cool"
echo $cmd15
time eval $cmd15
echo $cmd16
time eval $cmd16
echo $cmd17
time eval $cmd17
echo $cmd18
time eval $cmd18
echo $cmd19
time eval $cmd19
echo $cmd20
time eval $cmd20
echo $cmd21
time eval $cmd21
echo $cmd22
time eval $cmd22
echo $cmd23
time eval $cmd23
echo $cmd24
time eval $cmd24
echo $cmd25
time eval $cmd25
echo $cmd26
time eval $cmd26
echo $cmd27
time eval $cmd27
echo $cmd28
time eval $cmd28
echo $cmd29
time eval $cmd29
echo $cmd30
time eval $cmd30

######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> PAIRTOOLS VERSION:"
pairtools --version
echo -e "\n>>> COOLTOOLS VERSION:"
cooltools --version
echo -e "\n>>> PYTHON VERSION:"
python --version
echo -e "\n>>> BOWTIE2 VERSION:"
bowtie2 --version
echo -e ">>> END: Micro-C complete."
