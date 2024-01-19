# MicC_pipeline.sh
The following is a pipeline designed for usage in data analysis of Micro-C data. 

A couple packages are needed to run this code all of which are downloadable from https://anaconda.org/anaconda/conda including:
- python version 3.8+
- cooltools
- pairtools
- numpy
- matplotlib

-Note that if you are new to this I highly recommend creating an environment for package installations of specific projects. In other words an environment
for a chip-seq experiment and a different one for a micro-c experiment to prevent package conflicts. For CSU on alpine you need to run this command before
you are able to run conda commands "source /curc/sw/anaconda3/latest".

The main focal point of the pipeline begins with MicC_pipeline.sh and if properly working should be the only interface needed. 
This code was designed with simplicity for the user in mind, and in parallel with this running the MicC_pipeline_V3.sh with your .fastq files
should be all that is needed. 

example: sbatch MicC_pipeline_V3.sh your_microc_dataR1.fastq your_microc_dataR2.fastq name_for_newly_made_files

In order, for the java jar commands to work you need to download a jar file from the aiden lab github here -> https://github.com/aidenlab/juicer/wiki/Download 

There are a couple assumptions made by the MicC_pipeline.sh pipeline that should be changed in the batch script. These are lines 21-26. 

In line 21 you need your genome data. In line 22 you need your chromosome sizes data. In line 23 you need to point to where your jar file is In line 24 you need 
the distance you want to get from your distance decay plots.In line 25 you need to put in the prefix you would like for your new files to have. I usually set this 
to $3 and write it in the command to this script but it can be manually set if you would like. As an example for the fastq files LogWT_01.fastq and LogWT_02.fastq 
I made my prefix "LogWT". Line 26 records a preset mapq of the users choosing. I would recommend figuring out what mapq score is the best for your dataset 
without using a mapq score filter and then testing a lot of different mapq scores until you find one that removes the most amount of artifact but leaves behind as 
much real data as possible.

The MicC_pipeline.sh pipeline also contains some python scripts that need to be pointed to on line 56 and lines 73-76. You can do so like this:
python /path/to/python/script/separate_by_orientation.py ${sample_name}_output.pairs
These are part of this repository and are:
- separate_by_orientation.py
- distance_decay.py

distance_decay.py is called inside the pipeline and will output a graph containing the average short contacts of the different orientations.

distance_decay.py has additional functionality intended to be used outside of the main script. 
Essentially once you have multiple graphs you may want to see them on the same scale so you can directly compare them. 
To do this you can call the python file accordingly.

"py distance_decay.py name_of_your_pairs_file.pairs total_reads_in_pairs_file the_distance_you_want_graphed True x_limit1 x_limit2 y_limit1 y_limit2"

Note that in the script the "True" value above is set to False and therefore the distance decay plots generated are not in the same scale. 

The parameters are as follows a file where pairs should be located, a number that should indicated total reads in that specific reads orientation file, 
a number indicating how far you want the distance decay plot to graph, a string that needs to be either True or False (case sensitive) and it needs to 
be True if you want to set the limit, and the last 4 are all numbers indicating what limits you are setting. 

overlap_dist_decay.py is not called inside the pipeline and is there if you need to compare two different micro-c data sets. 

There are some additional scripts here that are available for other uses.

Starting with:

overlap_distance_decay.py

This script is intended to take two separate datasets of micro-c data in a .pairs file and simply creates an overlapping graph of them.

Can be used like this: python overlap_distance_decay.py sample1.pairs sample2.pairs

You can download the packages using these installation commands:

conda install -c bioconda bowtie2

conda install python=3.8

conda install -c conda-forge matplotlib-base

conda install -c bioconda cooltools

conda install -c conda-forge -c bioconda pairtools=0.3.0 python=3.8
