# swygertlab
The following is a pipeline designed for usage in data analysis of Micro-C data. 

A couple packages are needed to run this code all of which are downloadable from https://anaconda.org/anaconda/conda including:
- python version 3.10+
- cooltools
- pairtools
- numpy
- matplotlib
- deeptools

-Note that if you are new to this I highly recommend creating an environment for package installations of specific projects. In other words an environment
for a chip-seq experiment and a different one for a micro-c experiment to prevent package conflicts. For CSU on alpine you need to run this command before
you are able to run conda commands "source /curc/sw/anaconda3/latest".

The main focal point of the pipeline begins with MicC_pipeline.sh and if properly working should be the only interface needed. 
This code was designed with simplicity for the user in mind, and in parallel with this running the MicC_pipeline_V3.sh with your .fastq files
should be all that is needed. 

example: sbatch MicC_pipeline_V3.sh your_microc_dataR1.fastq your_microc_dataR2.fastq pick_a_title_for_your_data
-Note that a 3rd parameter is taken to name the future files that are being outputted by the script feel free to type in any title that makes sense to you.

distance_decay.py is called inside the pipeline and will output a graph containing the average short contacts of the different orientations.

distance_decay.py has additional functionality intended to be used outside of the main script. 
Essentially once you have multiple graphs you may want to see them on the same scale so you can directly compare them. 
To do this you can call the python file accordingly.

"py distance_decay.py name_of_your_pairs_file.pairs total_reads_in_pairs_file the_distance_you_want_graphed True x_limit1 x_limit2 y_limit1 y_limit2"

The parameters are as follows a file, a number, a number, a string that needs to be either True or False (case sensitive) and it needs to be True if you 
want to use set the limit, and the last 4 are all numbers indicating what limits you are setting. 

overlap_dist_decay.py is not called inside the pipeline and is there if you need to compare two different micro-c data sets. 

There are a couple assumptions made by the MicC_pipeline.sh pipeline mostly about the absolute paths of files, and can be changed in the 
batch script itself if desired or simply follow the file directory scheme.

Firstly, the genome of your organism is in a folder called Saccer3. # See line 17

Secondly, the size of each chromosome is in a folder called microc and the file is called sacCer3.chrSizes. # See lines 18

In order, for the java jar commands to work you need to download a jar file from the aiden lab github here -> https://github.com/aidenlab/juicer/wiki/Download 
this file is also assumed to be in the working directory. 

The MicC_pipeline.sh pipeline also assumes that the associated python scripts are in the working directory. 
These are part of this repository and are:
- filter_pairs.py
- distance_decay.py

There are some additional scripts here that are available for other uses.

Starting with:

overlap_distance_decay.py

This script is intended to take two separate datasets of micro-c data in a .pairs file and simply creates an overlapping graph of them.

Can be used like this: python overlap_distance_decay.py sample1.pairs sample2.pairs
