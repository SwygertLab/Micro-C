###Tutorial

Use this tutorial alongside the pipeline readme.md files to understand how to use this script. 

##Pipeline_V2

##Configuring dependencies

You will need an environment with the following packages which can be found on [anaconda.org](https://anaconda.org/anaconda/conda):

-python version 3.8+
-cooltools
-pairtools
-numpy
-matplotlib
-deeptools

Ideally you will have a conda environment specific to Micro-C experiments. In order to create a conda environment, you will need an active version of anaconda and type in conda create --name (your_conda_env_name).

These are the commands and versions used at the time of creating this script.

conda install -c bioconda bowtie2

conda install python=3.8

conda install -c conda-forge matplotlib-base

conda install -c bioconda cooltools

conda install -c conda-forge -c bioconda pairtools=0.3.0 python=3.8

##Downloading files

The main script that runs the rest is in pipeline_V2 folder called MicroC_V2.sh. You will need that file, distance_decay.py, and separate_by_orientation.py which are in the essential_files folder. 
The last file required that does not vary based on sample or experiment is called juicer_tools_1.22.01.jar. This file version can be found and downloaded from the Aiden lab github here https://github.com/aidenlab/juicer/wiki/Download. 

As for the files required for the Micro-C analysis that are specific to your organism and experiment. They include a reference genome (or simply reference to the region of chromatin that was sampled from) and a file containing the size of the chromosomes. 
-There is an example of what the chromosome sizes file organization should look like in essential_files folder called example_sacCer3.chrSizes. If your organism is S. cerevisiae then that example is also a usable file for chromosome sizes for yeast. If you are unsure about the reference genome, there should be many resources on not only how to make it if you just google the organism name reference genome, but also there should be already made reference genomes available. They are typically .fa files. 

There is no required schema for where or how your files are located as the script requires that you edit it to include where the things it calls are located.
Although there is one required schema being that distance_decay.py and separate_by_orientation.py, wherever they are, are located in the same directory. 
I will go through edits need to be made to the script after finishing download details.

In /test_data/ there is a pair end set of Micro-C data. This is Log Phase W303 S. cerevisiae from Molecular Cell paper, Swygert et al, 2019 and truncated to only the first 100,000 lines of some our sample data as the whole thing is too large to upload.
Download this test data into a directory as you would like.

##Script Details

I have included in the first couple lines of the script more tutorial that is already covered here or will be covered here.

Look for the line that says ####### MODIFY THIS SECTION #############
From there on onwards are a couple things that need to be changed. The first thing is you should change the conda activate command to the name of your environment with the required packages. 

This script also requires a metadata file. An example of what that should look like is provided in this folder called metadata_example.txt and it includes the first set of your paired end data, the second set of your paired end data, and the prefix you would like to add to all the new files that will be made in the script. I personally add more data in there such as the date and location of all important files for that experiment, but that is all that is required. 

Then you just need to provide the file path to your reference genome file, chromosome sizes file, and juicer jar file. Personally, I keep these in a subdirectory away from the experiments with their own folders so that I can always use the same path in my different experiments. 

Then I include a line for the date and a place to put a path to your input (sample) directory, scripts directory, and output directory. While no particular schema is required and you could put all the files in the same folder and point to the same location for all these directories, I will list how I set it up.

I typically create an experiment directory and in there create an input directory where I tend to store the Micro-C experiment .fastq files and the metadata file, scripts directory with distance_decay.py, separate_by_orientation.py, and MicroC_V2.sh, and a output directory where the outputs from the script will be output to. 
-Note that whatever the output directory points to is where sevaral new directories is set to be created to sort the different types of output made from this pipeline. 

Part of the output of this pipeline are distance decay plots, if you know to what distance you would like the graphs to be made to then you can change the distance_graphed variable to that distance just note that the data is binned into 10 base pair bins so you divide the distance desired by 10 and that python does not include the number you set it to so you need to add 1 to where you want to go. 
-The distancy_decay.py script can also handle setting the axis of the graphs so you may place all the sample you would like on the same scale, but that will be covered in the readme.md of pipeline_V2

Also the mapq_filter variable can be changed as you like, it is a bit arbitrary; however, mine is set to 2 as I found that to be the minimum required filtering score to remove artifact from transposon regions of the sample and maintain the most amount of data. If you would like to find out more about how mapq filtering works you can do so here: https://samtools.github.io/hts-specs/SAMv1.pdf

-These are all the required changes to the script itself. More changes can be made, but do so at your own knowledge and understanding of this analysis. 

##Running the script

Assuming you have all your data downloaded and everything is ready to go, then here is an example command of how to run this script. 
sbatch /path/to/MicroC_V2.sh /path/to/metadatafile.txt

The script is constructed in such a way that if it fails, the error message will include where in the pipeline it failed to hopefully help you figure out what is wrong. One of the version commands also makes it harder for the computer to recognize the file type as it thinks there is non-text in there, but there isn't anything that makes it unopenable, so you can skip the warning message or # out the version output command at the bottom of the script. 

The results should look like the following for the sample data provided:

Micro-C heatmaps opened through juicebox-

Distance Decay Plots-

