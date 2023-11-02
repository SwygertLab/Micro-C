###Tutorial

In /test_data/ there is a pair end set of Micro-C data. This is Log Phase W303 S. cerevisiae from Molecular Cell paper, Swygert et al, 2019.

You will need an environment with the following packages which can be found on [anaconda.org](https://anaconda.org/anaconda/conda):

-python version 3.10+
-cooltools
-pairtools
-numpy
-matplotlib
-deeptools

Download this test data into a directory as you would like but include filter_orientation_heading.py and distance_decay.py in the same directory.

So you should have a directory that includes the following files:

- LogPhase_W303_SCR1.fastq
- LogPhase_W303_SCR2.fastq
- filter_orientation_heading.py
- distance_decay.py

You do need MicroC_pipeline.sh as well, it just does not need to be in the same directory.

Once this is complete while in the directory with your fastq files run a command like this:

sbatch /path/to/MicroC_pipeline.sh LogPhase_W303_SCR1.fastq LogPhase_W303SCR2.fastq LogPhase_W303SC

The last parameter here is just the prefix given to the new files generated so feel free to change as you would like.

The results should look like the following:

Micro-C heatmaps opened through juicebox-




Distance Decay Plots-

