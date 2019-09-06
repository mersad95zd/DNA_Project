# Intro
The aim of this project is to recover the data stored on DNA sequences from millions of noisy reads. In this regard, we first encoded a sample zip file, then induce errors and erasures to the encoded DNA sequences, then cluster the reads using Locality-Sensitive Hashing method, then by performing Multiple Sequence Alignment and Majority Voting on each cluster we generate a number of candidates and by putting them into a decoder we receover the original file. A step-by-step approach for an example is provided below.

==============

# Setup and Installation
Assuming the experiment is being performed in a docker container, the following libraries and packages need to be installed.

        apt-get update
        apt-get install gcc
        apt-get install make
        apt-get install git
        apt-get install libboost-all-dev
	apt-get install python3.6

	apt-get install python-numpy
	apt-get install python-biopython
	apt-get install python-sklearn

 
The encoder/decoder I used is a Reed-Solomon encoder/decoder written by Reinhard Heckel. Download https://github.com/reinhardh/dna_rs_coding as a zip file and extract it, then put the LSH_clustering.ipynb jupyter notebook in the main directory where you extracted that zip file.

For multiple sequence alignment, I used MUSCLE command line. Please download the proper version of the software from http://www.drive5.com/muscle/downloads.htm and put it in the main directory. LSH_clustering.ipynb is written for linux, so if you are using another operating system, after downloading the proper version of MUSCLE software and putting it in the main directory, you need to change "muscle_exe" in the "multiple_alignment_muscle" function (cell #12 in LSH_clustering.ipynb) to the name of the file you dowloaded.

# Run and enjoy!
Now, all you need to do is to open and LSH_clustering.ipynb or simply run the following code in the command line.

`jupyter notebook LSH_clustering.ipynb`
