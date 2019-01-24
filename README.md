# Extracting electropherogram information and manipulating FASTQ files

*Project for the course Practical Bioinformatics for Biologists 2019*

*By Marcos Suárez Menéndez (S3684202) and Nadine Rijsdijk (S3813185)*

*Date: 24-01-2019*

*University of Groningen*

## Basic information to run the program

* The program contains two scripts: *Project.py* and *shell.sh*

* Remember to make the scripts executable for the user.

* Both of the scripts need to be run from the folder where you have the data (the folder containing the directories with
the data).

* We recommend removing the writing rights to your original data.

* Biopython v. 1.73 and FastQC need to be downloaded before running this script.

* Our sample names have 8 characters. If your sample name has a different number of characters, change the 8 in
'uniq -w 8' in the shell.sh script to the number of characters your sample has.

* Clean.sh will delete all the files that are created by the program (to be run from the folder where the data is).
**Do not use this script if there are more files in your folders than only the ab1 files.**

## Explanation of the script

The script provided is a Python (v. 2.7.15rc1) script that extracts information from ab1 files and converts it to fastq
files. If you want, it adds the name of the parent directory to the file name. Then, the sequences are cut at specific
points to obtain the part of the sequence you want to keep. This is done by searching the 5' and 3' extreme of the
sequence you want to keep, and cutting off the outsides of this part. If the program is unable to cut the sequence, it
will look for the reverse by looking for the complementary sequence of the 5' end extreme of the sequence. If it finds
this, it will reverse the whole sequence and the line with the corresponding quality scores in the fastq file. Then, it
will cut the reverse sequence.

Next, quality analysis is done on the files. If the sequences were cut, this analysis is only done on the cut files.
If the sequences were not cut, the analysis is done on all the files. The fastq files are first converted to qual files.
In this file, the quality score (in ASCII) for each nucleotide in the FASTQ file is converted to a Phred quality score
for each nucleotide. Then, the mean Phred quality score for each sequence is calculated. If a Phred threshold is entered
by the user at the prompt, the sequences with a mean Phred quality score lower than this threshold are not taken into
account for further analysis. A text file with the path to the samples and their mean Phred quality score is created.

For each file, a FastQC report is created (html), which also contains a graph with Phred quality score plotted against
nucleotide position.

The sequence with the highest mean Phred quality of each sample is filtered from all the samples (using a shell script;
*shell.sh*). A final FASTA file is created which contains the best sequences. Each sequence has a header with
information about the sequence.

All the steps and the progress is logged in a text file.

*If there is an error in the cutting points in the sequence itself (at the 5' and 3' extreme)*
*(for example 'W' instead of 'A') the program won't be able to cut the sequence. The sequence name will then be printed*
*in the log so it can be manually checked later. Also, the length of the sequence will be printed, so it is visible*
*whether the sequence is incomplete.*

**The following questions are asked to the user before the script is run:**

1. Whether you want to extract fastq files from ab1 files (if not the script will work with your existing fastq
files directly).
2. Whether you want to add the name of the parent directory to the file name (this should be done only once! Otherwise
your file names will become very long each time you run the script).
3. Whether you want to cut the fastq files to get a certain part of the sequence. If yes, the program will ask
you to specify the beginning and ending parts of the sequence, so the program can cut the sequence.
4. Whether you want to get the quality information of the fastq files. If yes, the program will ask you whether
you want to use a Phred threshold value for your sequences, and if so, it will ask you to specify this threshold.
5. Whether you want to create a fasta file that includes the sequence with the highest quality from each sample. If yes,
it will ask you to enter the name of the species you are analysing, and the region in the sequence you are analysing.

![rug](https://www.rug.nl/_definition/shared/images/logo--en.png)
