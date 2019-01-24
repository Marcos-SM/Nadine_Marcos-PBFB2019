# Data information

*Project for the course Practical Bioinformatics for Biologists 2019*

*By Marcos Suárez Menéndez (S3684202) and Nadine Rijsdijk (S3813185)*

*Date: 24-01-2019*

*University of Groningen*

This folder contains the data of different samples (one folder per sample) plus the output of a trial run of the program. Each sample contains several ab1 files that could either be a forward or a reverse sequence. In this case, the samples are random ab1 files from the internet to
use as an example for the program.

* Parameters set during this run:
  * Do you want to extract fastq files from ab1 files
  (if not the script will work with existing fastq files)?: Y
  * Do you want to add the parent directory to the file name?: Y
  * Do you want to cut the fastq files?: Y
  * Point of cut 5' extreme: TGRAGG
  * Point of cut 3' extreme: GCCCGG
  * Do you want to get quality information of the fastq files?: Y
  * Do you want a Phred threshold value?: Y
  * Phred value: 50
  * Do you want to create a Fasta file with the best sequence from each sample?: Y
  * Enter species name separated by "\_":Prueba_prueba
  * Enter region of sequences:Test

  
* It also contains the program scripts and its outputs:
  * Log.txt: Information of how the program ran
  * Quality.txt: Mean qualities of the analysed files
  * bestqs.txt: Fastq with best mean quality per samples

* Inside each samples folder several files were created:
  * \*.fastq: extracted Fastq from the ab1 file
  * \*.qual: Quality data with Phred values of each Fastq
  * \*.html: Summary of the quality analysis done by FASTQC
  * \*.zip: Zip containing more information about the FASTQC analysis



*Sources of our example data:*

[Example data 1](http://www.applied-maths.com/download/sample-data)

[Example data 2](https://github.com/labsquare/CutePeaks/tree/master/examples)
