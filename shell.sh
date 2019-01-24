#! /bin/bash

#This function is used to be able to make a list of the sequence with the highest quality for each sample
#This is done by first replacing the forward slashes in the lines by tabs (sed -e)
#Then, the files are sorted by the first and third column (sample name and mean Phred quality), and the qualities are sorted form high to low
#With the uniq command the first line of each sample remains. This line contains the highest quality score.
#The tabs are again changed to forward slashes and the best sequence of each sample is saved in a new file.
grep -E "\w+"/*. Quality.txt | sed -e 's/\//\t/g' | sort -t $'\t' -u -k1,1 -k3,3r | uniq -w 8 | sed -E 's/(\w+)\t(\w+\.qual)/\1\/\2/g' > bestqs.txt
