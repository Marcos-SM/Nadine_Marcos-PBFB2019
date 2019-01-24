#! /usr/bin/env python2
import Bio
import datetime
import re
from Bio import SeqIO
from Bio.Seq import Seq
import os
from os import listdir
from Bio.Alphabet import IUPAC
import subprocess



#Authors: Marcos Suarez Menendez and Nadine Rijsdijk
#Date: 24-01-2019
#This script was written for the course PBfB2019 at the University of Groningen

#This program allows you to:
#Extract a fastq file from each ab1 file (if you want to)
#Modify the names of ab1 files to include their parent directory (should be done only once!)
#Cut the fastq (both the sequence and the quality) in a specific point of the sequence (allows for polymorphisms in this sequence), target sequence included
#Specify a Phred threshold value for your sequences
#In case the program is not able to cut the sequence, it will check if it is the reverse sequence, reverse it and cut it
#Makes a FastaQC report
#Makes a fasta file with the sequence with the highest quality for each sample and adds information to the header

log = open("log.txt", "a")
now = datetime.datetime.now()
log.write("Date and time the script was run: %s\n" % now) #Date and time the script was executed

#Input questions of how to run the program and to specify several parameters
Abifiles = raw_input("Do you want to extract fastq files from ab1 files \n(if not the script will work with existing fastq files)?: (Y/N) ")
log.write("Do you want to extract fastq files from ab1 files \n(if not the script will work with existing fastq files)?: %s\n" % Abifiles)
if Abifiles is 'Y':
    Names = raw_input('\033[1m' + '\033[91m' + "Do you want to add the parent directory to the file name? \n[Warning: Do this only once]: (Y/N)"  + '\x1b[0m')
    log.write("Do you want to add the parent directory to the file name?: %s\n" % Names)
elif Abifiles is 'y': #Only shows this at the prompt when the user types a lowercase 'y'
    print "Answers should be in uppercase"
    os.remove("log.txt")
    quit()
elif Abifiles is 'n':
    print "Answers should be in uppercase"
    os.remove("log.txt")
    quit()

Cutting = raw_input("Do you want to cut the fastq files?:(Y/N) ")
log.write("Do you want to cut the fastq files?: %s\n" % Cutting)
if Cutting is 'Y':
    sequence5 = raw_input("Point of cut 5' extreme:  ").upper()
    sequence3 = raw_input("Point of cut 3' extreme:  ").upper()
    log.write("Point of cut 5' extreme: %s\nPoint of cut 3' extreme: %s\n" % (sequence5, sequence3))
Qualanalysis = raw_input("Do you want to get quality information of the fastq files?:(Y/N) ")
log.write("Do you want to get quality information of the fastq files?: %s\n" % Qualanalysis)
if Qualanalysis is 'Y':
    Threshold = raw_input("Do you want a Phred threshold value?:(Y/N) ")
    log.write("Do you want a Phred threshold value?: %s\n" % Threshold)
    if Threshold is 'Y':
        Phred = raw_input("Phred value: ")
        log.write("Phred value: %s\n" % Phred)
    Fasta = raw_input("Do you want to create a Fasta file with the best sequence from each sample?:(Y/N) ")
    log.write("Do you want to create a Fasta file with the best sequence from each sample?: %s\n" % Fasta)
    if Fasta is 'Y':
        Species = raw_input('Enter species name separated by "_": ')
        Region = raw_input('Enter region of sequences: ')
        log.write('Enter species name separated by "_":%s\nEnter region of sequences:%s\n' % (Species, Region))
exit = raw_input("Press Enter to continue... Type E to exit: ")
if exit is 'E':
    os.remove("log.txt") #The log file that was created will be deleted if the program is exited before running it
    quit()

log.write("\n----------------------\n\n")
print "\n\n----------------------\n\n"


#This modifies the file name of all the ab1 files to add the name of the parent directory (sample name):
if Names is 'Y':
    for root, dirs, files in os.walk("./"):
        if not files:
            continue
        prefix = os.path.basename(root)
        for f in files:
                path = os.path.relpath(os.path.join(root, f), "./")
                if path.endswith(".ab1"): #This specifies the extension of the files to be changed
                    os.rename(os.path.join(root, f), os.path.join(root, "{}_{}".format(prefix, f)))

#Function that creates a list with the path to the files with the specified extension and counts the number of files
def pathtofiles(Allfiles, arg1):
    count = 0
    for root, dirs, files in os.walk("./"):
        for f in files:
            path = os.path.relpath(os.path.join(root, f), "./")
            if path.endswith("%s" % arg1):
                Allfiles.append(path)
                count += 1
    log.write("Number of files to be processed:%d and paths to files:\n\n" % count)
    for f in Allfiles:
        log.write("%s\n" % f) #Write the path to the file to the log

if Abifiles is 'Y':
    print "Extracting fastq from ab1 files...\n"
    Allfiles = list() #Necessary for the next function
    pathtofiles(Allfiles, arg1='.ab1')
    #Creates a list with the output of the conversion from ab1 to fastq
    fastq = [ext.replace(".ab1", ".fastq") for ext in Allfiles] #Change extension
    log.write("\n---------\n\nExtraction of fastq from ab1 files:\n")
    converted = 0
    converted_errors = 0
    for filename, output in zip(Allfiles, fastq): #Loops through both lists at the same time
        try:
            SeqIO.convert(filename, "abi", output, "fastq-sanger") #Convert ab1 to fastq
            print ("Converted %s to %s") % (filename, output)
            log.write("Converted %s to %s\n" % (filename, output))
            converted += 1
        except: #If there is an error converting the ab1 files:
            print ("Error converting %s to %s") % (filename, output)
            log.write("Error converting %s to %s\n" % (filename, output))
            os.remove(output)
            converted_errors += 1
            pass
    log.write("\n\nSummary of conversion: \nNumber of files converted: %d\n" % converted)
    log.write("Number of errors during conversion: %d\n" % converted_errors)
    log.write("\n----------------------\n\n")
    print "\n\n----------------------\n\n"


Fastaqfiles = list() #This list is necessary for the function
pathtofiles(Fastaqfiles, arg1='.fastq')

log.write("\n----------------------\n\nCut fastq files:\n")

if Cutting is 'Y': #Cutting sequences in specific points (also cutting the quality in those points)
    print "Cutting fastq files...\n"
    #Cut sequences
    #Make a dictionary for the nucleotides. Also the ambiguity codes are put in the dictionary
    #For example, if you put a 'R' in your 5' cutting point, the program will look if there is either an 'A' or a 'G' in that position
    DIUPAC={
    'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'R':'[AG]',
    'Y':'[CT]',
    'S':'[GC]',
    'W':'[AT]',
    'K':'[GT]',
    'M':'[AC]',
    'B':'[CGT]',
    'D':'[AGT]',
    'H':'[ACT]',
    'V':'[ACG]',
    'N':'[ATGC]',
    }
    #Lists necessary for converting ambiguous nucleotides
    sequence5f = list() #5' cutting point
    sequence3f = list() #3' end cutting point
    sequence5c = list() #5' complementary cutting point

    #Creating the 5' and 3' cutting points
    for nucleotide in sequence5:
        sequence5f.append(DIUPAC[nucleotide]) #Changes the nucleotide key for its value
    sequence5f = "".join(sequence5f) #Joins the list into a string
    for nucleotide in sequence3:
        sequence3f.append(DIUPAC[nucleotide])
    sequence3f = "".join(sequence3f)
    #Obtaining complementary of 5' target sequence (cutting point) to be able to check whether a sequence is the forward or the reverse
    check = Seq(sequence5, IUPAC.IUPACUnambiguousDNA())
    checkseq =  check.reverse_complement()
    for nucleotide in checkseq:
        sequence5c.append(DIUPAC[nucleotide])
    sequence5c = "".join(sequence5c)

    log.write("Target string used to cut 5' end: %s\nTarget string used to cut 3' end: %s\n\n" % (sequence5f, sequence3f))
    Preverse = list() #Necessary for cutting
    Cutfiles = list()

    def cutfile(filename, sequence5f, sequence3f): #Function to cut the fastq files in the indicated points
        global cut #The variable can also be used outside the function
        global uncut
        with open(filename, "r") as filec:
            #The sequence you want to obtain is captured, and the parts that come before and after are captured
            SearchString ='(\w+)(%s\w+%s)(.+)' % (sequence5f, sequence3f)
            lines = filec.readlines()
            Result = re.search(SearchString, lines[1])
            if Result != None: #Only does it if it finds the first cutting point in the sequence
                region1 = len(Result.group(1))
                region2 = len(Result.group(2)) #The region that will remain after the sequence will be cut
                region3 = len(Result.group(3))
                print "Fastq %s cut, final size: %d " % (filename, region2)
                log.write("Fastq %s cut, final size: %d\n" % (filename, region2))
                with open(filename, "w") as filec : #Rewrites the file after cutting
                    seq = re.sub(SearchString, r'\2', lines[1])
                    qualitycut = re.sub(r'.{%d}(.+).{%d}\n' % (region1, region3), r'\1\n', lines[3])
                    text = lines[0] + seq + lines[2] + qualitycut
                    filec.write(text)
                    Cutfiles.append(filename)
                    cut += 1
            else:
                print "Sequence of %s could not be cut" % filename
                log.write("Sequence of %s could not be cut\n" % filename)
                Preverse.append(filename)
                uncut += 1
        filec.close()

    cut = 0
    uncut = 0
    for filename in Fastaqfiles:
        cutfile(filename, sequence5f, sequence3f)
    print "\nResult of initial cutting:\n%d files were cut" % cut
    log.write("\nResult of initial cutting:\n%d files were cut\n" % cut)
    log.write("%d files were not cut\n\n-----------\n\n" % uncut)
    print "%d files were not cut\n-----------\n" % uncut
    print "Checking if uncut files are reverse sequences and in that case, cutting them...\n"
    log.write("Checking if uncut files are reverse sequences and in that case, cutting them:\n\n")

    #This checks if a sequence is the reverse by looking for the complementary sequence of
    #the 5' target input and if it finds it, it reverses the sequence
    areverse = list()
    reverse = 0
    unreversed = 0
    uncomplete = list()
    for filename in Preverse:
                with open(filename, "r") as filec :
                    lines = filec.readlines()
                    sequence = lines[1] #Sequence is on the second line (lines start with 0)
                    qual = lines[3][::-1][1::] #Reverse the line with the qualities, start with '1' to skip the end-of-line character at the beginning
                    seq = Seq(sequence, IUPAC.IUPACUnambiguousDNA()) #Create object seq
                    seqreverse = seq.reverse_complement()
                    seqfinal = seqreverse[1:] #Skips end-of-line character
                    SearchString = sequence5c
                    Result = re.search(SearchString, str(lines[1]))
                    if Result != None:
                        with open(filename, "w") as filec :
                            text = "%s%s\n%s%s\n" %(lines[0], seqfinal, lines[2], qual) #Rewrite the file so the sequence and qualities are reversed
                            filec.write(text)
                            areverse.append(filename)
                            log.write("Sequence of file %s reversed\n" % filename)
                            print "Sequence of file %s reversed" % filename
                            reverse += 1
                    else:
                        print "Sequence %s not complete" % filename
                        unreversed += 1
                        uncomplete.append(filename)
                        log.write("Sequence %s not complete\n" % filename)
                filec.close()

    print "\n\nTrying to cut reversed sequences:\n"
    log.write("\nFiles reversed: %d\nFiles not reversed: %d\n\n" % (reverse, unreversed))
    log.write("Trying to cut reversed sequences:\n\n")

    #This loops tries to cut the fastq files that were reversed in the previous step
    Preverse = list() #Setting this list blank
    for filename in areverse: #Cutting the reverse files
        cutfile(filename, sequence5f, sequence3f)


    #This gets a list of the files containing a sequence that couldn't be cut anyway and provides the length of the sequence
    log.write("\nSequence length of files that couldn't be cut:\n\n")
    print "\nSequence length of files that couldn't be cut:"
    filerrors = uncomplete + Preverse
    for filename in filerrors:
        with open(filename, "r") as filec:
            lines = filec.readlines()
            lengthseq = len(lines[1])
            log.write("%s : %d\n" % (filename, lengthseq))
            print "%s : %d" % (filename, lengthseq)
        filec.close()


#If the files were cut the quality analysis is done on the cut files
#If they weren't cut, then quality analysis is done on all the files
if Qualanalysis is 'Y':

    if Cutting is 'Y':
        Listoffiles = Cutfiles
    else:
        Listoffiles = Fastaqfiles

    log.write("\nQuality analysis:\n\n")
    qual = [ext.replace(".fastq", ".qual") for ext in Listoffiles]

    for filename, output in zip(Listoffiles, qual):
        SeqIO.convert(filename, "fastq", output, "qual") #Converting fastq files to qual files

    #Define the list Qualfiles
    Qualfiles = list()
    pathtofiles(Qualfiles, arg1='.qual')

    #Open a file that holds the path to the sample \t meanquality
    #Want to append to the file, + means: if the file doesn't exist already, create it
    Quality = open("Quality.txt", "a+")

    def fastqc(file): #Function that runs the program FastQC in the shell for a specified file
        command = "fastqc %s" % file
        subprocess.check_call(command, shell=True)
        log.write("%s analysed\n" % file)

    #Creating a list of the mean Phred quality for each sequence
    print "\nQuality analysis:\n"
    log.write("\nQuality analysis:\n\n")
    for item,file in zip(Qualfiles, Listoffiles):
        for record in SeqIO.parse(item, "qual"):
            Quallist = (record.letter_annotations["phred_quality"])
        # Coverting the integers in the list to floats:
        FloatQual = [float(Element) for Element in Quallist]
        # Calculating mean quality for each sequence:
        MeanQual = [sum(FloatQual) / len(FloatQual)]
        if Threshold is 'Y': #If the user entered a Phred threshold value at the prompt
            if MeanQual[0] > int(Phred):
                Quality.write("%s\t%d\n" % (item, MeanQual[0]))
                print "%s\t%d" % (item, MeanQual[0])
                log.write("%s\t%d\n" % (item, MeanQual[0]))
                fastqc(file) #A FastQC report is made for the file
            else:
                print "Quality mean value of %s (%d) lower than threshold %s" % (item, MeanQual[0], Phred)
                log.write("Quality mean value of %s (%d) lower than threshold %s\n" % (item, MeanQual[0], Phred))
        else:
            Quality.write("%s\t%d\n" % (item, MeanQual[0]))
            print "%s\t%d" % (item, MeanQual[0])
            log.write("%s\t%d\n" % (item, MeanQual[0]))
            fastqc(file) #A FastQC report is made for the file

    Quality.close()


    if Fasta is 'Y':
        #Finds the sequence with the highest quality for each sample
        #Call the shell script
        subprocess.call(['./shell.sh'])

        #Fasta file with best sequences per sample
        #Joins sequences with best quality per sample into a fasta file
        with open("bestqs.txt", "r") as filec :
            content = filec.read()
            filec.close()
            log.write("\n\nBest fastq from each sample:\n\n%s\n" % content)
            fastqpath = re.sub(r'([^.]+).+', r'\1.fastq', content) #Gets paths to fastq files
            qualvalue = re.sub(r'.+\t(\d+)', r'\1', content) #Gets mean Phred quality values
        fastalist = [y for y in (x.strip() for x in fastqpath.splitlines()) if y] #Converts the previous variables to a list
        qualvaluef = [y for y in (x.strip() for x in qualvalue.splitlines()) if y]

        #Making the final fasta file
        final = open("final.fasta", "a")
        for filename, qualv in zip(fastalist, qualvaluef):
            with open(filename, "r") as filec:
                lines = filec.readlines()
                sample = re.sub(r'([^\/]+).+', r'\1', str(filename))
                header = re.sub(r'\@(.+)', r'>\1_%s_%s_%s_Q_%s' % (sample, Species, Region, qualv), lines[0]) #Modifies header
                final.write(header + lines[1]) #Writes the header and the sequence in the fasta file (on different lines)
            filec.close()
        print "\nBest sequence from each sample joined in one fasta file"
        log.write("\nFasta sequences of from these fastq joined in final.fasta\n\n")
log.write("Finished!")
log.close()
