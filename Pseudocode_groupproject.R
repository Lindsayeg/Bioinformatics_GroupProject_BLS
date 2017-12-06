###Pseudocode to translate nucleotide sequences to aminoacid sequences

#Read the aa sequence file - codonmap.txt file
#Read the nucleotide fasta files control 1&2 and obsese 1&2
#Translate nucleotide fasta files into aa fasta files
#Reading the nucleotide fasta into aa fasta - start with startcodon 
#and end at first stopcodon - using regex and translate


###Pseudocode for building HMM profiles
#script for muscle alignment
#forloop for muscle alignment
#script for building HMM profile
#for loop HMM profile Build
#script for hmmsearch
#for loop for hmmsearch of four RNAseq files and 6 transcript files


###Pseudocode for tophits table and identify 6 differentially expressed transcripts
#Use BLAST to identify the genes encoded by the 6 differentially expressed transcripts
#listed in uniquetranscripts.fasta.
#Using Unix commands, make a single table that includes the top hit for each
#transcript.
#Search the NCBI protein database for amino acid sequences corresponding to
#these 6 transcripts.