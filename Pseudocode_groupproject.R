###Pseudocode to translate nucleotide sequences to aminoacid sequences
##Load the required packages
#Read the aa sequence file - codonmap.txt file
#Read the nucleotide fasta files control 1&2 and obsese 1&2
#Translate nucleotide fasta files into aa fasta files
#Reading the nucleotide fasta into aa fasta - start with startcodon 
#and end at first stopcodon - using regex and translate (Bhavana)
#-using grep and while loop (Lindsay)
#convert the translated file into fasta format


###Pseudocode for building HMM profiles using Unix
#script for muscle alignment
#forloop for muscle alignment
#script for building HMM profile
#for loop HMM profile Build
#script for hmmsearch
#for loop for hmmsearch of four RNAseq files and 6 transcript files
#Tabulate the counts fer each protein



###Pseudocode for tophits table and identify 6 differentially expressed transcripts
#Use BLAST to identify the genes encoded by the 6 differentially expressed transcripts
#listed in uniquetranscripts.fasta.
#Using Unix commands, make a single table that includes the top hit for each
#transcript.
#Search the NCBI protein database for amino acid sequences corresponding to
#these 6 transcripts.

<<<<<<< HEAD
##Further exploration 1
#Use BLAST to identify the genes encoded by the 3 differentially expressed transcripts
#Using Unix commands, make a single table that includes the top hit for each
#transcript.
#Search the NCBI protein database for amino acid sequences corresponding to
#these 3 transcripts.

##Further exploration 2
#reused the hmm build script
#identify distantly related proteins (less conserved)
#identify the expression counts
#tabulate the differential counts

=======
###Psuedocode for Further Exporation 1
#Use the same methods to obtain BLAST hits as given in instructions and used above on the uniquetranscripts.fasta, except
#make the reported sequences restricted to mouse sequences to make the task more 
#manegable. Use all three algorithims (megablast, discontinuous megablast,
#and blastn) which can be found in the "Optimize for" option (as mentioned in the instructions).
>>>>>>> 5cd300383b07837cb85a1766f4cbeff5af0a0386
