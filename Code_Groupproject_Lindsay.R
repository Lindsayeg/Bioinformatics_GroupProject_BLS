#Pseudocode to translate nucleotide sequences to aminoacid sequences

#Read the aa sequence file - codonmap.txt file
#Read the nucleotide fasta files control 1&2 and obsese 1&2
#Translate nucleotide fasta files into aa fasta files
#Reading the nucleotide fasta into aa fasta - start with startcodon and end at first stopcodon - use sed?(manual?)

#Set working directory and read each of the original neucleotide fasta files.
setwd("~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/")
Codon_Map <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/codonmap.txt', what = , sep="\n", character())
Control1 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/control1.fasta', what = , sep="\n", character())
Control2 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/control2.fasta', what = , sep="\n", character())
Obese1 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/obese1.fasta', what = , sep="\n", character())
Obese2 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/obese2.fasta', what = , sep="\n", character())


#Load package
library('stringr')
library('Biostrings')
library('devtools')
library('seqinr')


#Read in fasta file, each line becomes an item in a vector
Control1=scan("control1.fasta",what=character(),sep="\n")

###Determine the ORF of each sequence file using a for loop.
for (Line in 1:length(Control1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(Control1[Line],">")){
    #Find the ORF in each sequence line and write to a new file
    Control1_ORF = str_extract(Control1[Line],"ATG([ATGC]{3})+(TAG|TAA|TGA)?")
  }
}

###Make a directory of amino acids and thier names.
#Make a new list called amino_acid_list with the length of Codon_Map.
amino_acid_list <- vector(mode="list", length=length(Codon_Map))

#Use a regex to put all nucleotide codon sequences into it. 
Codons <- str_extract(Codon_Map,"[A-Z]{3}")

#Set the names of the amino acid list using Codons. 
names(amino_acid_list) <- Codons

#Then use a regex to extract the amino acid names and put it in a new variable.
amino_acids <- str_extract(Codon_Map,"[A-Za-z]{1,4}")

#Use a for loop to combine the nucleotide sequences with the amino acid names.
for (i in 1:length(amino_acid_list)){
  amino_acid_list[[i]] <- amino_acids[i]
}


  
### Stil need to apply amino_acid_list to the mRNA sequences to get the protein sequences.
### Use a while loop (while i is less than or equal to the length of the reading frame, i+2)

Control1_Protein = character(length(Control1_ORF))
i <- 1
while (i < length(Control1_ORF)){
  Control1_Protein = c(Control1_Protein, amino_acid_list+(Control1_ORF[i:i+2]))
  print[i]
  i=i+3
}






###This is all in UNIX on the command line.
###Align reference sequence files using muscle and then build an HMM profile using hmmbuild.

#Set working directory
cd /Users/lou/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Protein_FASTA_Files

###Align reference sequence files using muscle and a for loop.
for file in ./*.fasta
do
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in $file -out $file.align
done

###Build a profile HMM using alignment sequences and hmmbuild using a for loop.
for file in ./*.align
do
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild $file.hmm $file
done








####Boneyard: Please IGNORE!
###Align reference sequence files using muscle either format below is okay! Repeat alignment for all files.

#Transcript 1
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript1_Gsta2.fasta -out Transcript1_Gsta2.align
#Transcript 2
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript2_Slc7a12.fasta -out Transcript2_Slc7a12.align
#Transcript 6
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript6_Ptpn5.fasta -out Transcript6_Ptpn5.align
#Transcript 8
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript8_Atp12a.fasta -out Transcript8_Atp12a.align
#Transcript 9
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript9_Lhx2.fasta -out Transcript9_Lhx2.align
#Transcript 10
/Users/lou/Desktop/muscle3.8.31_i86darwin64 -in Transcript10_Synpr.fasta -out Transcript10_Synpr.align



###Build a profile HMM using alignment sequences and hmmbuild.
#Transcript 1
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript1_Gsta2.hmm Transcript1_Gsta2.align
#Transcript 2
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript2_Slc7a12.hmm Transcript2_Slc7a12.align
#Transcript 6
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript6_Ptpn5.hmm Transcript6_Ptpn5.align
#Transcript 8
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript8_Atp12a.hmm Transcript8_Atp12a.align
#Transcript 9
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript9_Lhx2.hmm Transcript9_Lhx2.align
#Transcript 10
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmbuild Transcript10_Synpr.hmm Transcript10_Synpr.align
