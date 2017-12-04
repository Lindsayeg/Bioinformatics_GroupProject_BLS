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

#install Bioconductor and biostrings
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

install.packages("devtools")
install.packages("dataMeta")

#Load package
library('stringr')
library('Biostrings')
library('devtools')
library('dataMeta')


#Read in fasta file, each line becomes an item in a vector
Control1=scan("control1.fasta",what=character(),sep="\n")

#Loop over vector of lines
for (Line in 1:length(Control1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(Control1[Line],">")){
    #Find the ORF in each sequence line and write to a new file
    Control1.txt = str_extract(Control1[Line],"ATG([ATGC]{3})+(TAG|TAA|TGA)?")
  }
}

#Make a new list called amino_acid_list and use a regex to put only nucleotide codons into it. Then use a regex to extract the amino acid names and use a for loop to combine the nucleotide sequences with the amino acid names.
amino_acid_list <- vector(mode="list", length=length(Codon_Map))
Codons <- str_extract(Codon_Map,"[A-Z]{3}")
names(amino_acid_list) <- Codons
amino_acids <- str_extract(Codon_Map,"[A-Za-z]{1,4}")
for (i in 1:length(amino_acid_list)){
  amino_acid_list[[i]] <- amino_acids[i]
}
  
### Stil need to apply amino_acid_list to the mRNA sequences to get the protein sequences
### Use a while loop (while i is less than or equal to the length of the reading frame, i+2)
New_Protein_Seq <- incorpoate_attr  

#Translate nucelotide sequence into its corresponding amino acid sequence.
translate(Control1.txt, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")


