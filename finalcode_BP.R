#Read the aa sequence file - codonmap.txt file
codondata=scan(file = "codonmap.txt", what = character(), sep = "\n")
codondict=read.table("codonmap.txt", header=FALSE, sep="")

#Assign codons and aminoacid values from codon dictionary
codons=codondict[,2]
aa=codondict[,1]
match1=character(length(294))

#Read the nucleotide fasta files control 1&2 and obese 1&2
control1=scan(file = "control1.fasta.txt", what = character(), sep = "\n")
control2=scan(file = "control2.fasta.txt", what = character(), sep = "\n")
obese1=scan(file = "obese1.fasta.txt", what = character(), sep = "\n")
obese2=scan(file = "obese2.fasta.txt", what = character(), sep = "\n")

#Load stringr/i and seqinr packages
library(stringr)
library(stringi)
library(seqinr)

#Assign regex for openreading frame (ORF)
ORF="ATG([ATCG]{3})+(TAA|TAG|TGA)?"

#forloop for control1 seq file orf extract
for (Line in 1:length(control1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(control1[Line],">")==TRUE){
    #Find the ORF in each sequence line
    match1 = str_extract(control1[Line],"ATG([ATCG]{3})+(TAA|TAG|TGA)?")
    #Print the ORF to standard out
    print(match1)
  }
}

#use translate function to translate nucleotides into aa
control1aa=(translate(s2c(match1)))
#convert vector of aa into a dataframe      
seqc1<-as.data.frame(control1aa)

