setwd("~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/")

#Read in fasta file, each line becomes an item in a vector
UniqueT=scan("uniquetranscripts.fasta",what=character(),sep="\n")


#Load package
library('stringr')
library('Biostrings')
library('devtools')
library('seqinr')


#Create a vector one half the length of Control1, Control2, Obese1, and Obese2.
UniqueT_ORF = character(length(UniqueT)/2)
Codon_Map=scan("codonmap.txt",what=character(),sep="\n")

###Determine the ORF of each sequence file using a for loop.
for (Line in 1:length(UniqueT)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(UniqueT[Line],">")){
    #Find the ORF in each sequence line and write to a new file
    #Search for the start codon (ATG) in each line of Control1 (with grepexpr) and assign it a new variable called atg_positions.
    atg_positions <- gregexpr('ATG', UniqueT[Line])[[1]]
    #Determine which codons are in the first reading frame (where atg_position is divisible by 3), and assign the position of the first to a new variable called start_codon_position.
    start_codon_position <- min(atg_positions[atg_positions%%3==1])
    #Extract sequences from each line of Control 1 from the first start codon position and place into a new variable called ORF_Temp.
    ORF_temp <- substring(UniqueT[Line], start_codon_position)
    #Repeat this process for the stop codon, by first searching for each of the three stop codons in ORF_Temp.
    end_positions <- gregexpr("(TAA)|(TAG)|(TGA)", ORF_temp)[[1]]
    #Determine which stop codons are in the first reading frame (where stop_codon_position is divisible by 3), and assign the position of the first to a new variable called end_positions.
    end_codon_position <- min(end_positions[end_positions%%3==1])
    #Extract the sequence from each line from the start codon in ORF_temp to the first stop codon and assign it a final varible.
    UniqueT_ORF[Line/2] = substring(ORF_temp, 1, end_codon_position+2)
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
### Create a new vector the length of *_ORF and use a for loop look at each sequnce. 
#Set a variable to one and while that variable is less then the number of characters in each sequence of *_.ORF, 
#replace each set of 3 nucleotides with the corresponding amino acid. Shfit the value of the nucleotide up 3 and repeat.
UniqueT_Protein = character(length(UniqueT_ORF))
for (sequence in 1:length(UniqueT_ORF)){
  nucleotide <- 1
  while (nucleotide < nchar(UniqueT_ORF[sequence])){
    UniqueT_Protein[sequence] = paste(UniqueT_Protein[sequence], amino_acid_list[substr(UniqueT_ORF[sequence], nucleotide,nucleotide+2)], sep="")
    nucleotide=nucleotide+3
  }
}

#Write each file out as a fasta file.
write.fasta(sequences = as.list(UniqueT_Protein), names = rep("UniqueT", times=length(UniqueT_Protein)),file.out = "UniqueT_Protein.fasta")





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



#Create a duplicate vector of Control1 called Control1_FASTA_Export.
Control1_FASTA_Export <- Control1
#Use a for loop to change the even enteries of Control1_FASTA_Export to the output amino acid sequences in Control1_Protein.
for (Line in seq(2, length(Control1_FASTA_Export), 2)){
  Control1_FASTA_Export[Line] <- Control1_Protein[Line/2] 
}


#Read in files
Codon_Map <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/codonmap.txt', what = , sep="\n", character())
Control1 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/control1.fasta', what = , sep="\n", character())
Control2 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/control2.fasta', what = , sep="\n", character())
Obese1 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/obese1.fasta', what = , sep="\n", character())
Obese2 <- scan('~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/obese2.fasta', what = , sep="\n", character())

