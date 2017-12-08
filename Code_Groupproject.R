#Pseudocode to translate nucleotide sequences to aminoacid sequences

#Read the aa sequence file - codonmap.txt file
#Read the nucleotide fasta files control 1&2 and obsese 1&2
#Translate nucleotide fasta files into aa fasta files
#Reading the nucleotide fasta into aa fasta - start with startcodon and end at first stopcodon - use sed?(manual?)

#UNIX  code to generate Top Hits
for file in *.csv
do
cat $file | head -1 $file >> TopHits2.csv 
done


#Set working directory and read each of the original neucleotide fasta files.
setwd("~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Original_Files/")

#Read in fasta file, each line becomes an item in a vector
Codon_Map=scan("codonmap.txt",what=character(),sep="\n")
Control1=scan("control1.fasta",what=character(),sep="\n")
Control2=scan("control2.fasta",what=character(),sep="\n")
Obese1=scan("obese1.fasta",what=character(),sep="\n")
Obese2=scan("obese2.fasta",what=character(),sep="\n")

#Load package
library('stringr')
library('Biostrings')
library('devtools')
library('seqinr')


#Create a vector one half the length of Control1, Control2, Obese1, and Obese2.
Control1_ORF = character(length(Control1)/2)
Control2_ORF = character(length(Control2)/2)
Obese1_ORF = character(length(Obese2)/2)
Obese2_ORF = character(length(Obese2)/2)

###Determine the ORF of each sequence file using a for loop.
for (Line in 1:length(Control1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(Control1[Line],">")){
    #Find the ORF in each sequence line and write to a new file
    #Search for the start codon (ATG) in each line of Control1 (with grepexpr) and assign it a new variable called atg_positions.
    atg_positions <- gregexpr('ATG', Control1[Line])[[1]]
    #Determine which codons are in the first reading frame (where atg_position is divisible by 3), and assign the position of the first to a new variable called start_codon_position.
    start_codon_position <- min(atg_positions[atg_positions%%3==1])
    #Extract sequences from each line of Control 1 from the first start codon position and place into a new variable called ORF_Temp.
    ORF_temp <- substring(Control1[Line], start_codon_position)
    #Repeat this process for the stop codon, by first searching for each of the three stop codons in ORF_Temp.
    end_positions <- gregexpr("(TAA)|(TAG)|(TGA)", ORF_temp)[[1]]
    #Determine which stop codons are in the first reading frame (where stop_codon_position is divisible by 3), and assign the position of the first to a new variable called end_positions.
    end_codon_position <- min(end_positions[end_positions%%3==1])
    #Extract the sequence from each line from the start codon in ORF_temp to the first stop codon and assign it a final varible.
    Control1_ORF[Line/2] = substring(ORF_temp, 1, end_codon_position+2)
    #Control1_ORF[Line/2] = str_extract(ORF_temp,"ATG([ATGC]{3})+((TAG)|(TAA)|(TGA))?")
  }
}


#Repeat for each file: Control2, Obese1, and Obese2
#Control2
for (Line in 1:length(Control2)){
  if (!str_detect(Control2[Line],">")){
    atg_positions2 <- gregexpr('ATG', Control2[Line])[[1]]
    start_codon_position2 <- min(atg_positions2[atg_positions2%%3==1])
    ORF_temp2 <- substring(Control2[Line], start_codon_position2)
    end_positions2 <- gregexpr("(TAA)|(TAG)|(TGA)", ORF_temp2)[[1]]
    end_codon_position2 <- min(end_positions2[end_positions2%%3==1])
    Control2_ORF[Line/2] = substring(ORF_temp2, 1, end_codon_position2+2)
  }
}

#Obese1
for (Line in 1:length(Obese1)){
  if (!str_detect(Obese1[Line],">")){
    atg_positions3 <- gregexpr('ATG', Obese1[Line])[[1]]
    start_codon_position3 <- min(atg_positions3[atg_positions3%%3==1])
    ORF_temp3 <- substring(Obese1[Line], start_codon_position3)
    end_positions3 <- gregexpr("(TAA)|(TAG)|(TGA)", ORF_temp3)[[1]]
    end_codon_position3 <- min(end_positions3[end_positions3%%3==1])
    Obese1_ORF[Line/2] = substring(ORF_temp3, 1, end_codon_position3+2)
  }
}

#Obese2
for (Line in 1:length(Obese2)){
  if (!str_detect(Obese2[Line],">")){
    atg_positions4 <- gregexpr('ATG', Obese2[Line])[[1]]
    start_codon_position4 <- min(atg_positions4[atg_positions4%%3==1])
    ORF_temp4 <- substring(Obese2[Line], start_codon_position4)
    end_positions4 <- gregexpr("(TAA)|(TAG)|(TGA)", ORF_temp4)[[1]]
    end_codon_position4 <- min(end_positions4[end_positions4%%3==1])
    Obese2_ORF[Line/2] = substring(ORF_temp4, 1, end_codon_position4+2)
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

### Create a new vector the length of *_ORF and use a for loop look at each sequnce. 
    #Set a variable to one and while that variable is less then the number of characters in each sequence of *_.ORF, 
    #replace each set of 3 nucleotides with the corresponding amino acid. Shfit the value of the nucleotide up 3 and repeat.
Control1_Protein = character(length(Control1_ORF))
for (sequence in 1:length(Control1_ORF)){
  nucleotide <- 1
    while (nucleotide < nchar(Control1_ORF[sequence])){
    Control1_Protein[sequence] = paste(Control1_Protein[sequence], amino_acid_list[substr(Control1_ORF[sequence], nucleotide,nucleotide+2)], sep="")
    nucleotide=nucleotide+3
  }
}

#Repeat for Control2, Obese1, and Obese2.
#Control2
Control2_Protein = character(length(Control2_ORF))
for (sequence in 1:length(Control2_ORF)){
  nucleotide <- 1
  while (nucleotide < nchar(Control2_ORF[sequence])){
    Control2_Protein[sequence] = paste(Control2_Protein[sequence], amino_acid_list[substr(Control2_ORF[sequence], nucleotide,nucleotide+2)], sep="")
    nucleotide=nucleotide+3
  }
}

#Obese1
Obese1_Protein = character(length(Obese1_ORF))
for (sequence in 1:length(Obese1_ORF)){
  nucleotide <- 1
  while (nucleotide < nchar(Obese1_ORF[sequence])){
    Obese1_Protein[sequence] = paste(Obese1_Protein[sequence], amino_acid_list[substr(Obese1_ORF[sequence], nucleotide,nucleotide+2)], sep="")
    nucleotide=nucleotide+3
  }
}

#Obese2
Obese2_Protein = character(length(Obese2_ORF))
for (sequence in 1:length(Obese2_ORF)){
  nucleotide <- 1
  while (nucleotide < nchar(Obese2_ORF[sequence])){
    Obese2_Protein[sequence] = paste(Obese2_Protein[sequence], amino_acid_list[substr(Obese2_ORF[sequence], nucleotide,nucleotide+2)], sep="")
    nucleotide=nucleotide+3
  }
}


#Write each file out as a fasta file.
write.fasta(sequences = as.list(Control1_Protein), names = rep("Control1", times=length(Control1_Protein)),file.out = "Control1_Protein.fasta")
write.fasta(sequences = as.list(Control2_Protein), names = rep("Control2", times=length(Control2_Protein)),file.out = "Control2_Protein.fasta")
write.fasta(sequences = as.list(Obese1_Protein), names = rep("Obese1", times=length(Obese1_Protein)),file.out = "Obese1_Protein.fasta")
write.fasta(sequences = as.list(Obese2_Protein), names = rep("Obese2", times=length(Obese2_Protein)),file.out = "Obese2_Protein.fasta")




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

###Search all 4 translated “RNAseqfiles” for each of the 6 HMM protein models using hmmsearch.
for file in *_Protein.fasta
do 
for filename in *.hmm
do 
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmsearch --tblout $file.$filename.Results.txt $filename $file 
done 
done


#Cat all the results with a file name ending in Results.txt. Remove the lines beginning with #, and then print only columns 1 and 3.
#Use the unique funtion to count the number of hits for each file. Use awk one more time to remove space and append to a new file called All.Results.txt.
cat *Results.txt | grep -v '#' | awk '{print $1,$3}' | uniq -c | awk '{print $1,$2,$3}' >> All.Results.txt



##Back in R
setwd("~/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Protein_FASTA_Files/")

#Load library
library(ggplot2)

#Read in fasta file, each line becomes an item in a vector
All_Results = read.table(file = "All.Results.txt", col.names = 1:3)

#Subset data by Transcript
Control1_Counts = subset(All_Results, X2 == "Control1", c(1-3))
Control2_Counts = subset(All_Results, X2 == "Control2", c(1-3))
Obese1_Counts = subset(All_Results, X2 == "Obese1", c(1-3))
Obese2_Counts = subset(All_Results, X2 == "Obese2", c(1-3))

#Plot
ggplot(data = Control1_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = Control2_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = Obese1_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = Obese2_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))




#Further Exploration 2
###This is all in UNIX on the command line.
###Align reference sequence files using muscle and then build an HMM profile using hmmbuild.
#Set working directory
cd /Users/lou/Desktop/data-shell/BioInfo_GroupProject/Bioinformatics_GroupProject_BLS/Further_Exploration_2/

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

###Search all 4 translated “RNAseqfiles” for each of the 6 HMM protein models using hmmsearch.
for file in *_Protein.fasta
do 
for filename in *.hmm
do 
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmsearch --tblout $file.$filename.Results.txt $filename $file 
done 
done

#Cat all the results with a file name ending in Results.txt. Remove the lines beginning with #, and then print only columns 1 and 3.
#Use the unique funtion to count the number of hits for each file. Use awk one more time to remove space and append to a new file called All.Results.txt.
cat *Results.txt | grep -v '#' | awk '{print $1,$3}' | uniq -c | awk '{print $1,$2,$3}' >> All.Results.txt

