##Translating RNAseq files into aminoacid sequences
##
##Read the aa sequence file - codonmap.txt file
codonmap=scan(file = "codonmap.txt", what = character(), sep = "\n")
##
#Assign character vectors for file to extract open reading frame into them
match1=character(length(control1)/2)
match2=character(length(control2)/2)
match3=character(length(obese1)/2)
match4=character(length(obese2)/2)
##
##
#Read the nucleotide fasta files control 1&2 and obese 1&2
control1=scan(file = "control1.fasta.txt", what = character(), sep = "\n")
control2=scan(file = "control2.fasta.txt", what = character(), sep = "\n")
obese1=scan(file = "obese1.fasta.txt", what = character(), sep = "\n")
obese2=scan(file = "obese2.fasta.txt", what = character(), sep = "\n")
##
##
#Load stringr/i and seqinr packages
library(stringr)
library(stringi)
library(seqinr)
##
##
##
##
#forloop for control1 seq file orf extract
for (Line in 1:length(control1)){
#Operate only on lines that do not include >, skips header lines
#Note that str_detect returns a logical T/F
if (!str_detect(control1[Line],">")==TRUE){
#Find the ORF in each sequence line
match1[Line/2] = str_match(control1[Line],"([ATCG]{3})*?(ATG([ATCG]{3})+?(TAA|TAG|TGA))")[,3]
#Print the ORF to standard out
print(match1)
  }
}
##
#create a empty vector of length equal to match1
control1aa=character(length(match1))
control1protein=c(length(match1))
#use translate function to translate nucleotides into aa through a for loop
for(Line in 1:length(match1)){
  control1aa=(translate(s2c(match1[Line])))
  print(control1aa)
  control1protein[Line]=paste(control1aa, sep="", collapse="")
  print(control1protein)
  write.fasta(sequences = as.list(control1protein), names = rep("control1", times=length(control1protein)),file.out = "control1protein.fasta")
}
##
##
##Repeat the same for all the original files(control2,obese1,and obese2)
##
##
#forloop for control2 seq file orf extract
for (Line in 1:length(control2)){
#Operate only on lines that do not include >, skips header lines
#Note that str_detect returns a logical T/F
if (!str_detect(control2[Line],">")==TRUE){
#Find the ORF in each sequence line
match2[Line/2] = str_match(control2[Line],"([ATCG]{3})*?(ATG([ATCG]{3})+?(TAA|TAG|TGA))")[,3]
#Print the ORF to standard out
print(match2)
}
}
##
#create an empty vector of length equal to match2 each for storing translated aminoacid before and after formating
control2aa=character(length(match2))
control2protein=c(length(match2))
#use translate function to translate nucleotides into aa through a for loop
for(Line in 1:length(match2)){
  control2aa=(translate(s2c(match2[Line])))
  print(control2aa)
  control2protein[Line]=paste(control2aa, sep="", collapse="")
  print(control2protein)
  write.fasta(sequences = as.list(control2protein), names = rep("control2", times=length(control2protein)),file.out = "control2protein.fasta")
  
}
##
##
#forloop for obese1 seq file orf extract
for (Line in 1:length(obese1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(obese1[Line],">")==TRUE){
    #Find the ORF in each sequence line
    match3[Line/2] = str_match(obese1[Line],"([ATCG]{3})*?(ATG([ATCG]{3})+?(TAA|TAG|TGA))")[,3]
    #Print the ORF to standard out
    print(match3)
  }
}
##
##
#create an empty vector of length equal to match3 each for storing translated aminoacid before and after formating
obese1aa=character(length(match3))
obese1protein=c(length(match3))
#use translate function to translate nucleotides into aa through a for loop
for(Line in 1:length(match3)){
  obese1aa=(translate(s2c(match3[Line])))
  print(obese1aa)
  obese1protein[Line]=paste(obese1aa, sep="", collapse="")
  print(obese1protein)
  write.fasta(sequences = as.list(obese1protein), names = rep("obese1", times=length(obese1protein)),file.out = "obese1protein.fasta")
  
  }
##

##
#forloop for obese2 seq file orf extract
for (Line in 1:length(obese2)){
#Operate only on lines that do not include >, skips header lines
#Note that str_detect returns a logical T/F
if (!str_detect(obese2[Line],">")==TRUE){
#Find the ORF in each sequence line
match4[Line/2] = str_match(obese2[Line],"([ATCG]{3})*?(ATG([ATCG]{3})+?(TAA|TAG|TGA))")[,3]
    #Print the ORF to standard out
    print(match4)
}
  }
##
##
#create an empty vector of length equal to match4 each for storing translated aminoacid before and after formating
obese2aa=character(length(match4))
obese2protein=c(length(match4))
#use translate function to translate nucleotides into aa
for(Line in 1:length(match4)){
obese2aa=(translate(s2c(match4[Line])))
print(obese2aa)
obese2protein[Line]=paste(obese2aa, sep="", collapse="")
print(obese2protein)
write.fasta(sequences = as.list(obese2protein), names = rep("obese2", times=length(obese2protein)),file.out = "obese2protein.fasta")

  }

###This is all in UNIX on the command line.
###Align reference sequence files using muscle and then build an HMM profile using hmmbuild.

##forloop muscle alignment
for file in ./*.fasta
do
./muscle3.8.31_i86win32.exe -in $file -out $file.align
done
##
##
##for loop HMM profile Build
for file in ./*.align
do
./hmmbuild.exe $file.hmm $file
done



###Search all 4 translated “RNAseqfiles” for each of the 6 HMM protein models using hmmsearch.
for file in *protein.fasta
do 
for filename in *.hmm 
do 
./hmmsearch.exe --tblout $file.$filename.Results.txt $filename $file
done
done

#Cat all the results with a file name ending in Results.txt. Remove the lines beginning with #, and then print only columns 1 and 3.
#Use the unique funtion to count the number of hits for each file. Use awk one more time to remove space and append to a new file called All.Results.txt.
cat *Results.txt | grep -v '#' | awk '{print $1,$3}' | uniq -c | awk '{print $1,$2,$3}' >> All.Results.txt
##
###Plotting in R
#Load library
library(ggplot2)
##
##
#Read in fasta file, each line becomes an item in a vector
All_Results = read.table(file = "All.Results.txt", col.names = 1:3)
##
##
#Subsetting data by Transcript
Control1_Counts = subset(All_Results, X2 == "control1", c(1-3))
Control2_Counts = subset(All_Results, X2 == "control2", c(1-3))
Obese1_Counts = subset(All_Results, X2 == "obese1", c(1-3))
Obese2_Counts = subset(All_Results, X2 == "obese2", c(1-3))
##
##
#Plotting the expression levels
ggplot(data = Control1_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts")

ggplot(data = Control2_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "black", color = "black") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts")

ggplot(data = Obese1_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "blue", color = "blue") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts")

ggplot(data = Obese2_Counts)+geom_bar(aes(x =as.factor(X3), y = X1), stat = "summary",fun.y = "mean", fill = "blue", color = "blue") +
  theme_classic() +xlab("Transcript") +ylab("Transcript Counts")