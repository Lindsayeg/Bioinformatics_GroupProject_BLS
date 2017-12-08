#Pseudocode to translate nucleotide sequences to aminoacid sequences

#Read the aa sequence file - codonmap.txt file
#Read the nucleotide fasta files control 1&2 and obsese 1&2
#Translate nucleotide fasta files into aa fasta files
#Reading the nucleotide fasta into aa fasta - start with startcodon and end at first stopcodon - use sed?(manual?)

#Set working directory and read each of the original neucleotide fasta files.
1
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
for file in ./../Original_Files/*Protein.fasta
do
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmsearch --tblout Results.txt $file.hmm $file
cat Results.txt >> All_Results.txt
done

#Moved the generated protein fasta files into the smaee folder as the  
for file in *_Protein.fasta
do
for filename in *.hmm
do
/Users/lou/Desktop/hmmer-3.1b2-macosx-intel/binaries/hmmsearch --tblout $file.$filename.Results.txt $filename $file
done
done


cat *Results.txt | grep -v '#' | awk '{print $1,$3}' | uniq -c | awk '{print $1,$2,$3}' 

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

