#Read the aa sequence file - codonmap.txt file
codondata=scan(file = "codonmap.txt", what = character(), sep = "\n")

#Read the nucleotide fasta files control 1&2 and obese 1&2
control1=scan(file = "control1.fasta.txt", what = character(), sep = "\n")
control2=scan(file = "control2.fasta.txt", what = character(), sep = "\n")
obese1=scan(file = "obese1.fasta.txt", what = character(), sep = "\n")
obese2=scan(file = "obese2.fasta.txt", what = character(), sep = "\n")

#Translate nucleotide fasta files into aa fasta files

#Load package
library('stringr')
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

#Read in fasta file, each line becomes an item in a vector

#Loop over vector of lines
for (Line in 1:length(control1)){
  #Operate only on lines that do not include >, skips header lines
  #Note that str_detect returns a logical T/F
  if (!str_detect(control1[Line],">")){
    #Find the ORF in each sequence line
    match1 = str_extract(control1[Line],"ATG([ATCG]{3})+(TAA|TAG|TGA)")
    #Print the ORF to standard out
    print(match1)
  }
}

library(Biostrings)
library(BiocInstaller)
translate(match1, genetic.code = codondata, if.fuzzy.codon = "error")

for ( Line in 1:length(match1)){
  if (str_detect(match1[Line], "ATT")){
  control1aa = str_replace_all(match1, "ATT", "I")
print (control1aa)
  }
  if (str_detect(match1[Line], "ATC")){
    control1aa = str_replace_all(match1, "ATC", "I")
  }
  if (str_detect(match1[Line], "ATA")){
    control1aa = str_replace_all(match1, "ATA", "I")
  }
  if (str_detect(match1[Line], "CTT")){
    control1aa = str_replace_all(match1, "CTT", "L")
  }
  if (str_detect(match1[Line], "CTC")){
    control1aa = str_replace_all(match1, "CTC", "L")
  }
  if (str_detect(match1[Line], "CTA")){
    control1aa = str_replace_all(match1, "CTA", "L")
  }
    if (str_detect(match1[Line], "CTG")){
      control1aa = str_replace_all(match1, "CTG", "L")
    }
  if (str_detect(match1[Line], "TTA")){
    control1aa = str_replace_all(match1, "TTA", "L")
  }
  if (str_detect(match1[Line], "TTG")){
    control1aa = str_replace_all(match1, "TTG", "L")
  }
  if (str_detect(match1[Line], "GTT")){
    control1aa = str_replace_all(match1, "GTT", "V")
  }
  if (str_detect(match1[Line], "GTC")){
    control1aa = str_replace_all(match1, "GTC", "V")
  }
  if (str_detect(match1[Line], "GTA")){
    control1aa = str_replace_all(match1, "GTA", "V")
  }
  if (str_detect(match1[Line], "GTG")){
    control1aa = str_replace_all(match1, "GTG", "V")
  }
  if (str_detect(match1[Line], "TTT")){
    control1aa = str_replace_all(match1, "TTT", "F")
  }
  if (str_detect(match1[Line], "TTC")){
    control1aa = str_replace_all(match1, "TTC", "F")
  }
  if (str_detect(match1[Line], "ATG")){
    control1aa = str_replace_all(match1, "ATG", "M")
  }
  if (str_detect(match1[Line], "TGT")){
    control1aa = str_replace_all(match1, "TGT", "C")
  }
  if (str_detect(match1[Line], "TGC")){
    control1aa = str_replace_all(match1, "TGC", "C")
  }
  if (str_detect(match1[Line], "GCT")){
    control1aa = str_replace_all(match1, "GCT", "A")
  }
  if (str_detect(match1[Line], "GCC")){
    control1aa = str_replace_all(match1, "GCC", "A")
  }
  if (str_detect(match1[Line], "GCA")){
    control1aa = str_replace_all(match1, "GCA", "A")
  }
  if (str_detect(match1[Line], "GCG")){
    control1aa = str_replace_all(match1, "GCG", "A")
  }
  if (str_detect(match1[Line], "GGT")){
    control1aa = str_replace_all(match1, "GGT", "G")
  }
  if (str_detect(match1[Line], "GGC")){
    control1aa = str_replace_all(match1, "GGC", "G")
  }
  if (str_detect(match1[Line], "GGA")){
    control1aa = str_replace_all(match1, "GGA", "G")
  }
  if (str_detect(match1[Line], "GGG")){
    control1aa = str_replace_all(match1, "GGG", "G")
  }
  if (str_detect(match1[Line], "CCT")){
    control1aa = str_replace_all(match1, "CCT", "P")
  }
  if (str_detect(match1[Line], "CCA")){
    control1aa = str_replace_all(match1, "CCA", "P")
  }
  if (str_detect(match1[Line], "CCG")){
    control1aa = str_replace_all(match1, "CCG", "P")
  }
  if (str_detect(match1[Line], "ACT")){
    control1aa = str_replace_all(match1, "ACT", "T")
  }
  if (str_detect(match1[Line], "ACC")){
    control1aa = str_replace_all(match1, "ACC", "T")
  }
  if (str_detect(match1[Line], "ACA")){
    control1aa = str_replace_all(match1, "ACA", "T")
  }
  if (str_detect(match1[Line], "ACG")){
    control1aa = str_replace_all(match1, "ACG", "T")
  }
  if (str_detect(match1[Line], "TCT")){
    control1aa = str_replace_all(match1, "TCT", "S")
  }
  if (str_detect(match1[Line], "TCC")){
    control1aa = str_replace_all(match1, "TCC", "S")
  }
  if (str_detect(match1[Line], "TCA")){
    control1aa = str_replace_all(match1, "TCA", "S")
  }
  if (str_detect(match1[Line], "TCG")){
    control1aa = str_replace_all(match1, "TCG", "S")
  }
  if (str_detect(match1[Line], "AGT")){
    control1aa = str_replace_all(match1, "AGT", "S")
  }
  if (str_detect(match1[Line], "AGC")){
    control1aa = str_replace_all(match1, "AGC", "S")
  }
  if (str_detect(match1[Line], "TAT")){
    control1aa = str_replace_all(match1, "TAT", "Y")
  }
  if (str_detect(match1[Line], "TAC")){
    control1aa = str_replace_all(match1, "TAC", "Y")
  }
  if (str_detect(match1[Line], "TGG")){
    control1aa = str_replace_all(match1, "TGG", "W")
  }
  if (str_detect(match1[Line], "CAA")){
    control1aa = str_replace_all(match1, "CAA", "Q")
  }
  if (str_detect(match1[Line], "CAG")){
    control1aa = str_replace_all(match1, "CAG", "Q")
  }
  if (str_detect(match1[Line], "AAT")){
    control1aa = str_replace_all(match1, "AAT", "N")
  }
  if (str_detect(match1[Line], "AAC")){
    control1aa = str_replace_all(match1, "AAC", "N")
  }
  if (str_detect(match1[Line], "CAT")){
    control1aa = str_replace_all(match1, "CAT", "H")
  }
  if (str_detect(match1[Line], "CAC")){
    control1aa = str_replace_all(match1, "CAC", "H")
  }
  if (str_detect(match1[Line], "GAA")){
    control1aa = str_replace_all(match1, "GAA", "E")
  }
  if (str_detect(match1[Line], "GAG")){
    control1aa = str_replace_all(match1, "GAG", "E")
  }
  if (str_detect(match1[Line], "GAT")){
    control1aa = str_replace_all(match1, "GAT", "D")
  }
  if (str_detect(match1[Line], "GAC")){
    control1aa = str_replace_all(match1, "GAC", "D")
  }
  if (str_detect(match1[Line], "GAC")){
    control1aa = str_replace_all(match1, "GAC", "D")
  }
  if (str_detect(match1[Line], "GAC")){
    control1aa = str_replace_all(match1, "GAC", "D")
  }
  }




D	GAC
K	AAA
K	AAG

print (control1aa)















for (i in 1:length(match1)){
  if (str_detect(match1[i],'ATG')==TRUE){
    control1aa[1]=str_replace_all(match1[i],"ATT","I")
  }
  
  
  
  grep("I", control1aa)
  

grep("ATG", match1)

"sed -i -e 's/ATG/I' match1"
for (i in 1:length(vari)){
  if (str_detect(vari[i],'CHROM')==TRUE){
    vari[2]=str_replace_all(vari[i],texasnames,"Cf.Sfa.")
  }

translate()




for ( Line in 1:length(match1)){
  if (str_detect(match1,  "ATG([ATGC]{3})+(TAA|TAG|TGA)")){
    control1aa = str_replace_all(string, pattern, replacement)
    (data[>], "ATG([ATGC]{3})+(TAA|TAG|TGA)")
    print match
  }
}

# grep("ATG([ATGC]{3})+(TAA|TAG|TGA)", control1)
# str_extract_all("control1", "ATG([ATGC]{3})+(TAA|TAG|TGA)")


#/location="(.*?)"/



library(stringi)
library(stringr)
str_replace(control1, ATT, "I")

#df$aa = translate(readDNAStringSet(control1))

#Reading the nucleotide fasta into aa fasta - start with startcodon and end at first stopcodon - use sed? (manual?)