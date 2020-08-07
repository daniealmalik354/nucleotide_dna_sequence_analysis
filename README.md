---
Title: "nucleotide_dna_sequence_analysis" 
Author: "daniealmalik"
---
## I will try to download nucleotide sequence from GenBank such as NCBI, and analysis through R

```{r loading_packages, message = FALSE}
#Packages 
library(ape) 
library(seqinr) 
library(dplyr) 
library(ggtree) 
library(ggplot2)

#list of accesion number of sequence files from genebank that we want to download 
list <- c("KP298953", "KC313384", "AY640303", "KC313382", "AY640302")

#download the data sequences 
seq_list <- read.GenBank(list); seq_list

#explore the information of sequence file 
str(seq_list) 
attributes(seq_list) 
names(seq_list) 
attr(seq_list, "species")

#we want to set the atribute information sequence data is just name and accesion number 
seq_list_ID <- paste(attr(seq_list, "species"), names(seq_list)); seq_list_ID

#write sequence files as fasta format 
write.dna(seq_list, file = "seq_list.fasta", format = "fasta")

#read a fasta file using the seqinr package 
seq_list_fasta <- read.fasta(file = "seq_list.fasta", seqtype = "DNA", as.string = TRUE) 

#write the fasta file with names format 
write.fasta(sequences = seq_list_fasta, names = seq_list_ID, file.out = "seq_list_fasta_format.fasta")

#Align sequence files using 
library(msa) 
library(DECIPHER) 
seq_list_fasta_format <- readDNAStringSet("seq_list_fasta_format.fasta"); seq_list_fasta_format 

seq_list_align <- AlignSeqs(seq_list_fasta_format); seq_list_align

writeXStringSet(seq_list_align, filepath = "seq_list_align.fasta", format = "fasta")

#Make a tree with Neighbor Joining method 
data_align <- read.alignment(file = "seq_list_align.fasta", format = "fasta")

data_align <- read.dna(file = "seq_list_align.fasta", format = "fasta") class(data_align)

data_align_distace <- dist.dna(data_align, model = "TN93")

tree <- nj(data_align_distace)

tree2 <- root(tree, out = 1) boots <- boot.phylo(tree2, data_align, function(e) root(nj(dist.dna(e, model = "TN93")), 1))

write.tree(tree2, "tree.tre")

plot(tree2, edge.width = 5) 
axisPhylo() 
nodelabels(boots, cex = 1, frame = "none", adj = c(1.3, 1.5))
```
