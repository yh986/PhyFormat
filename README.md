# PhyFormat
Helper functions in R to help with some phylogenetic data formatting.

## Installation
```
install.packages("devtools")
devtools::install_github("yh986/phyformat")
```
## Usage examples

```
x <- ReadFasta("sequences.fa") #Creates a physeq object from an existing FASTA file
x$Sequences[[1]] <- "PQHLIRV" #Edit a sequence in the physeq object
WriteFasta(x, "sequences2.fa") #Write to a fasta file
WritePhylip(x, "sequences2.phy") #'Convert'/write to a PHYLIP file
```

## Functions

`ReadFasta(filename)` and `ReadPhylip(filename)` are functions to read fasta or phylip files into a physeq object.


`WriteFasta(x, filename, multiline = 60)` and `WritePhylip(x, filename, blocks = 5, blocksize = 10)` are functions to write data in a list as described above, to a specified file, in either fasta or phylip format. Used in conjunction with ReadFasta or ReadPhylip, this may be used to convert files between the two formats. The parameters 'multiline' and 'blocksize' specify the number of characters to be used in each 'paragraph' in the respective functions.


`SubsetNewick(tipnames, newicktree)` will remove branches and nodes from the 'newicktree' (string in Newick format) such that all that is left of the tree is the branches named in the 'tipnames' array, and the internal nodes connecting them. When branches are merged (since useless internal nodes are removed), they are simply added together, which may not actually be phylogenetically accurate, depending on the distance model.
