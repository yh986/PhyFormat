# PhyFormat
Helper functions in R to help with some phylogenetic data formatting.

**FastaPhylip.R**


For formatting sequence data in the form of fasta or phylip.


ReadFasta(filename) and ReadPhylip(filename) are functions to read fasta or phylip files into an R list.
The list contains 2 chararacter arrays: "Headings", with the sequence names as indicated in the fasta or phylip files; and "Sequences".


WriteFasta(y, filename, multiline = 60) and WritePhylip(y, filename, blocks = 5, blocksize = 10) are functions to write data in a list as described above, to a specified file, in either fasta or phylip format. Used in conjunction with ReadFasta or ReadPhylip, this may be used to convert files between the two formats. The parameters 'multiline' and 'blocksize' specify the number of characters to be used in each 'paragraph' in the respective functions.


**SubsetNewick.R**

I couldn't find an R package that could easily subset a Newick-format phylogenetic tree! SubsetNewick(tipnames, newicktree) will remove branches and nodes from the 'newicktree' (string in Newick format) such that all that is left of the tree is the branches named in the 'tipnames' array, and the internal nodes connecting them. When branches are merged (since useless internal nodes are removed), they are simply added together, which may not actually be phylogenetically accurate, depending on the distance model.
