#Input a phylogenetic tree in the Newick format.
#If a 'leaf' node has no siblings, then the parent of the leaf is a useless 'singleton' internal node.
#Then, merge the branches together to remove the node, until leaves have siblings.
  
RemoveSingletonTips <- function(Tree) {
  Iterativetree <- Tree
  while (TRUE) {
    #Find and merge one singleton
   
    singletoninfo <- unlist(
      strsplit(
        sub("^.*([(]([-A-Za-z0-9_<> ]+):([.0-9]+)[)]:([.0-9]+)).*$", 
            "\\2 \\3 \\4 \\1", Iterativetree)
        , 
          " ")
        #Matches ... (Xyz_0:3.43):1.235 ...
        #So that the 3.43 and 1.235 are then combined 
      )
      
      #Note: singletoninfo has the structure
      #3.43; 1.235; ...; Xyz_0
      
    if (length(singletoninfo) != 4)
    {break}
    
    #New branch length is a simple addition of the branch lengths before and after the internal node to be removed
    newbranchlength <- 
      as.numeric(singletoninfo[2]) + as.numeric(singletoninfo[3])
    newbranch <- paste0(singletoninfo[1], ":", newbranchlength)
    oldbranch <- singletoninfo[4]
    Iterativetree <- sub(oldbranch, newbranch, Iterativetree, fixed = T)
  }
  return (Iterativetree)
}

#Removes deep internal singletons, such that all subtrees have siblings.
RemoveSingletons <- function(Tree) {
  
  i <- 1
  subtrees <- NULL
  replacements <- NULL
  subtreepattern <- "^.*([(][-A-Za-z0-9_<> ]+:[.0-9]+,[-A-Za-z0-9_,.:<>]+[)]).*$"
  
  #First round of removing singleton tips, in case there are no subtrees at all
  Iterativetree <- RemoveSingletonTips(Tree)
  
  while(grepl(pattern = subtreepattern, x = Iterativetree)) {
    #Subtree substitution
    while(grepl(pattern = subtreepattern, x = Iterativetree)) {
      subtree_i <- sub(subtreepattern, "\\1", Iterativetree)
      replacement_i <- paste0("<", i, ">")
      subtrees <- c(subtrees, subtree_i)
      replacements <- c(replacements, replacement_i)
      Iterativetree <- sub(subtree_i, replacement_i, Iterativetree, fixed = TRUE)
      i <- i + 1
    }
    #Collapse the next layer
    Iterativetree <- RemoveSingletonTips(Iterativetree)
  }
  #End of the double loop of replacing subtrees and collapsing singletons.
  #Now, restore the subtrees back to its original form
  if (length(replacements) > 0) {
    for (i in length(replacements):1) {
      Iterativetree <- sub(pattern = replacements[i], 
                           replacement = subtrees[i],
                           x = Iterativetree)
    }
  }
  return(Iterativetree)
}

#Subset a newick tree by a specified set of tips;
#singleton branches are merged by simple addition (by RemoveSingletons).
#Leaf node names of the newicktree can only have
#Alphanumeric characters, spaces, dashes, or underscores.
#NO other punctuation marks! This will mess up the parsing.

#Bug: if there is just one header provided, it will fail to
#remove singletons (because of lack of commas)
SubsetNewick <- function(tipnames, newicktree) {
  tree_protected <- newicktree
  
  #1.
  #Get rid of any tips that are not specified in tipnames
  
  #Protect wanted tips
  for (i in 1:length(tipnames)) {
    header_i <- tipnames[[i]]
    protection_i <- paste0("<", i, ">")
    tree_protected <- 
      sub(pattern = header_i, replacement = protection_i, 
          x = tree_protected, fixed = T)
  }
  #Remove all non-protected tips
  tree_purged <- gsub("[-A-Za-z0-9_ ]+:[.0-9]+,?", "", tree_protected)
  #Remove hanging commas
  tree_purged <- gsub(",+)", ")", tree_purged)
  #Restore headers
  for (i in 1:length(tipnames)) {
    header_i <- tipnames[[i]]
    protection_i <- paste0("<", i, ">")
    tree_purged <- sub(pattern = protection_i, 
                       replacement = header_i, 
                       x = tree_purged,
                       fixed = T)
  }
  
  #2.
  #Remove exposed internal nodes
  uselessinternalnode <- "[(][)][:][.0-9]+,?"
  while(grepl(pattern = uselessinternalnode, x = tree_purged)) {
    #Removes "():0.123456," (regardless of presence of comma)
    tree_purged <- gsub(pattern = uselessinternalnode, 
                        replacement = "", 
                        x = tree_purged)
    #Remove hanging commas
    tree_purged <- gsub(pattern = ",+)", 
                        replacement = ")", 
                        x = tree_purged)
  }
  
  #3.
  #If a branch only has one child, then merge the branches together
  #to remove the "singleton" internal node.
  
  return(RemoveSingletons(tree_purged))
  
}


