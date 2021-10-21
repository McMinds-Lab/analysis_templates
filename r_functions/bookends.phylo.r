
bookends.phylo <- function(tree, target, depth_nodes=2, max_sisters=1, keep_trivial=FALSE) {

  inds <- which(tree$tip.label %in% target)
  keepers <- NULL
  for(current_tip in inds) {
    current_node <- current_tip
    putative_keepers <- NULL
    new_keepers <- NULL
    for(y in 1:depth_nodes) {
      putative_keepers <- phangorn::Descendants(tree,phangorn::Siblings(tree,current_node)[[1]])[[1]] ## find all tips that are sisters (y==1) or cousins (y>1) to the current tip
      if(keep_trivial | !any(putative_keepers %in% inds)) {
        putative_keepers <- putative_keepers[!putative_keepers %in% c(new_keepers,inds)] ## subset tips to those that have not already been selected for keeping
        new_keepers <- c(new_keepers, as.numeric(sample(as.character(putative_keepers), min(length(putative_keepers), max_sisters)))) ## sample max of two new tips and add to keepers list
      }
      current_node <- phangorn::Ancestors(tree,current_node,'parent')
    }
    keepers <- unique(c(keepers,new_keepers))
  }
  
  return(ape::keep.tip(tree,c(keepers,inds)))

}
