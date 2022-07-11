genetree_paths <- list.files(path='outputs/primates/filtered_ensembl_trees', pattern='*.phy$', full.names=T)
names(genetree_paths) <- sapply(genetree_paths, function(x) paste(strsplit(basename(x),'.', fixed=T)[[1]][[1]], sep='_') )
genetrees <- lapply(genetree_paths, function(x) ape::read.tree(x))

species_strings <- c(Callithrix_jacchus='ENSCJAG',Homo_sapiens='ENSG',Microcebus_murinus='ENSMICG',Macaca_mulatta='ENSMMUG',Papio_hamadryas='ENSPANG',Pongo_abelii='ENSPPYG')

species_tree <- ape::read.tree('data/primate_allometry/Phylo tree files/pruned_elevenspecies.phy')
species_tree$tip.label[species_tree$tip.label == 'Pongo_pygmaeus'] <- 'Pongo_abelii'

species_tree_norm <- species_tree
species_tree_norm$edge.length <- species_tree_norm$edge.length / max(phytools::nodeHeights(species_tree_norm))

tol <- 1e-8
genetrees_chronos <- lapply(names(genetrees), \(tree) {

  cat(round(which(tree==names(genetrees)) / length(genetrees),2) * 100, '%: ', tree,'\n')

  tryCatch({
    genetree_species <- genetrees[[tree]]
    genetree_species$edge.length[genetree_species$edge.length < tol * max(genetree_species$edge.length)] <- tol * max(genetree_species$edge.length)
    genetree_species$tip.label <- sapply(genetrees[[tree]]$tip.label, \(x) names(species_strings)[startsWith(x,species_strings)])

    spec_nodes <- unique(genetree_species$edge[,1])[sapply(phangorn::Children(genetree_species, unique(genetree_species$edge[,1])),
                                                           \(x) length(intersect(genetree_species$tip.label[phangorn::Descendants(genetree_species, x[[1]])[[1]]],
                                                                                 genetree_species$tip.label[phangorn::Descendants(genetree_species, x[[2]])[[1]]])) == 0)]
    calib <- ape::makeChronosCalib(genetree_species,
                                   node = spec_nodes,
                                   age.min = sapply(phangorn::Descendants(genetree_species, spec_nodes),
                                                    \(x) 1 - phytools::findMRCA(species_tree_norm,
                                                                                unique(genetree_species$tip.label[x]),
                                                                                'height')))

    keep <- logical(nrow(calib))
    for(row in 1:nrow(calib)) {
      des <- phangorn::Descendants(genetree_species, calib$node[[row]], 'all')
      redund <- sapply(des, \(d) {
        if(d %in% calib$node) {
          calib$age.min[calib$node == d] >= calib$age.min[row]
        } else {
          FALSE
        }
      })
      keep[row] <- !any(redund)
    }


    gt1 <- ape::chronos(genetree_species,
                        lambda = 0.1,
                        calibration = calib[keep,])

    return(ape::di2multi(gt1))
  }, error = \(x) NA)

})

