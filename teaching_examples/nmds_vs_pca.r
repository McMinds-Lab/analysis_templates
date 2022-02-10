rdat <- matrix(rnorm(30*4),ncol=4)
rdat[1:10,1] <- rdat[1:10,1] + 25
rdat[11:20,2] <- rdat[11:20,2] + 25
rdat[21:30,3] <- rdat[21:30,3] + 25

pcs <- princomp(rdat)$scores

nmds <- vegan::monoMDS(vegan::vegdist(rdat, method = "euclidean"), pc=F)$points
