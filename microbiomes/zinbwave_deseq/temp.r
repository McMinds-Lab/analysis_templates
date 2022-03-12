load('ASVs.RData')

samplenames <- read.table('~/data/gainsbury_soil/Info on tags per sample.txt',sep='\t',header=T)

metadat <- read.table('~/data/gainsbury_soil/Soil eDNA_18S.txt',sep='\t',header=T)

id_conversion <- read.table('~/data/gainsbury_soil/argaly_id_conversion.txt',sep='\t',header=T)

rownames(seqtab) <- sapply(rownames(seqtab), function(x) samplenames$Sample[samplenames$Tag == sub('-',':',x)])

seqtab_F <- seqtab[grep('UFL',rownames(seqtab)),]

seqtab_M <- t(sapply(unique(sub('_.*','',rownames(seqtab_F))), function(x) apply(seqtab_F[grep(x,rownames(seqtab_F)),], 2, sum)))

m2 <- metadat[match(rownames(seqtab_M),metadat$id_argaly),]

library(ecodist)
library(vegan)

vstd <- DESeq2::varianceStabilizingTransformation(t(seqtab_M))

euc <- vegdist(t(vstd), method = "euclidean") 
pcoaVS <- pco(euc) 
plot(pcoaVS$vectors[,1], pcoaVS$vectors[,2], xlab = "PCA1", ylab = "PCA2",axes = TRUE, main = "PCA of normalized log relative abundance of ASVs", col=as.factor(m2$env.features), pch=16)
legend(x='topleft',legend=levels(as.factor(m2$env.features)),col=1:nlevels(as.factor(m2$env.features)),lty=1)

dna <- Biostrings::DNAStringSet(dada2::getSequences(seqtab))
load("~/data/ref/SILVA_SSU_r138_2019.RData")
ids <- DECIPHER::IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) 
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") 
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- dada2::getSequences(seqtab)
