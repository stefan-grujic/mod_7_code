
#Converts 16S_XString.fasta to a newick format tree.

rootloc <- readline(prompt="Enter Project Directory: ")
treeloc <- paste(rootloc,'/TREE', sep='')
scoloc <- paste(rootloc,'/SCOARY/INPUT/tree.newick', sep='')

library(pacman)
p_load('ape', 'phangorn', 'msa')

seqs <- readDNAStringSet(paste(treeloc, '/16S_XString.fasta',sep=''))
aln <- msaClustalOmega(seqs)

aln_ape <- msaConvert(aln,'ape::DNAbin')

distdna<-dist.dna(aln_ape,model ="N",pairwise.deletion = TRUE,as.matrix = TRUE)
af_dist_transformed<-dist(distdna)
tre<-nj(af_dist_transformed)

ws<-midpoint(tre)
ws$edge.length<-ifelse(ws$edge.length>0,ws$edge.length,-ws$edge.length)
write.tree(ws, file = scoloc)
write.tree(ws, file = paste(treeloc,"/tree.newick",sep=''))
