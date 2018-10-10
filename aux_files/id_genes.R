dir <- "genes_in_hg19"

refGene <- read.table(paste(, sep=""))

library(org.Hs.eg.db)
# For the reverse map:
x <- org.Hs.egREFSEQ2EG
# Get the RefSeq identifier that are mapped to an entrez gene ID
mapped_seqs <- mappedkeys(x)
# Convert to a list
xx <- as.data.frame(x[mapped_seqs])

new_refGene <- refGene
new_refGene$Entrez_ID <- NA
new_refGene$other <- NA
rows_refGene <-  unique(which(as.character(refGene[,1])%in%as.character(xx[,2])))
rows_ID <- unique(which(as.character(xx[,2])%in%as.character(refGene[,1])))

new_refGene[rows_refGene,c("Entrez_ID", "other")] <- xx[rows_ID,]

