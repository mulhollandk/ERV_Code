library(IRanges)
library(GenomicRanges)

ERV_Coords <- read.delim("hg19_rmsk_TE.gtf.txt", header=FALSE)
herv_coords <- ERV_Coords
herv_coords2 <- separate(data = herv_coords, col = V9, into = c("gene_id", "transcript_id", "family_id","class_id"), sep = "\\;")

Sig_DE_ERVs <- read.csv("LM_logfib_11covs_fdrsig.csv", header=TRUE, row.names=1)

Sig_ERV_Coords <- merge(Sig_DE_ERVs[,1:2], herv_coords2,by=0)
rownames(Sig_ERV_Coords) <- Sig_ERV_Coords[,1]
Sig_ERV_Coords <- Sig_ERV_Coords[,-1]

sig.herv.gr=GRanges(seqnames= Sig_ERV_Coords[,3],
ranges=IRanges(start= Sig_ERV_Coords[,6],
end= Sig_ERV_Coords[,7]),
				ERV_ID = Sig_ERV_Coords[,12],
				ERV_Class = Sig_ERV_Coords[,14],
				ERV_Family = Sig_ERV_Coords[,13],
				coeff= Sig_ERV_Coords[,1],
				pval= Sig_ERV_Coords[,2])

meth_coords <-read.delim("EPIC_hg19_CpG_sorted.bed",header=FALSE)
colnames(meth_coords) <- c("meth_chr","meth_start","meth_end","meth_cpg_name")

meth.gr=GRanges(seqnames= meth_coords[,1],
ranges=IRanges(start= meth_coords[,2],
end= meth_coords[,3]), CpG= meth_coords[,4])

sig.overlap <-mergeByOverlaps(sig.herv.gr,meth.gr)
sig.overlap <- as.data.frame(sig.overlap)
rownames(sig.overlap) <- make.names(sig.overlap$CpG, unique = TRUE)
overlap_sigherv_unqiue <- unique(sig.overlap[c("sig.herv.gr.ERV_ID")])
