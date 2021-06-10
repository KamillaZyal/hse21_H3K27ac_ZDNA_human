source('lib.R')

###

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
 BiocManager::install("ChIPseeker")
 BiocManager::install("clusterProfiler")
 BiocManager::install("org.Hs.eg.db")

 source("http://bioconductor.org/biocLite.R")
 ## biocLite("BiocUpgrade") ## you may need this
biocLite("ChIPseeker") 
 
 # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
 install.packages("libcurl4-openssl-dev") 
 library(ChIPseeker)
 library(TxDb.Hsapiens.UCSC.hg19.knownGene)
 #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
 library(clusterProfiler)
 library(org.Hs.eg.db)
###

#NAME <- 'H3K4me3_A549.intersect_with_DeepZ'
#NAME <- 'DeepZ'
#NAME <- 'H3K4me3_A549.ENCFF573MUH.hg19.filtered'
NAME191 <- 'H3K27ac_A549.ENCFF389RXK.hg19.filtered'
NAME192 <- 'H3K27ac_A549.ENCFF926NKP.hg19.filtered' 
NAMEZDNA <- 'DeepZ'
NAMEZDI <- 'H3K27ac_A549.intersect_with_DeepZ'
NAME <- NAMEZDI
BED_FN <- paste0(DATA_DIR, NAME, '.bed')

###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(BED_FN, tssRegion=c(-2000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.pdf'))
png(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.png'))
plotAnnoPie(peakAnno)
dev.off()

# peak <- readPeakFile(BED_FN)
# pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.covplot.pdf'))
# covplot(peak, weightCol="V5")
# dev.off()
# 