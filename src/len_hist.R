source('lib.R')

###

NAMED <- 'DeepZ'
NAMEI<-'H3K27ac_A549.intersect_with_DeepZ'
#NAME <- 'H3K4me3_A549.ENCFF573MUH.hg19'
#NAME <- 'H3K4me3_A549.ENCFF573MUH.hg38'
#NAME <- 'H3K4me3_A549.ENCFF832EOL.hg19'
#NAME <- 'H3K4me3_A549.ENCFF832EOL.hg38'

###
NAME<-NAMEI
bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
#colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
colnames(bed_df) <- c('chrom', 'start', 'end')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.png'), path = OUT_DIR)
dev.off()