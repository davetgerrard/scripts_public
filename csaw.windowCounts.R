

library(csaw)
library(getopt)

### this file has been setup to use a config script specifying the input and output files. 

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'bam.file'  , 'b', 1, "character",
  'out.file' , 'o', 2, "character",
  'width' , 'w', 2, "integer"
), byrow=TRUE, ncol=4);
opt <- getopt(spec);
if ( is.null(opt$width ) ) { opt$width <- 10000}
if( is.null(opt$verbose)) {opt$verbose <- FALSE}
outFile <- paste(opt$bam.file,opt$width, "binCount.Rdata", sep=".") 
if( is.null(opt$out.file)) {opt$out.file <- outFile}

param <- readParam(minq=254)   # should catch all STAR uniquely mapped
print(param)
#bam.files <- c("STAR_hg38_FULL/STAR_hg38_FULL_H3K4me3_liver_AB_20140320/STAR_hg38_FULL_H3K4me3_liver_AB_20140320Aligned.sortedByCoord.out.bam",
#		"STAR_hg38_FULL/STAR_hg38_FULL_H3K4me3_brain_AB_20151208/STAR_hg38_FULL_H3K4me3_brain_AB_20151208Aligned.sortedByCoord.out.bam")
#data <- windowCounts(bam.files, ext=110, width=10, param=param)
data <- windowCounts(opt$bam.file, width=opt$width, bin=TRUE, param=param,filter=0)
save(data, file=opt$out.file)
