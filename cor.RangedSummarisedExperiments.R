#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2016
#
######################################################


########FUNCTIONS


correlate.RSEs <- function(rse1, rse2, col.map=NULL, window.size=1*1000*1000, use1=TRUE, use2=TRUE, chroms=seqlevels(rse1), withinChrom=TRUE, 
                           cor.method="pearson")  {
  
  # take two RangedSummarisedExperiment objects and calculate correlations between their matched columns.
  # rse1
  # rse2
  # col.map
  # window.size
  # use1
  # use2
  # chroms
  # withinChrom - only compare regions on same chromosome (default TRUE, and also TRUE if window.size > 0)
    
  require(SummarizedExperiment)
  
  # check consistency of objects 
  commonCols <- intersect(row.names(colData(rse1)), row.names(colData(rse2)))
  n.common <- length(commonCols)
  print(paste(n.common, "common columns"))
  
  if(is.null(rownames(rse1)) ) rownames(rse1) <- 1:length(rse1)
  if(is.null(rownames(rse2)) ) rownames(rse2) <- 1:length(rse2)
  
  
  resultTable <- data.frame()
  # what size output is likely, do we need to add to a file?
  for(thisChrom in chroms) {

  rse1.valid <- rse1[seqnames(rse1) == thisChrom ]
  rse2.valid <- rse2[seqnames(rse2) == thisChrom ]
  print(paste(thisChrom, "Valid regions:", length(rse1.valid)))
  # check if anything left...
  # i <- 1
  for(i in 1:length(rse1.valid)) {
    # i <- i + 1
    cat(".")
    if(i %% 10 == 0) print(i)
    #focal.chr <- as.character(enhancers[i,"chr"])
   # focal.start <- enhancers[i,"start"]
    #focal.end <- enhancers[i,"end"]
   # focal.mid <- focal.start + floor(focal.end + 1 - focal.start)
   # en.values <- enhancers[i,samples.shared]
    
  # begin finding correlations
  
  # get index of rse2 elements within window.
  
    ov <- findOverlaps(rse1.valid[i], rse2.valid, maxgap = window.size)
    
    #subjectHits(ov)
    
    #assay(rse1[i, commonCols])   
    #apply(exp.subset, 1, FUN=function(x) cor(as.numeric(x[samples.shared]), as.numeric(en.values), method=cor.method))
    hitCors <- apply(assay(rse2.valid[subjectHits(ov), commonCols]) , 1, FUN=function(x) cor(as.numeric(x), as.numeric(assay(rse1.valid[i, commonCols]) ), method=cor.method))
    
    hitDists <- start(rse2.valid[subjectHits(ov)]) - start(rse1.valid[i])
    
    hitIds <- rownames(rse2.valid[subjectHits(ov)])
    
    
    iRow <- data.frame(i=rownames(rse1.valid[i]), hitIds = paste(hitIds, collapse=","),
                       hitCors = paste(signif(hitCors, digits=3), collapse=","),
                       hitDists = paste(hitDists, collapse=","))
  # 
  # paste(signif(exp.subset$en.scores, digits=3), collapse=",")
    resultTable <- rbind(resultTable,iRow)
  }
  
  }
  # return something
  return(resultTable)
}


###########################################


### OR use getopt:-

require(getopt)

#Col3:  0=noargument, 1=required argument, 2=optional argument
spec <- matrix(c(
    'rse1'     , 'a', 1, "character", "first rse rdata file",
     'rse2'     , 'b', 1, "character", "second rse rdata file",
    'output'     , 'o', 1, "character", "outFile",
    'chr' 	, 'c', 1, "character", "chromosome to score", 
    'verbose'    , 'v',0, "logical"  , "verbose output [default:FALSE]",
    'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$output)  ) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}


if(is.null(opt$verbose))  opt$verbose <- FALSE

outFilename <- opt$out


resultTable <- correlate.RSEs(rse1 =opt$rse1, rse2=opt$rse2, col.map=NULL, window.size=1*1000*1000, use1=TRUE, use2=TRUE, chroms=opt$chr, withinChrom=TRUE, cor.method="pearson")

print("finished correlations!")
 
save(resultTable, file=outFilename)

