

# ---TAG--- getopt
# ---TAG--- general
# ---TAG--- summarise


require(getopt)

#Col3:  0=noargument, 1=required argument, 2=optional argument
spec <- matrix(c(
  'dir'     , 'd', 2, "character", "search directory [default: getwd()]",
  'out'     , 'o', 1, "character", "outFile",
  'pattern'     , 'p', 1, "character", "searchPattern",
  'labCol'    , 'r',2, "integer"  , "index of row labels column [default: 1]",
  'dataCol'    , 'c',2, "integer"  , "index of data colum to extract [default: 2]",
  'labName'     , 'n', 2, "character", "row label column name [default: Statistic]", 
  'strip'     , 's', 1, "character", "pattern to strip from all sample names [default: '']",
  'verbose'    , 'v',0, "logical"  , "verbose output [default:FALSE]",
  'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$out) || is.null(opt$pattern) ) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# if not specified, provide defaults
if(is.null(opt$dir))  opt$dir <- getwd()
if(is.null(opt$verbose))  opt$verbose <- FALSE
if(is.null(opt$labName))  opt$labName <- "Statistic"
if(is.null(opt$labCol))  opt$labCol <- 1
if(is.null(opt$dataCol))  opt$dataCol <- 2
if(is.null(opt$strip))  opt$strip <- ""

outFilename <- opt$out
searchDir <- opt$dir
searchPattern <- opt$pattern

descriptionColumn <- opt$labCol
descriptionColumnName <- opt$labName
dataColumn <- opt$dataCol
strip.text <- opt$strip    # can be several items in a list (how to enter on command line?)

if(opt$verbose) { print(format(opt)) }

# how to input filenames?
# 1. give a directory and do a recursive search on a pattern?
# 2. give a list of files
# 3. give a file containing a list of files.

#searchDir <- "C:/Users/Dave/HanleyGroup/BerryRNAseq_multiTissue/data/starMappingHydra"

filesList <- list.files(path=searchDir, pattern= searchPattern,full.names = TRUE, recursive = TRUE)

#i <- 1
#descriptionColumn <- 1
#descriptionColumnName <- "Statistic"
#dataColumn <- 2
#strip.text <- "Log.final.out"    # can be several items in a list (how to enter on command line?)



#fileName <- "C:/Users/Dave/HanleyGroup/BerryRNAseq_multiTissue/data/starMappingHydra/STAR_hg38_G23_Thyroid_1Log.final.out"


for(i in 1:length(filesList)) {
  fileName <- filesList[i]
  if(opt$verbose) { print(fileName) }
  fileData <- read.delim(fileName, header=F)
  
  dim(fileData)
  head(fileData)
  
  sampleName <- basename(fileName)
  for(thisPattern in strip.text)  {   # implemented as multiple arguments but practically can only be single value using getopt package
    sampleName <- sub(thisPattern, "", sampleName )
  }
  
    
  if(i ==1) {
    combinedTable <- fileData[,c(descriptionColumn, dataColumn)]
    names(combinedTable) <- c(descriptionColumnName, sampleName)
  }  else {
    combinedTable[, sampleName] <- fileData[, dataColumn]
    
    
  }
}


write.table(combinedTable, file=outFilename, sep="\t", quote=F, row.names=F)


