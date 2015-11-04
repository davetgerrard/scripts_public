#!/usr/bin/env Rscript

# ---SUMMARY--- report the tags in scripts


# search a directory and report tags with counts
# search for scripts with a specific tag (or set of tags?)
# ---TAG--- getopt
# ---TAG--- general
# ---TAG--- summarise

# TODO implement searching for a specific tag
# TODO implementing searching for multiple tags. (possible with getopt?)

require(getopt)

#Col3:  0=noargument, 1=required argument, 2=optional argument
spec <- matrix(c(
  'dir'     , 'd', 2, "character", "search directory [default: getwd()]",
  'out'     , 'o', 1, "character", "outFile",
  'pattern'     , 'p', 1, "character", "searchPattern",
  'tags' , 't', 2, "character", "list of tags (separated by commas)",
  'verbose'    , 'v',0, "logical"  , "verbose output [default:FALSE]",
  'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help)  ) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

listMode <- FALSE
# if not specified, provide defaults
if(is.null(opt$dir))  opt$dir <- getwd()
if(is.null(opt$verbose))  opt$verbose <- FALSE
if(is.null(opt$tags)) listMode <- TRUE

outFilename <- opt$out
searchDir <- opt$dir

tagPattern <- paste0('---','TAG---')   # so that this line is picked up as a TAG!
remPattern <- paste0(".*", tagPattern)

file.list <- dir(path = searchDir)

if(listMode) {    # find all tags and report counts
count.vec <- integer()
for(thisFile in file.list) {
  hits <- grep(tagPattern, readLines(thisFile), value=TRUE)
  
  if(length(hits > 0))  {
    hits <- sub(remPattern," ", hits)
    hits <- gsub("^[[:space:]]*","", hits)  # remove leading whitespace
    hits <- gsub("[[:space:]]*$","", hits)  # remove trailing whitespace
    new.hits <- setdiff(hits, names(count.vec))
    count.vec[new.hits] <- 1
    existing.hits <- setdiff(hits, new.hits)
    count.vec[existing.hits] <- count.vec[existing.hits] + 1

    
  }
  }
  

count.vec <- sort(count.vec, decreasing=TRUE)

print(data.frame(TAG=names(count.vec), FILE_COUNT=count.vec), row.names=F)

} else {  # find specific tags and report file with these tags.  (either/both?)
  tag.vec <- unlist(strsplit(x=opt$tags, split=","))
  tag.vec <- gsub("^[[:space:]]*","", tag.vec)  # remove leading whitespace
  tag.vec <- gsub("[[:space:]]*$","", tag.vec)  # remove trailing whitespace
  print(paste("Searching for these tags:",opt$tags ))
  for(thisFile in file.list) {
    hits <- grep(tagPattern, readLines(thisFile), value=TRUE)
    
    if(length(hits > 0))  {
      hits <- sub(remPattern," ", hits)
      hits <- gsub("^[[:space:]]*","", hits)  # remove leading whitespace
      hits <- gsub("[[:space:]]*$","", hits)  # remove trailing whitespace
      if(length(setdiff(tag.vec, hits)) == 0) {
       print(thisFile) 
      }
      }
    }
  }


