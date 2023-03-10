#!/opt/R/4.0.5/bin/Rscript

# Experiment basic statistics
## Tom van Schaik, 171030

# - Functionality:
# Gather all statistics files into one table

# - Input
# (required)
# * Basename
# * Statistics documents
# * Output directory
#

# - Output
# * [basename].statistics.txt

# - Method
# Just read all the produced statistics files and combine into one table

# - Version
# 1.0 - Initial version

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(basename, output_dir, parsed, mapped, counts, 
                 mapped_strainspecific) {
  
  ## Set-up
  dir.create(output_dir, showWarnings = FALSE)
  output <- file.path(output_dir, paste0(basename, ".statistics.txt"))
  
  ## Gather all statistics
  
  # 1. Parsing
  parsed <- strsplit(parsed, ",")[[1]]
  
  reads.total <- 0
  reads.parsed <- 0
  
  for (p in parsed) {
    p <- read.table(p, sep = "\t", stringsAsFactors = FALSE, header = T)
    
    reads.total <- reads.total + p[, "reads"]
    reads.parsed <- reads.parsed + p[, "reads_written"]
  }
  
  # 2. Mapping 
  mapped <- strsplit(mapped, ",")[[1]]
  
  reads.mapped <- 0
  
  for (m in mapped) {
    m <- read.table(m, sep = "\t", stringsAsFactors = FALSE)
    m.mapped <- as.integer(sub(" .*", "", m[5, 1]))
    
    reads.mapped <- reads.mapped + m.mapped
  }
  
  # 3. Counts
  counts <- read.table(counts, sep = "\t")
  
  reads.counted <- counts[1, 2]
  reads.used <- counts[2, 2]
  
  # 4. Mapped - strain specific
  mapped_strainspecific <- strsplit(mapped_strainspecific, ",")[[1]]
  
  reads.mapped_ss <- 0
  reads.used_ss <- 0
  reads.129S1 <- 0
  reads.CAS <- 0
  reads.mismatch <- 0
  
  for (m in mapped_strainspecific) {
    lines <- readLines(m)

    TOTAL <- as.integer(gsub(".* ", "", lines[1]))
    REF <- as.integer(gsub(".* ", "", lines[2]))
    ALT <- as.integer(gsub(".* ", "", lines[3]))
    MIS <- as.integer(gsub(".* ", "", lines[4]))
    
    reads.mapped_ss <- reads.mapped_ss + TOTAL
    reads.used_ss <- reads.used_ss + REF + ALT
    reads.129S1 <- reads.129S1 + REF
    reads.CAS <- reads.CAS + ALT
    reads.mismatch <- reads.mismatch + MIS
    
  }
  
  ## Combine statistics and write 
  statistics <- data.frame(basename = basename,
                           total = reads.total,
                           parsed = reads.parsed,
                           mapped = reads.mapped,
                           counted = reads.counted,
                           used = reads.used,
                           strainspecific_total = reads.mapped_ss,
                           strainspecific_used = reads.used_ss,
                           strainspecific_129S1 = reads.129S1,
                           strainspecific_CAS = reads.CAS,
                           strainspecific_mismatch = reads.mismatch)
  
  # Write output table
  write.table(statistics, output, sep = "\t", row.names = F, quote = F)
  
}

##############################################################################
## Supporting functions ######################################################
##############################################################################

# None

##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.0"

option_list <- list(
  make_option(c("-b", "--basename"), type="character", default=NULL, 
              help="basename (required)", 
              metavar="character"),
  make_option(c("-p", "--parsed"), type="character", default=NULL, 
              help="parsed statistics, comma separated (required)", 
              metavar="character"),
  make_option(c("-m", "--mapped"), type="character", default=NULL, 
              help="mapped statistics, comma separated (required)", 
              metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="counts statistics, just one (required)", 
              metavar="character"),
  make_option(c("-s", "--mapped_strainspecific"), type="character", default=NULL, 
              help="mapped strain specific statistics (required)", 
              metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="output directory (required)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$basename) | is.null(opt$output_dir) | is.null(opt$parsed) |
    is.null(opt$mapped) | is.null(opt$counts) | is.null(opt$mapped_strainspecific)) {
  print_help(opt_parser)
  stop("Not all required arguments given\n", call.=FALSE)
} else {
  main(opt$basename, opt$output_dir, opt$parsed, opt$mapped, opt$counts,
       opt$mapped_strainspecific)
}