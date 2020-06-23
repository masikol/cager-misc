# |=== pub.R ===|

# Name stands for "Pick Up Barcodes".
# The script is designed for automatic picking up sequenceing barcodes.
# pub.R uses pam algorithm (from 'cluster package') to cluster barcodes and select ones that
#    are the most dissimilar from others. 


######################
# Parse CL arguments #
######################


# barcodes.fpath <- commandArgs(T)[1]
# num.samples <- commandArgs(T)[2]
# out.fpath <- commandArgs(T)[3]
# 
# if (! file.exists(barcodes.fpath)) {
#   stop(paste("File '", barcodes.fpath, "' does not exist!", sep = ""))
# }
# 
# if (! grepl("^[0-9]+$", num.samples)) {
#   stop(paste("Invalid number of samples:", num.samples))
# } else {
#   num.samples <- as.numeric(num.samples)
# }
# 
# if (is.na(out.fpath)) {
#   out.fpath <- "idx-pick-up_result"
# }

if (! require(cluster)) {
  cat("Package 'cluster' is not installed!\n")
  cat("You can install it with 'install.packages(\"cluster\")'\n")
  stop()
}

# Test setup

barcodes.fpath <- "/home/dell/Documents/barcodes.csv"
num.samples <- 3
out.fpath <- "/home/dell/Documents/test.out"

#########################################
# Define functions and data strunctures #
#########################################

# There cannot be other IUPAC characters in barcode sequences:
bases <- "AGCT"


# == Function 'check.structure' ==
# Function checks if input file is "single" or "double".
# Moreover, it checks if there is a header in the table.
# :param fpath: path to input file;
# :type fpath: character;
# Returns vector of logical values:
#    c(<IS_THERE_A_HEADER_IN_FILE>, <IS_FILE_IN_DOUBLE_FORMAT>)
check.structure <- function(fpath) {
  
  # Read the first line of file
  conn <- file(fpath)
  # Get rid of leading and trailing whitespaces
  line <- trimws(readLines(conn, n = 1), which = "both")
  close(conn)
  
  # Count commas (1 comma - "single" format, 3 commas -- "double" format)
  num.commas <- length( unlist(regmatches(line, gregexpr(',', line))) )
  if (num.commas != 1 && num.commas != 3) {
    stop("Unrecognized file structure. Please, see manual.")
  }
  
  # Prepare regexp patterns for format recognition
  single.patt <- paste("^.+,[", bases, "]+$", sep = "")
  double.patt <- paste("^.+,[", bases, "]+,.+,[", bases, "]+$", sep = "")
  
  # Detect format
  if ( grepl(single.patt, line, ignore.case = T) ) {
    return( c(FALSE, FALSE) ) # no header, "single format"
  } else if ( grepl(double.patt, line, ignore.case = T) ) {
    return( c(FALSE, TRUE) )  # no header, "double format"
  } else {
    double <- num.commas == 3 # TRUE if file in "double" format
    return( c(TRUE, double) ) # there is a header
  }
}

# == Function 'read.barcodes' ==
# Function reads a dataframe from input file.
# :param barcodes.fpath: path to input file;
# :type fpath: character;
# :param header: logical value indicating whether there is
#       a header in input file. TRUE if there is one;
# :type header: logical;
# :param double: logical value indicating whether input file
#       is in "double" format. TRUE if it is in "double" format;
# :type double: logical;
# Returns read dataframe with trimmed values.
read.barcodes <- function(barcodes.fpath, header, double) {
  
  # Define column names
  if (! double) {
    col.names <- c("Idx_name", "Idx_seq")
  } else {
    col.names <- c("Idx_name_F", "Idx_seq_F", "Idx_name_R", "Idx_seq_R")
  }
  
  # Read file and get dataframe
  barcodes.df <- read.csv(barcodes.fpath,
                   header = header,
                   col.names = col.names,
                   colClasses = "character",
                   comment.char = '#') # we want to skip commented lines
  
  # Get rid of whitespaces
  for (i in 1:length(barcodes.df)) {
    barcodes.df[[i]] <- trimws(barcodes.df[[i]], which = "both")
  }
  
  # Prepare regexp pattern for invalid barcodes recognition
  bases.patt <- paste("^[", bases, "]+$", sep = "")
  
  # Check if add barcodes are DNA sequences.
  # And bring them to upper case just in case.
  for (i in seq(2, length(barcodes.df), by = 2)) {
    
    barcodes.df[[i]] <- toupper(barcodes.df[[i]]) # to upper case
    
    if (! all( grepl(bases.patt, barcodes.df[[i]]) )) {
      stop("Your barcodes do not look line DNA sequences.")
    }
  }
  
  return(barcodes.df)
}

# == Quite obvious distance matrix ==
base.dists <- matrix(
  c(
#      A     C     G     T
     0.0,  1.0,  1.0,  1.0, # A
     1.0,  0.0,  1.0,  1.0, # C
     1.0,  1.0,  0.0,  1.0, # G
     1.0,  1.0,  1.0,  0.0  # T
  ),
  byrow = TRUE, nrow = nchar(bases), ncol = nchar(bases),
  dimnames = list(unlist(strsplit(bases, "")), unlist(strsplit(bases, "")))
)

# == Function 'get.base.dist' ==
# Function returns distance between two bases passed to it.
# :param base1: the first base to compare;
# :type base1: character;
# :param base2: the second base to compare;
# :type base2: logical;
# Returns distance between bases.
get.base.dist <- function(base1, base2) {
  return(base.dists[base1, base2])
}

# == Function 'seqdist' ==
# Function calculates distance between two barcode sequences.
# This distance if defined as sum of distances between pairwise bases in sequences.
# :param seq1: the first sequence to compare;
# :type seq1: character;
# :param seq2: the second sequence to compare;
# :type seq2: character;
# Returns distance between sequences.
seqdist <- function(seq1, seq2) {
  # Sequences will be cropped to min length
  min.len <- min(nchar(seq1), nchar(seq2))
  
  # Crop and split sequences (because R cannot iterate over strings)
  seq1 <- unlist(strsplit( substr(seq1, 1, min.len), "" ))
  seq2 <- unlist(strsplit( substr(seq2, 1, min.len), "" ))
  
  # Calculate distance
  dist <- sum( mapply(get.base.dist, seq1, seq2) )
  return(dist)
}

# == Function 'calc.dist.matr' ==
# Function calculates dissimilarity matrix for set of barcodes.
# :param barcodes.df: dataframe with names of barcodes and their sequences;
# :type barcodes.df: dataframe;
# Returns dissimilarity matrix.
calc.dist.matr <- function(barcodes.df) {
  
  # Dimensions
  n <- nrow(barcodes.df)
  n.sq <- n * n

  # Initialize with zeros
  dist.matr <- matrix(rep(0, n.sq),
                nrow = n, ncol = n,
                # Names of barcodes will be rownames and colnames
                dimnames = list(barcodes.df[[1]], barcodes.df[[1]]))
  
  # Shameful C-style loop for calculating distances.
  # Dissimilarity matrix is symmetric -- we won't recalculate already calculated values.
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      dist.matr[i, j] <- seqdist(barcodes.df$Idx_seq[i], barcodes.df$Idx_seq[j])
      dist.matr[j, i] <- dist.matr[i, j] # just copy -- matrix is symmertic
    }
  }
  
  return(dist.matr)
}


###########
# Proceed #
###########

# Check file struncture
header.double <- check.structure(barcodes.fpath)

# Get rid of vector -- we need just variables
header <- header.double[1]
double <- header.double[2]
rm(header.double) # remove!

# Print some message
if (double) {
  message("Input file is in \"double\" format.")
} else {
  message("Input file is in \"single\" format.")
}

# Read barcode dataframe
barcodes.df <- read.barcodes(barcodes.fpath, header, double)

# Check if there is enough barcodes in file
if (nrow(barcodes.df) < num.samples) {
  cat("\n")
  err.msg <- "Specified number of samples is greater than number of barcodes in file!\n"
  err.msg <- paste(err.msg, sprintf("Number of samples: %d;\n", num.samples),
                  sprintf("Number of barcodes: %d;\n", nrow(barcodes.df)))
  stop(err.msg)
}

# Calculate dissimilarity matrix
# 'barcodes.df[,1:2]' is passed to 'calc.dist.matr'
#   in order to have same syntax with "double"-formatted work mode.
message("Calculating distances...")
dist.matr <- calc.dist.matr(barcodes.df[,1:2])

# If input file is in "double" format -- calculate dissimilarity matrix
#   for reverse barcodes. Final matrix will be calculated as Euclidean distance
#   between the first and the second matrices.
if (double) {
  reverse.matr <- calc.dist.matr(barcodes.df[,3:4])
  dist.matr <- sqrt(dist.matr ^ 2 + reverse.matr ^ 2) # calculate final matrix
  rm(reverse.matr) # remove!
  
  # Rewrite rownames and colnames for double mode
  rownames(dist.matr) <- paste(barcodes.df$Idx_name_F,
                               barcodes.df$Idx_name_R, sep = ' -- ')
  colnames(dist.matr) <- paste(barcodes.df$Idx_name_F,
                               barcodes.df$Idx_name_R, sep = ' -- ')
}

# Cluster barcodes
message("Clustering...")
pam.clusters <- pam(dist.matr, num.samples, diss = T)

# print(summary(pam.clusters))

# Get representative barcodes
medoids <- pam.clusters$medoids

# Print the result
if (! double) {
  cat("\nBest indices:\n")
} else {
  cat("\nBest pairs of indices:\n")
}

for (barcode.name in medoids) {
  cat(sprintf("%s\n", barcode.name))
}



