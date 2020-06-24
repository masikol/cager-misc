# |=== pub.R ===|

# Name stands for "Pick Up Barcodes".
# The script is designed for automatic picking up sequenceing barcodes.
# pub.R uses pam algorithm (from 'cluster' package) to cluster barcodes and select ones that
#    are the most dissimilar from others.

script.version <- "1.0.a"
# Last modified 2020.06.24

######################
# Parse CL arguments #
######################

usage.msg <- "pub.R\nUsage:\n  Rscript pub.R <csv_file_with_barcodes> <number_of_samples>"

# Print help message
if (commandArgs(T)[1] %in% c("-h", "-help", "--help")) {
  cat(paste0("pub.R; Version ", script.version, '\n'))
  cat(paste0(usage.msg, '\n'))
  cat("Example:\n")
  cat("Pick up 28 barcodes from file `my_favorite_barcodes.csv`.\n")
  cat("\n  Rscript pub.R my_favorite_barcodes.csv 28\n")
  
  quit(save = "no")
}

if (length(commandArgs(T)) != 2) {
  cat(paste0(usage.msg, '\n'))
  quit(save = "no", status = 1)
}

barcodes.fpath <- commandArgs(T)[1]
num.samples <- commandArgs(T)[2]

if (! file.exists(barcodes.fpath)) {
  stop(paste0("File '", barcodes.fpath, "' does not exist!"))
}

if (! grepl("^[0-9]+$", num.samples)) {
  stop(paste0("Invalid number of samples: '", num.samples, "'"))
} else {
  num.samples <- as.numeric(num.samples)
}

cat(paste0("pub.R; Version ", script.version, '\n\n'))

if (! require(cluster)) {
  cat("Package 'cluster' is not installed!\n")
  cat("You can install it with 'install.packages(\"cluster\")'\n")
  stop()
}


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
  single.patt <- paste0("^.+,[", bases, "]+$")
  double.patt <- paste0("^.+,[", bases, "]+,.+,[", bases, "]+$")
  
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
  bases.patt <- paste0("^[", bases, "]+$")
  
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

# == Function 'calc.diss.matr' ==
# Function calculates dissimilarity matrix for set of barcodes.
# :param barcodes.df: dataframe with names of barcodes and their sequences;
# :type barcodes.df: dataframe;
# Returns dissimilarity matrix.
calc.diss.matr <- function(barcodes.df) {
  
  # Dimensions
  n <- nrow(barcodes.df)
  n.sq <- n * n

  # Initialize with zeros
  diss.matr <- matrix(rep(0, n.sq),
                nrow = n, ncol = n,
                # Names of barcodes will be rownames and colnames
                dimnames = list(barcodes.df[[1]], barcodes.df[[1]]))
  
  # Shameful C-style loop for calculating distances.
  # Dissimilarity matrix is symmetric -- we won't recalculate already calculated values.
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      diss.matr[i, j] <- seqdist(barcodes.df$Idx_seq[i], barcodes.df$Idx_seq[j])
      diss.matr[j, i] <- diss.matr[i, j] # just copy -- matrix is symmertic
    }
  }
  
  return(diss.matr)
}

# |=== Functions for checking color ballance ===|

red.group  <- c('C', 'A') # "red" bases
green.group <- c('G', 'T') # "green" bases

# == Function is.ballanced ==
# Function checks if set of sequences are color-ballanced,
# See this page for details: https://international.neb.com/protocols/2016/03/23/index-pooling-guidelines-neb-e6609
# :param seqs: vector of sequences to check for color ballance;
# :type seqs: vector<character>;
# Returns TRUE if set of sequences is color-ballanced, else FALSE.
is.ballanced <- function(seqs) {
  
  # Find min length among all sequences
  min.len <- min(as.vector(sapply(X = seqs, FUN = nchar)))
  
  for (i in 1:min.len) {
    # Select i-th characters of all sequences
    chars <- as.vector(sapply(X = seqs, FUN = substr,
                       start = i, stop = i)) # arguments for substr
    
    # Check if there are at least one nucleotide from each group
    red   <- any(red.group %in% chars)
    green <- any(green.group %in% chars)

    if (! (red && green)) {
      return(FALSE) # color ballance is corrupted
    }
  }
  
  return(TRUE)
}


# == Function get.forw.names ==
# Returns returns name of forward barcode (first in pair)
#   given my ad hoc naming scheme: "<forw_name> -- <rev_name>"
# :param combined.name: name of "double" barcodes divided by ' -- ';
# :type combined.name: character;
# Returns returns name of forward barcode.
get.forw.names <- function(combined.name) {
  return(strsplit(combined.name, " -- ")[[1]][1])
}


# == Function get.rev.names ==
# Returns returns name of reverce barcode (second in pair)
#   given my ad hoc naming scheme: "<forw_name> -- <rev_name>"
# :param combined.name: name of "double" barcodes divided by ' -- ';
# :type combined.name: character;
# Returns returns name of reverce barcode.
get.rev.names <- function(combined.name) {
  return(strsplit(combined.name, " -- ")[[1]][2])
}


# == Function get.seqs ==
# Function finds sequences of barcodes in 'barcodes.df'
#   given their names.
# :param barcode.names: vector of barcode names;
# :type barcode.names: vector<character>;
# :param barcodes.df: dataframe with names of barcodes and their sequences;
# :type barcodes.df: dataframe;
# :param name.i: integer indicating whether it is forward or reverse barcode.
#    1 for forward (and "single"), 3 for reverse;
# :type name.i: numeric;
get.seqs <- function(barcode.names, barcodes.df, name.i) {
  
  # Prepare vector for sequences
  seqs <- rep(NA, length(barcode.names))
  found.seqs <- 0 # just counter
  
  # Find corresponding sequences
  for (i in 1:nrow(barcodes.df)) {
    if (barcodes.df[[name.i]][i] %in% barcode.names) {
      found.seqs <- found.seqs + 1
      seqs[found.seqs] <- barcodes.df[[name.i+1]][i]
    }
  }

  return(unique(seqs))
}

# == Function check.color.ballance ==
# Function checks color ballance of picked up barcodes.
# :param barcodes.df: dataframe with names of barcodes and their sequences;
# :type barcodes.df: dataframe;
# :param medoids: vector of named of picked up barcodes;
# :type medoids: vector<character>;
# :param double: logical value. TRUE if input file is "double". Else TRUE;
# :type double: logical;
check.color.ballance <- function(barcodes.df, medoids, double) {
  
  # Naive assumption :-)
  ballance <- TRUE
  
  if (! double) {
    # It is all simple in "single" mode:
    
    seqs <- get.seqs(medoids, barcodes.df, 1) # find sequences
    ballance <- is.ballanced(seqs) # validate them for color ballance

  } else {
    
    # Ad hoc trick: name.i is 1 for "single" and 3 for "double" mode.
    # We'll just skip the 2-nd element in this vector in
    #   order not to create new variable and use name.i for indexing.
    get.names.funcs <- c(get.forw.names, NA, get.rev.names)

    for (name.i in c(1, 3)) {
      # Parse names of barcodes
      barcode.names <- as.vector(sapply(X = medoids,
                                      FUN = get.names.funcs[[name.i]]))
      # Fin sequences
      seqs <- get.seqs(barcode.names, barcodes.df, name.i)
      # Validate them for color ballance
      ballance <- ballance && is.ballanced(seqs)
    }
  }
  
  return(ballance)
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
# 'barcodes.df[,1:2]' is passed to 'calc.diss.matr'
#   in order to have same syntax with "double"-formatted work mode.
message("Calculating distances...")
diss.matr <- calc.diss.matr(barcodes.df[,1:2])

# If input file is in "double" format -- calculate dissimilarity matrix
#   for reverse barcodes. Final matrix will be calculated as Euclidean distance
#   between the first and the second matrices.
if (double) {
  reverse.matr <- calc.diss.matr(barcodes.df[,3:4])
  diss.matr <- sqrt(diss.matr ^ 2 + reverse.matr ^ 2) # calculate final matrix
  rm(reverse.matr) # remove!
  
  # Rewrite rownames and colnames for double mode
  rownames(diss.matr) <- paste(barcodes.df$Idx_name_F,
                               barcodes.df$Idx_name_R, sep = ' -- ')
  colnames(diss.matr) <- paste(barcodes.df$Idx_name_F,
                               barcodes.df$Idx_name_R, sep = ' -- ')
}

# Cluster barcodes
message("Clustering...")
pam.clusters <- pam(diss.matr, num.samples, diss = T)

# Get representative barcodes
medoids <- pam.clusters$medoids


# == Color ballance ==
message("Checking color ballance...")
# Check color ballance
ballance <- check.color.ballance(barcodes.df, medoids, double)

# If ballance is corrupted -- inform a user and backup medoids.
if (! ballance) {
  message("Color ballance is corrupted.")
  message("Attepmting to reach color ballance...")
  
  # Backup
  init.medoids <- medoids
  
  # Convert this matrix to data frame -- it will be more convinient
  silh.df <- as.data.frame(pam.clusters$silinfo$widths)
  # Add a column with names
  silh.df$name <- rownames(pam.clusters$silinfo$widths)
}

# while ballance is not reached do:
# 1. Find medoid with the smallest silhouette score. It will be called "old" barcode.
#    It's cluster will have name N.
# 2. Find barcode in the N cluster with the next greatest sulhouette score.
# 3. Replace old barcode with new one.
# 4. Check ballance again.
# end
while (! ballance) {

  # If no more barcodes left -- inform a user, break and use backuped barcodes.
  if (nrow(silh.df) < num.samples) {
    message("Color ballance cannot be reached.")
    medoids <- init.medoids
    break
  }

  # Find medoid with the smallest silhouette score.
  med.silh.scores <- subset(silh.df, silh.df$name %in% medoids)
  min.silh.idx <- which.min(med.silh.scores[,3])
  min.clust.idx <- med.silh.scores$cluster[min.silh.idx]
  
  # Extract N cluster as separate dataframe
  min.cluster <- subset(silh.df, silh.df$cluster == min.clust.idx)
  # Remove old barcode from this dataframe
  min.cluster <- subset(min.cluster, min.cluster$name != med.silh.scores$name[min.silh.idx])
  
  # Get name of old barcode
  old.barcode <- med.silh.scores$name[min.silh.idx]
  # Get name of new barcode
  next.barcode <- subset(min.cluster,
                         min.cluster$sil_width == max(min.cluster$sil_width))$name
  
  # Replace old barcode with new one
  medoids <- replace(medoids, medoids == old.barcode, next.barcode)
  # Remove old barcode from 'silh.df'
  silh.df <- subset(silh.df, silh.df$name != old.barcode)
  
  # Check ballance once again
  ballance <- check.color.ballance(barcodes.df, medoids, double)
}

if (ballance) {
  message("OK")
}

# Print the result
if (! double) {
  cat("\nBest indices:\n")
  for (barcode.name in medoids) {
    cat(sprintf("%s\t%s\n",
                barcode.name,
                barcodes.df$Idx_seq[barcodes.df$Idx_name == barcode.name]))
  }
} else {
  cat("\nBest pairs of indices:\n")
  for (barcode.name in medoids) {
    forw.name <- get.forw.names(barcode.name)
    forw.seq <- unique(barcodes.df$Idx_seq_F[barcodes.df$Idx_name_F == forw.name])
    rev.name <- get.rev.names(barcode.name)
    rev.seq <- unique(barcodes.df$Idx_seq_R[barcodes.df$Idx_name_R == rev.name])
    
    cat(sprintf("%s\t%s\t%s\t%s\n",
                forw.name,forw.seq,
                rev.name, rev.seq))
  }
}
