suppressMessages(library(ShortRead))
suppressMessages(library(stringi))
suppressMessages(library(stringr))
suppressMessages(library(stringdist))
suppressMessages(library(dplyr))
suppressMessages(library(insect))
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))
suppressMessages(library(seqinr))

# Set up command line arguments
option_list = list(make_option(c("-f", "--file"), type = "character",
                               default = NULL, action = "store",
                               help = "FASTQ containing target site reads"),
                   make_option(c("-p", "--pegID"), type = "character",
                               default = NULL, action = "store",
                               help = "pegRNA ID string"),
                   make_option(c("-c","--cellline"), type = "character",
                               default = NULL, action = "store",
                               help = "Sample cell line string"),
                   make_option(c("-d","--day"), type = "integer",
                               default = NULL, action = "store",
                               help = "Day sample was collected [integer]"),
                   make_option(c("-r","--replicate"), type = "integer",
                               default = NULL, action = "store",
                               help = "Replicate number [integer]"),
                   make_option(c("-y","--maxcycle"), type = "integer",
                               default = 80, action = "store",
                               help = "Maximum cycle number where a barcode can end"),
                   make_option(c("-b","--barcode"), type = "character",
                               default = NULL, action = "store",
                               help = "Expected barcode for pegRNA"),
                   make_option(c("-s","--bcstartpos"), type = "integer",
                               default = 54, action = "store",
                               help = "Expected cycle number where barcode begins in read"),
                   make_option(c("-m","--maxdist"), type = "integer",
                               default = 2, action = "store",
                               help = "Max mismatches allowed in barcode before filtering [default = %default]"),
                   make_option(c("-w","--targetregion"), type = "character",
                               default = NULL, action = "store",
                               help = "Expected target region for pegRNA"),
                   make_option(c("-n","--trorientation"), type = "character",
                               default = "forward", action = "store",
                               help = "Whether target site in read is same orientation as in file (forward) or (reverse)"),
                   make_option(c("-t","--trstartpos"), type = "integer",
                               default = 7, action = "store",
                               help = "Expected cycle number where target region begins in the read"),
                   make_option(c("-q","--minfreq"), type = "double",
                               default = 0.001, action = "store",
                               help = "Minimum frequncy for outcome filter [default = %default]"),
                   make_option(c("-a","--altminreads"), type = "integer",
                               default = 10, action = "store",
                               help = "Alternative outcome filter for low coverage pegs [default = %default]"),
                   make_option(c("-l","--basequal"), type = "integer",
                               default = 20, action = "store",
                               help = "Minimum average base quality to not be corrected [default = %default]"),
                   make_option(c("-o","--output"), type = "character",
                               default = ".", action = "store",
                               help = "Directory path for output files [default = current directory]"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Pull info from arguments
fq.path <- opt$file #path containing reads to be analyzed
pegID <- opt$pegID #pegRNA ID string
cl <- opt$cellline #cell line
day <- paste0("Day ", opt$day) #collection day
repl <- paste0("R", opt$replicate) #replicate
max.pos <- opt$maxcycle #latest cycle in a read where a barcode is allowed to begin
exp.bc <- opt$barcode #expected barcode sequence
bc.start <- opt$bcstartpos #expected cycle in read where barcode begins
max.dist <- opt$maxdist #max distance allowed in barcode
exp.wtr <- opt$targetregion #expected wide target region
wtr.start <- opt$trstartpos #expected cycle in read where target region begins
orient <- opt$trorientation #orientation of target region in read
min.freq <- opt$minfreq #minimum observed frequency to keep an outcome, %/100
min.reads <- opt$altminreads #alternative outcome filter for low coverage (<10000x) pegs
min.bq <- opt$basequal #minimum base quality to not correct to expected
out.path <- opt$output #directory for storing outputs

# Check all required components provided
if(is.null(opt$file) | is.null(opt$pegID) | is.null(opt$cellline) |
   is.null(opt$day) | is.null(opt$replicate) | is.null(opt$barcode) |
   is.null(opt$targetregion)){
  print_help(opt_parser)
  stop("Please provide required arguments.", call. = FALSE)
  quit(save = "no")
}

# Extract intended edit from pegID
int.edit <- gsub("^[^_]*_[^_]*_","", pegID)
int.edit <- gsub("_.*","", int.edit)
des.pos <- as.numeric(gsub("[ACGT]to[ACGT]","", int.edit))
des.edit <- gsub("^[0-9+]","", int.edit)
des.edit <- gsub("to",">", des.edit)
# Extract gene name from pegID
gene.name <- gsub("_.*", "", pegID)

print(paste0("Analyzing file ", fq.path))
print(paste0("Expected barcode for peg ", pegID, " is ", exp.bc))
print(paste0("Expected target region for peg ", pegID, " is ", exp.wtr))
print(paste0("Output files are in directory ", out.path))

##################################
### STEP 1 - Barcode Filtering ###
##################################

# Updated code from package author to support indels with vmatchPattern
# Taken from: https://support.bioconductor.org/p/58350/
vmatchPattern2 <- function(pattern, subject,
                           max.mismatch=0, min.mismatch=0,
                           with.indels=FALSE, fixed=TRUE,
                           algorithm="auto")
{
  if (!is(subject, "XStringSet"))
    subject <- Biostrings:::XStringSet(NULL, subject)
  algo <- Biostrings:::normargAlgorithm(algorithm)
  if (Biostrings:::isCharacterAlgo(algo))
    stop("'subject' must be a single (non-empty) string ",
         "for this algorithm")
  pattern <- Biostrings:::normargPattern(pattern, subject)
  max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
  min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch,
                                                  max.mismatch)
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  algo <- Biostrings:::selectAlgo(algo, pattern,
                                  max.mismatch, min.mismatch,
                                  with.indels, fixed)
  C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
                  max.mismatch, min.mismatch,
                  with.indels, fixed, algo,
                  "MATCHES_AS_RANGES",
                  PACKAGE="Biostrings")
  unlisted_ans <- IRanges(start=unlist(C_ans[[1L]],
                                       use.names=FALSE),
                          width=unlist(C_ans[[2L]],
                                       use.names=FALSE))
  relist(unlisted_ans, C_ans[[1L]])
}

print("Importing reads...")

# Import reads
reads <- readFastq(fq.path)

print(paste0(length(reads), " reads imported. Performing barcode filtering..."))

# Fuzzy match to expacted barcode
bc.matches <- as.data.frame(as(vmatchPattern2(rc(exp.bc), sread(reads), max.mismatch = floor(nchar(exp.bc)/2)),
                            "CompressedIRangesList"))
# Extract barcode match closest to expected position (bases 54 - 70)
# DOUBLECHECK: CAN END BE end[which.min(abs(start - bc.start))]??
bcs.sum <- bc.matches %>% 
                group_by(group) %>% 
                summarise(start = start[which.min(abs(start - bc.start))], end = start + nchar(exp.bc) - 1)

# Identify reads with no barcode match
no.bc <- setdiff(1:length(reads), bcs.sum$group)
# Identify reads with barcode match outside bounds of read
too.deep <- setdiff(1:length(reads), bcs.sum$group[bcs.sum$end <= max.pos])
too.early <- setdiff(1:length(reads), bcs.sum$group[bcs.sum$start >= 1])
# Remove reads with no barcode match
reads <- reads[setdiff(1:length(reads), c(no.bc, too.deep, too.early))]
bcs.sum <- bcs.sum[bcs.sum$start >= 1 & bcs.sum$end <= max.pos,]

print(paste0("Removed ", length(no.bc), " reads with no barcode match"))
print(paste0("Removed ", length(too.deep), " reads with barcode match too late in read (end position > max.pos)"))
print(paste0("Removed ", length(too.early), " reads with barcode match too early in read (start position < 1)"))


# Get barcode
bcs.sum$bc <- rc(str_sub(sread(reads), bcs.sum$start,bcs.sum$end))
# Get distance to expected barcode, use optimal string alignment to allow for indels
bcs.sum$dist <- stringdist(bcs.sum$bc, exp.bc)

print("Summarizing barcode quality scores...")

bcs.sum$min <- NA
bcs.sum$quant1 <- NA
bcs.sum$med <- NA
bcs.sum$mean<- NA
bcs.sum$quant3 <- NA
bcs.sum$max <- NA

# Summarize quality scores in identified barcode
for(i in 1:nrow(bcs.sum)){
  if(is.na(bcs.sum$start[i])){
    next;
  }
  
  tmp <- summary(as(quality(reads)[i],"matrix")[,bcs.sum$start[i]:bcs.sum$end[i]])
  bcs.sum$min[i] <- tmp[1]
  bcs.sum$quant1[i] <- tmp[2]
  bcs.sum$med[i] <- tmp[3]
  bcs.sum$mean[i] <- tmp[4]
  bcs.sum$quant3[i] <- tmp[5]
  bcs.sum$max[i] <- tmp[6]
}

# Used for plotting - finds longest set of distances with max qual <= 33
dist.maxes <- bcs.sum %>% group_by(dist) %>% summarise(M = max(mean, na.rm = TRUE))
dist.ndx <- dist.maxes$dist[dist.maxes$M <= 33]
dist.seqs <- split(dist.ndx, cumsum(c(1, diff(dist.ndx) != 1)))
longest.subseq <- which(unlist(lapply(dist.seqs, length)) == max(unlist(lapply(dist.seqs, length))))
if(length(longest.subseq) > 1){
  longest.subseq <- max(longest.subseq)
}
bc.max.dist <- max(dist.seqs[[longest.subseq]])

# Filter reads with edit distance > max.dist
reads.filt <- reads[bcs.sum$dist <= max.dist]

# Plot and save edit distance vs mean quality
mean.dist.cor <- round(cor(bcs.sum$dist, bcs.sum$mean), 3)
bcs.sum$dist <- as.factor(bcs.sum$dist)

# Helper function for adding read count annotations
dist_length <- function(x){
  return(data.frame(y=min(bcs.sum$mean) - 0.5, label = paste0("n = ", length(x))))
}

bc.qual.plt <- ggplot(bcs.sum, aes(dist, mean)) + geom_boxplot() +
                  theme_bw() + xlab("Number of mismatches in barcode") +
                  ylab("Average base sequence quality in barcode") +
                  ggtitle(paste0("Barcode Quality for ", gene.name, 
                                 " - Edit +",des.pos,":",des.edit, " - ", cl,
                                 " - ", day, " - ", repl)) +
                  stat_summary(aes(x = dist),
                               fun.data = dist_length, 
                               geom = "text",
                               size = 3) + 
                  annotate(geom = "text", x = as.character(bc.max.dist), y = 35, 
                           label = paste0("PCC = ", mean.dist.cor),
                           size = 6, hjust = 1) +
                  theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())

ggsave(paste0(out.path, "/barcode_quality_",pegID,".pdf"),
       bc.qual.plt,
       device = "pdf",
       height = 6,
       width = 8)

print("Step 1 - Barcode filtering - complete...")

###############################################
### STEP 2 -  Correct and Collapse Outcomes ###
###############################################

# Fuzzy match to expected wide target region (WTR)
if(orient == "forward"){
  wtr.matches <- as.data.frame(as(vmatchPattern2(rc(exp.wtr), 
                                                 sread(reads.filt),
                                                 with.indels = TRUE,
                                                 max.mismatch = floor(.4*nchar(exp.wtr)), #at least 60% seq identity 
                                                 "CompressedIRangesList")))  
} else {
  wtr.matches <- as.data.frame(as(vmatchPattern2(exp.wtr, 
                                                 sread(reads.filt),
                                                 with.indels = TRUE,
                                                 max.mismatch = floor(.4*nchar(exp.wtr)), #at least 60% seq identity 
                                                 "CompressedIRangesList")))
}


# Extract WTR
wtr.sum <- wtr.matches %>% 
  group_by(group) %>% 
  summarise(start = start[which.min(abs(start - wtr.start))], end = end[which.min(abs(start - wtr.start))])

# Identify reads with no target region match
no.wtr <- setdiff(1:length(reads.filt), wtr.sum$group)
# Remove reads with no target region match
reads.filt <- reads[setdiff(1:length(reads.filt), no.wtr)]

print(paste0("Removed ", length(no.wtr), " reads with no match to target region..."))

# Get WTR
if(orient == "forward"){
  wtr.sum$wtr <- rc(str_sub(sread(reads.filt), wtr.sum$start,wtr.sum$end))  
} else {
  wtr.sum$wtr <- str_sub(sread(reads.filt), wtr.sum$start, wtr.sum$end)
}


# Get frequencies of unique outcomes
wtr.unique <- wtr.sum %>% group_by(wtr) %>% summarise(N = n())
wtr.unique$freq <- wtr.unique$N/sum(wtr.unique$N)

# Set min freq to either min.freq (>10000x coverage) or 10 reads (<10000x)
min.freq <- max(min.freq, min.reads/sum(wtr.unique$N))

# Filter low frequency outcomes
wtr.filt <- wtr.unique[wtr.unique$freq >= min.freq,]

if(nrow(wtr.filt) == 0){
  print("ERROR: NO OUTCOMES SURVIVE MIN FREQUENCY FILTER. LOW COVERAGE PEG.")
  print("SAVING OUTCOMES AND FILTERS...")
  # Save detected outcomes
  wtr.unique$pegID <- pegID
  wtr.unique$CellLine <- cl
  wtr.unique$Day <- day
  wtr.unique$Replicate <- repl
  write.table(wtr.unique[order(wtr.unique$freq, decreasing = TRUE),
                         c("pegID","CellLine","Day","Replicate","wtr","N","freq")],
              paste0(out.path, "/all_outcomes_", pegID, ".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
  # Save filtering statistics
  write.table(data.frame(pegID = pegID,
                         CellLine = cl,
                         Day = day,
                         Replicate = repl,
                         Qual.Dist.Cor = mean.dist.cor,
                         Coverage = length(reads) + length(no.bc),
                         No.BC.Detected = length(no.bc),
                         N.BC.Filter = length(reads) - length(reads.filt),
                         N.BC.Survive = length(reads.filt),
                         N.Outcomes = nrow(wtr.unique),
                         N.LowFreqOutcomes = nrow(wtr.unique) - nrow(wtr.filt),
                         N.Outcomes.Survive = nrow(wtr.filt),
                         N.CorrectedBases = 0,
                         N.CorrectedReads = 0,
                         N.CorrectedOutcomes = 0),
              paste0(out.path, "/filter_results_",pegID,".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
  # Save error to flag file
  writeLines(c(paste(pegID, cl, day, repl, "LOW COVERAGE - NO OUTCOMES PASS MIN FREQ", sep = "\t")),
             paste0(out.path,"/flag_",pegID,".txt"))
  
  print(paste0(pegID, " analyzed. Script stopped early due to no outcomes surviving minimum frequency filter (low coverage). Script complete."))
  quit(save = "no")
}

print("Correcting low quality bases in target region...")

# Identify and correct low quality bases
# Tracking variables for reporting
total.corr.lq <- 0
total.corr.reads <- 0
for(i in 1:nrow(wtr.filt)){
  # Identify reads mapping to each outcome
  read.matches <- which(wtr.sum$wtr == wtr.filt$wtr[i])
  # Extract quality scores for full read
  tmp.qual <- as(quality(reads.filt)[read.matches], "matrix")
  # Extract WTR position in read
  tmp.wtr.pos <- wtr.sum[read.matches,]
  # If positions are all the same, save work
  if(nrow(unique(tmp.wtr.pos[,c("start", "end")])) == 1){
    new.qual <- tmp.qual[,seq(tmp.wtr.pos$end[1], tmp.wtr.pos$start[1], -1)]
  } else {
    # Otherwise do work
    for(j in 1:nrow(tmp.qual)){
      # Going row-by-row, extract quality scores by position in each read
      # Note sequence goes backwards because the WTR has been reverse
      # complemented, but the quality score orientation has not
      if(j == 1){
        new.qual <- tmp.qual[j,seq(tmp.wtr.pos$end[j], tmp.wtr.pos$start[j], -1)]
      } else {
        new.qual <- rbind(new.qual,
                          tmp.qual[j,seq(tmp.wtr.pos$end[j], tmp.wtr.pos$start[j], -1)])
      }
    }
  }
  # Calculate mean quality at each position
  pos.means <- colMeans(new.qual)
  # Identify low quality bases
  low.qual.pos <- which(pos.means <= min.bq)
  # If low quality bases exist, check to correct
  if(length(low.qual.pos) != 0){
    # Check if low quality bases match expected bases, correct to expected base if not
    if(sum(str_split(wtr.filt$wtr[i], "")[[1]][low.qual.pos] != str_split(exp.wtr, "")[[1]][low.qual.pos]) > 0){
      # Split outcome string into bases
      tmp.corr.outcome <- str_split(wtr.filt$wtr[i], "")[[1]]
      # Add to statistics - must check before correcting
      total.corr.lq <- total.corr.lq + sum(tmp.corr.outcome[low.qual.pos] != str_split(exp.wtr, "")[[1]][low.qual.pos])
      # Correct bases
      tmp.corr.outcome[low.qual.pos] <- str_split(exp.wtr, "")[[1]][low.qual.pos]
      wtr.filt$wtr[i] <- paste(tmp.corr.outcome, collapse = "")
      # Document bases corrected
      cat(paste0("Corrected ", length(low.qual.pos), 
                 " low quality base(s) with mean base sequence quality scores:\n",
                 paste(round(pos.means[low.qual.pos], 2), collapse = " "),"\n",
                 "for ", wtr.filt$N[i], " reads...\n"))
      # Add to statistics
      total.corr.reads <- total.corr.reads + wtr.filt$N[i]
    }
  }
  
}
print(paste0("Total corrected low quality bases: ", total.corr.lq))
print(paste0("Total corrected reads: ", total.corr.reads))

# Re-collapse outcomes and re-calculate frequencies
corrected.wtr <- wtr.filt %>% group_by(wtr) %>% summarise(N = sum(N))
corrected.wtr$freq <- corrected.wtr$N/sum(corrected.wtr$N)

# Order for plotting
changed.wtr <- corrected.wtr[corrected.wtr$wtr != exp.wtr,]
changed.wtr <- changed.wtr[order(changed.wtr$freq, decreasing = TRUE),]
corrected.wtr <- rbind(corrected.wtr[corrected.wtr$wtr == exp.wtr,],
                       changed.wtr)

# Add meta info before saving
wtr.unique$pegID <- pegID
wtr.unique$CellLine <- cl
wtr.unique$Day <- day
wtr.unique$Replicate <- repl
corrected.wtr$pegID <- pegID
corrected.wtr$CellLine <- cl
corrected.wtr$Day <- day
corrected.wtr$Replicate <- repl
corrected.wtr$id <- paste0("Outcome_",1:nrow(corrected.wtr))

# Save outcomes
write.table(wtr.unique[order(wtr.unique$freq, decreasing = TRUE),
                       c("pegID","CellLine","Day","Replicate","wtr","N","freq")],
            paste0(out.path, "/all_outcomes_", pegID, ".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(corrected.wtr[,c("pegID","CellLine","Day","Replicate","id","wtr","N","freq")], 
            paste0(out.path, "/corrected_outcome_frequencies_",pegID,".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.fasta(as.list(corrected.wtr$wtr), 
            paste(corrected.wtr$pegID,
                  corrected.wtr$CellLine,
                  corrected.wtr$Day,
                  corrected.wtr$Replicate,
                  corrected.wtr$id, sep = ":"),
            paste0(out.path,"/corrected_outcomes_", pegID, ".fasta"))

# Save filtering statistics
write.table(data.frame(pegID = pegID,
                       CellLine = cl,
                       Day = day,
                       Replicate = repl,
                       Qual.Dist.Cor = mean.dist.cor,
                       Coverage = length(reads) + length(no.bc),
                       No.BC.Detected = length(no.bc),
                       N.BC.Filter = length(reads) - length(reads.filt),
                       N.BC.Survive = length(reads.filt),
                       N.Outcomes = nrow(wtr.unique),
                       N.LowFreqOutcomes = nrow(wtr.unique) - nrow(wtr.filt),
                       N.Outcomes.Survive = nrow(wtr.filt),
                       N.CorrectedBases = total.corr.lq,
                       N.CorrectedReads = total.corr.reads,
                       N.CorrectedOutcomes = nrow(corrected.wtr)),
            paste0(out.path, "/filter_results_",pegID,".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

print("Step 2 - Correcting and collapsing outcomes - is complete...")

#######################################
### STEP 3 -  Characterize Outcomes ###
#######################################

print("Characterizing outcomes...")

for(i in 1:nrow(corrected.wtr)){
  # Align outcome to expected wtr
  tmp.align <- pairwiseAlignment(corrected.wtr$wtr[i], exp.wtr)
  # Identify mismatches
  mismatches <- mismatchSummary(tmp.align)$subject[,1:3]
  # Correct position to be relative to the cut site
  mismatches$SubjectPosition[mismatches$SubjectPosition < 22] <- mismatches$SubjectPosition[mismatches$SubjectPosition < 22] - 22
  mismatches$SubjectPosition[mismatches$SubjectPosition >= 22] <- mismatches$SubjectPosition[mismatches$SubjectPosition >= 22] - 21
  
  # Add meta info
  if(nrow(mismatches) != 0){
    mismatches$id <- corrected.wtr$id[i]
    mismatches$freq <- corrected.wtr$freq[i]  
  }
  
  if(!exists("mismatch.res")){
    mismatch.res <- mismatches
  } else {
    mismatch.res <- rbind(mismatch.res, mismatches)
  }
  
  # Identify insertions
  insertions <- as.data.frame(unlist(insertion(tmp.align)))
  
  
  if(nrow(insertions) != 0){
    for(j in 1:nrow(insertions)){
      insertions$edit[j] <- paste("+",str_sub(pattern(tmp.align), insertions$start[j], insertions$end[j]), 
                                  collapse = "")
    }
    # Correct position to be relative to the cut site
    insertions$start[insertions$start < 22] <- insertions$start[insertions$start < 22] - 22
    insertions$start[insertions$start >= 22] <- insertions$start[insertions$start >= 22] - 21
    insertions$end[insertions$end < 22] <- insertions$end[insertions$end < 22] - 22
    insertions$end[insertions$end >= 22] <- insertions$end[insertions$end >= 22] - 21
    # Add meta info
    insertions$id <- corrected.wtr$id[i]
    insertions$freq <- corrected.wtr$freq[i]  
  }
  
  if(!exists("insertion.res")){
    insertion.res <- insertions
  } else {
    insertion.res <- rbind(insertion.res, insertions)
  }
  
  # Identify deletions
  deletions <- as.data.frame(unlist(deletion(tmp.align)))
  
  
  if(nrow(deletions) != 0){
    for(j in 1:nrow(deletions)){
      deletions$edit[j] <- paste("-",str_sub(subject(tmp.align), deletions$start[j], deletions$end[j]), 
                                 collapse = "")
    }
    # Correct position to be relative to the cut site
    deletions$start[deletions$start < 22] <- deletions$start[deletions$start < 22] - 22
    deletions$start[deletions$start >= 22] <- deletions$start[deletions$start >= 22] - 21
    deletions$end[deletions$end < 22] <- deletions$end[deletions$end < 22] - 22
    deletions$end[deletions$end >= 22] <- deletions$end[deletions$end >= 22] - 21
    # Add meta info
    deletions$id <- corrected.wtr$id[i]
    deletions$freq <- corrected.wtr$freq[i]  
  }
  
  if(!exists("deletion.res")){
    deletion.res <- deletions
  } else {
    deletion.res <- rbind(deletion.res, deletions)
  }
}

print("Characterization complete. Combining edits...")

if(nrow(mismatch.res) > 0){
  # Make consistent edit name
  mismatch.res$edit <- paste0(mismatch.res$Subject, ">", mismatch.res$Pattern)
  # Prep characterizations to combine
  colnames(mismatch.res)[1] <- "position"
  mismatch.res$type <- "Mismatch"
}

if(nrow(insertion.res) > 0){
  colnames(insertion.res)[1] <- "position"
  insertion.res$type <- "Insertion"
}

if(nrow(deletion.res) > 0){
  colnames(deletion.res)[1] <- "position"
  deletion.res$type <- "Deletion"
}

# Combine all characterizations
if(nrow(mismatch.res) > 0 & nrow(insertion.res) > 0 & nrow(deletion.res) > 0){
  all.edits <- rbind(mismatch.res[,c("id","position", "edit", "freq","type")],
                     insertion.res[,c("id","position", "edit", "freq", "type")],
                     deletion.res[,c("id","position", "edit", "freq", "type")])
} else if(nrow(mismatch.res) > 0 & nrow(insertion.res) > 0 & nrow(deletion.res) == 0){
  all.edits <- rbind(mismatch.res[,c("id","position", "edit", "freq","type")],
                     insertion.res[,c("id","position", "edit", "freq", "type")])
} else if(nrow(mismatch.res) > 0 & nrow(insertion.res) == 0 & nrow(deletion.res) > 0){
  all.edits <- rbind(mismatch.res[,c("id","position", "edit", "freq","type")],
                     deletion.res[,c("id","position", "edit", "freq", "type")])
} else if(nrow(mismatch.res) == 0 & nrow(insertion.res) > 0 & nrow(deletion.res) > 0){
  all.edits <- rbind(insertion.res[,c("id","position", "edit", "freq","type")],
                     deletion.res[,c("id","position", "edit", "freq", "type")])
} else if(nrow(mismatch.res) > 0 & nrow(insertion.res) == 0 & nrow(deletion.res) == 0){
  all.edits <- mismatch.res[,c("id","position", "edit", "freq","type")]
} else if(nrow(mismatch.res) == 0 & nrow(insertion.res) > 0 & nrow(deletion.res) == 0){
  all.edits <- insertion.res[,c("id","position", "edit", "freq","type")]
} else if(nrow(mismatch.res) == 0 & nrow(insertion.res) == 0 & nrow(deletion.res) > 0){
  all.edits <- deletion.res[,c("id","position", "edit", "freq","type")]
} else if(nrow(mismatch.res) == 0 & nrow(insertion.res) == 0 & nrow(deletion.res) == 0){
  print("ERROR: No editing detected...all reads are wildtype.")
  print("Saving statistics and flag to file...")
  
  write.table(data.frame(pegID = pegID,
                         CellLine = cl,
                         Day = day,
                         Replicate = repl,
                         TotalReads = sum(corrected.wtr$N),
                         Freq.WildType = corrected.wtr$freq[1],
                         Freq.TotalEdited = 0,
                         Freq.TotalOffTarget = 0,
                         Freq.EditOnly = 0,
                         Freq.EditPlusOtherMM = 0,
                         Freq.OtherMMOnly = 0,
                         Freq.EditPlusInsertion = 0,
                         Freq.InsertionOnly = 0,
                         Freq.EditPlusDeletion = 0,
                         Freq.DeletionOnly = 0),
              paste0(out.path,"/summarized_outcomes_",pegID,".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
  # Save error to flag file
  writeLines(c(paste(pegID, cl, day, repl, "NO EDITING DETECTED - ALL READS WILDTYPE", sep = "\t")),
             paste0(out.path,"/flag_",pegID,".txt"))
  
  print(paste0(pegID, " analyzed. Script stopped early due to no edits passing filters. Script complete."))
  
  quit(save = "no")
}

all.edits <- all.edits[order(as.numeric(gsub("Outcome_","",all.edits$id)), all.edits$position),]

# Identify editing events which match desired edit
all.edits$IsEdit <- (all.edits$position == des.pos & all.edits$edit == des.edit)
# Set column indicating whether edit co-occurs with desired edit
all.edits$cooccur <- FALSE
for(i in 1:nrow(all.edits)){
  # Get all edits for an outcome
  tmp <- all.edits[all.edits$id == all.edits$id[i],]
  # Go next if only one edit
  if(nrow(tmp) == 1){
    next
  }
  # Go next if no edit is intended edit
  if(sum(tmp$IsEdit) == 0){
    next
  }
  # Remove intended edit
  tmp <- tmp[!tmp$IsEdit,]
  # Set cooccur to TRUE for remaining edits
  for(j in 1:nrow(tmp)){
    all.edits$cooccur[all.edits$position == tmp$position[j] &
                        all.edits$edit == tmp$edit[j] &
                        all.edits$id == tmp$id[j]] <- TRUE
  }
}

print("Summarizing editing outcomes...")

# Summarize across outcomes
edits.sum <- all.edits %>% 
              group_by(position, edit, type) %>% 
              summarise(freq = sum(freq), 
                        N = n(), 
                        cooccur = any(cooccur),
                        only.co = sum(cooccur) == n())
edits.sum$type <- factor(edits.sum$type, 
                         levels = unique(edits.sum$type)[order(unique(edits.sum$type), decreasing = TRUE)]) 

# Summarize types of edits for each outcome
outcome.char <- all.edits %>% 
                  group_by(id) %>% 
                  summarise(HasEdit = sum(IsEdit), 
                            HasOtherMM = sum(!IsEdit & type == "Mismatch"), 
                            HasIns = sum(type == "Insertion") > 0, 
                            HasDel = sum(type == "Deletion") > 0, 
                            freq = freq)
outcome.char <- unique(outcome.char)
outcome.char$EditOnly <- outcome.char$HasEdit > 0 & 
                          outcome.char$HasOtherMM == 0 & 
                          !outcome.char$HasIns & !outcome.char$HasDel
outcome.char$EditOtherMM <- outcome.char$HasEdit > 0 & 
                              outcome.char$HasOtherMM > 0 &
                              !outcome.char$HasIns & !outcome.char$HasDel
outcome.char$OtherOnly <- outcome.char$HasEdit == 0 &
                            !outcome.char$HasIns & !outcome.char$HasDel

outcome.char <- outcome.char[order(as.numeric(gsub("Outcome_","", outcome.char$id))),]
outcome.char <- rbind(data.frame(id = "Outcome_1", HasEdit = 0, HasOtherMM = 0,
                                 HasIns = FALSE, HasDel = FALSE, freq = corrected.wtr$freq[1],
                                 EditOnly = FALSE, EditOtherMM = FALSE, OtherOnly = FALSE),
                      outcome.char)

# Calculate frequencies of different editing event combinations
### Edit only
eonly <- sum(outcome.char$freq[outcome.char$EditOnly])
### Edit + other mm
eother <- sum(outcome.char$freq[outcome.char$EditOtherMM])
### Other mm only
oonly <- sum(outcome.char$freq[outcome.char$OtherOnly])
### Edit + ins
insedit <- sum(outcome.char$freq[outcome.char$HasIns & outcome.char$HasEdit > 0])
### Ins only
insonly <- sum(outcome.char$freq[outcome.char$HasIns & outcome.char$HasEdit == 0])
### Edit + del
deledit <- sum(outcome.char$freq[outcome.char$HasDel & outcome.char$HasEdit > 0])
### Del only
delonly <- sum(outcome.char$freq[outcome.char$HasDel & outcome.char$HasEdit == 0])
### Frequency of intended edit
etotal <- sum(outcome.char$freq[outcome.char$HasEdit > 0])
### Frequency of unintended edits alone
uetotal <- sum(outcome.char$freq[outcome.char$HasEdit == 0 & outcome.char$id != "Outcome_1"])
### Frequency of WT (no edit)
wttotal <- corrected.wtr$freq[1]

# Prepare data for plotting
# Set column with base after editing
edits.sum$alt.base <- NA
edits.sum$alt.base[edits.sum$type == "Mismatch"] <- gsub(".>","",edits.sum$edit[edits.sum$type == "Mismatch"])
edits.sum$alt.base[is.na(edits.sum$alt.base)] <- edits.sum$edit[is.na(edits.sum$alt.base)]

# Add meta info before saving
all.edits$pegID <- pegID
all.edits$CellLine <- cl
all.edits$Day <- day
all.edits$Replicate <- repl
edits.sum$pegID <- pegID
edits.sum$CellLine <- cl
edits.sum$Day <- day
edits.sum$Replicate <- repl
outcome.char$pegID <- pegID
outcome.char$CellLine <- cl
outcome.char$Day <- day
outcome.char$Replicate <- repl

print("Saving outcomes...")

# Save outputs
write.table(all.edits[,c("pegID","CellLine","Day","Replicate","id", "position", "edit",
                         "type","freq", "IsEdit","cooccur")], 
            paste0(out.path, "/characterized_outcomes_",pegID,".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(edits.sum[,c("pegID","CellLine","Day","Replicate", "position", "edit",
                         "type", "freq", "N", "cooccur", "only.co", "alt.base")], 
            paste0(out.path, "/characterized_outcomes_",pegID,"_summary.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(outcome.char[,c("pegID","CellLine","Day","Replicate","id","HasEdit","HasOtherMM",
                            "HasIns","HasDel","freq","EditOnly","EditOtherMM","OtherOnly")],
            paste0(out.path, "/outcome_characteristics_",pegID,"_summary.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(data.frame(pegID = pegID,
                       CellLine = cl,
                       Day = day,
                       Replicate = repl,
                       TotalReads = sum(corrected.wtr$N),
                       Freq.TotalWildType = wttotal,
                       Freq.TotalEdited = etotal,
                       Freq.TotalOffTarget = uetotal,
                       Freq.EditOnly = eonly,
                       Freq.EditPlusOtherMM = eother,
                       Freq.OtherMMOnly = oonly,
                       Freq.EditPlusInsertion = insedit,
                       Freq.InsertionOnly = insonly,
                       Freq.EditPlusDeletion = deledit,
                       Freq.DeletionOnly = delonly),
            paste0(out.path,"/summarized_outcomes_",pegID,".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

if(sum(edits.sum$type == "Deletion") > 10 | sum(edits.sum$type == "Insertion") > 10){
  # Save error to flag file
  writeLines(c(paste(pegID, cl, day, repl, "HIGH INDELS ( > 10 ) - EDITS NOT PLOTTED", sep = "\t")),
             paste0(out.path,"/flag_",pegID,".txt"))
  
  print(paste0(pegID, " analyzed. Script stopped early due more than 10 insertions or 10 deletions present; edits not plotted. Script complete."))
  
  quit(save = "no")
}

print("Plotting outcomes...")

# Colors for x-axis
# Blue for protospacer, red for target base, green for RT template
# Note that the RT template starts at position +1, overlapping the protospacer
axis.colors <- c(rep("black", 5), 
                 rep("dodgerblue4", 19), 
                 "forestgreen", 
                 "firebrick", 
                 rep("forestgreen",as.numeric(gsub(".*_","",pegID))-5))
axis.colors <- c(axis.colors, rep("black",47-length(axis.colors)))

# Colors for edit types
# Sourced from wesanderson palettes r package
mismatch.colors <- c(A = "#00A08A", C = "#046C9A", G = "#F98400", T = "#B40F20")
insertion.colors <- c("#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000",
                      "#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
deletion.colors <- c("#E1BD6D", "#EABE94", "#0B775E", "#35274A" ,"#F2300F",
                     "#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")

# Co-occurence dfs
# These are needed to place the red/gray stars underneath outcomes which
# indicate co-occurence of an edit with the intended edit, little messy
if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Mismatch") > 0 &
   sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Mismatch") > 0){
  cooccur.mm <- rbind(data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Mismatch","position"],
                                 y = -0.02,
                                 label = '*'),
                      data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Mismatch","position"],
                                 y = -0.02,
                                 label = '+'))
} else if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Mismatch") > 0){
  cooccur.mm <- data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Mismatch","position"],
                           y = -0.02,
                           label = '*')
} else if(sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Mismatch") > 0){
  cooccur.mm <- data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Mismatch","position"],
                           y = -0.02,
                           label = '+')
} else {
  cooccur.mm <- NA
}
if(!is.null(nrow(cooccur.mm))){
  cooccur.mm$color <- ifelse(cooccur.mm$label == "*", "red", "#999999")  
}


if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Insertion") > 0 &
   sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Insertion") > 0){
  cooccur.ins <- rbind(data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Insertion","position"],
                                  y = -0.02,
                                  label = '*'),
                       data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Insertion","position"],
                                  y = -0.02,
                                  label = '+'))
} else if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Insertion") > 0){
  cooccur.ins <- data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Insertion","position"],
                            y = -0.02,
                            label = '*')
} else if(sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Insertion") > 0){
  cooccur.ins <- data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Insertion","position"],
                            y = -0.02,
                            label = '+')
} else {
  cooccur.ins <- NA
}
if(!is.null(nrow(cooccur.ins))){
  cooccur.ins$color <- ifelse(cooccur.ins$label == "*", "red", "#999999")
}

if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Deletion") > 0 &
   sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Deletion") > 0){
  cooccur.del <- rbind(data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Deletion","position"],
                                  y = -0.02,
                                  label = '*'),
                       data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Deletion","position"],
                                  y = -0.02,
                                  label = '+'))
} else if(sum(edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Deletion") > 0){
  cooccur.del <- data.frame(x = edits.sum[edits.sum$cooccur & edits.sum$only.co & edits.sum$type == "Deletion","position"],
                            y = -0.02,
                            label = '*')
} else if(sum(edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Deletion") > 0){
  cooccur.del <- data.frame(x = edits.sum[edits.sum$cooccur & !edits.sum$only.co & edits.sum$type == "Deletion","position"],
                            y = -0.02,
                            label = '+')
} else {
  cooccur.del <- NA
}
if(!is.null(nrow(cooccur.del))){
  cooccur.del$color <- ifelse(cooccur.del$label == "*", "red", "#999999")
}


# Creates a list of plots, one for each edit type found
ggList <- lapply(split(edits.sum, edits.sum$type), function(i){
            ggplot(i, aes(position, freq, group = 1, fill = alt.base)) +
              geom_bar(stat = "identity", width = 1) +
              theme_bw() +
              facet_grid(type ~ .) +
              geom_hline(yintercept = .01, color = "black", linetype = "dotted") +
              geom_text(aes(label = N), size = 3, vjust = -.4) +
              expand_limits(y = max(edits.sum$freq) + 0.05) +
              scale_x_continuous(breaks = c(-21:-1,1:26),
                                 labels = paste(c(-21:-1,1:26),str_split(exp.wtr, "")[[1]], sep = "\n"),
                                 limits = c(-21.5,26.5)) +
              theme(axis.text.x = element_text(size = 6, color = axis.colors),
                    panel.grid = element_blank(),
                    panel.background = element_blank()) +
              xlab("") + ylab("") 
})

# Modify plots generated by lapply depending on which ones are present

if(!is.null(ggList$Mismatch)){
  # Mismatch plot gets plot title
  ggList$Mismatch <- ggList$Mismatch + 
    annotate(geom = "text", x = 25, y = .55, 
             label = paste0("Desired edit (+", 
                            des.pos, ":", des.edit, 
                            ") only: ", round(eonly* 100,2), "%"), 
             size = 3, hjust = 1) + 
    annotate(geom = "text", x = 25, y = .48, 
             label = paste0("Desired edit + other mismatch: ", round(eother*100,2), "%"), 
             size = 3, hjust = 1) + 
    annotate(geom = "text", x = 25, y = .41, 
             label = paste0("Other mismatch only: ", round(oonly*100,2), "%"), 
             size = 3, hjust = 1) + 
    annotate(geom = "text", x = -20, y = .55,
             label = paste0("Total with desired edit: ", round(etotal*100,2), "%"),
             size = 3, hjust = 0) +
    annotate(geom = "text", x = -20, y = .48,
             label = paste0("Total off-target edit: ", round(uetotal*100,2), "%"),
             size = 3, hjust = 0) +
    annotate(geom = "text", x = -20, y = .41,
             label = paste0("Total wildtype: ", round(wttotal*100,2), "%"),
             size = 3, hjust = 0) +
    annotate(geom = "text", x = -20, y = .34,
             label = paste0("Total reads: ", sum(corrected.wtr$N)),
             size = 3, hjust = 0) +
    labs(fill = "Base after\nEditing") +
    scale_fill_manual(values = mismatch.colors, limits = force) +
    ggtitle(paste0("Editing Outcomes for ", gene.name, 
                   " - Edit +",des.pos,":",des.edit, " - ", cl,
                   " - ", day, " - ", repl)) 
  
  if(!is.null(nrow(cooccur.mm))){
    ggList$Mismatch <- ggList$Mismatch + annotate(geom = "text", x = cooccur.mm$position, y = -.04, 
                                                  color = cooccur.mm$color, label = "*") 
  }
}

if(!is.null(ggList$Insertion)){
  # Insertion plot gets y-axis label
  ggList$Insertion <- ggList$Insertion + 
    ylab("Edit Frequency") + 
    annotate(geom = "text", x = 25, y = .55, 
             label = paste0("Desired edit + insertion: ", round(insedit*100,2), "%"), 
             size = 3, hjust = 1) + 
    annotate(geom = "text", x = 25, y = .48, 
             label = paste0("Insertion only: ", round(insonly*100,2), "%"), 
             size = 3, hjust = 1) +
    labs(fill = "Insertion") +
    scale_fill_manual(values = insertion.colors)
  
  if(!is.null(nrow(cooccur.ins))){
    ggList$Insertion <- ggList$Insertion + annotate(geom = "text", x = cooccur.ins$position, y = -.04, 
                                                    color = cooccur.ins$color, label = "*")
  }
}

if(!is.null(ggList$Deletion)){
  # Deletion plot gets x-axis label
  ggList$Deletion <- ggList$Deletion + 
    xlab("Position relative to nicking site") + 
    annotate(geom = "text", x = 25, y = .55, 
             label = paste0("Desired edit + deletion: ", round(deledit*100,2), "%"), 
             size = 3, hjust = 1) + 
    annotate(geom = "text", x = 25, y = .48, 
             label = paste0("Deletion only: ", round(delonly*100,2), "%"), 
             size = 3, hjust = 1) +
    labs(fill = "Deletion") +
    scale_fill_manual(values = deletion.colors)
  
  if(!is.null(nrow(cooccur.del))){
    ggList$Deletion <- ggList$Deletion + annotate(geom = "text", x = cooccur.del$position, y = -.04, 
                                                  color = cooccur.del$color, label = "*") 
  }
}

### For plots with only two edit categories
# Mismatch + deletion
if(is.null(ggList$Insertion) & !is.null(ggList$Deletion) & !is.null(ggList$Mismatch)){
  ggList$Deletion <- ggList$Deletion + ylab("Edit Frequency")
}
# Mismatch + insertion
if(is.null(ggList$Deletion) & !is.null(ggList$Insertion) & !is.null(ggList$Mismatch)){
  ggList$Insertion <- ggList$Insertion + xlab("Position relative to nicking site")
}
# Insertion + deletion
if(!is.null(ggList$Insertion) & !is.null(ggList$Deletion) & is.null(ggList$Mismatch)){
  ggList$Insertion <- ggList$Insertion + ggtitle(paste0("Editing Outcomes for ", gene.name, 
                                                        " - Edit +",des.pos,":",des.edit, " - ", cl,
                                                        " - ", day, " - ", repl))
}

### For plots with only a single edit category
# Mismatch only
if(is.null(ggList$Insertion) & is.null(ggList$Deletion) & !is.null(ggList$Mismatch)){
  ggList$Mismatch <- ggList$Mismatch + 
                      xlab("Position relative to nicking site") +
                      ylab("Edit Frequency")
}
# Insertion only
if(!is.null(ggList$Insertion) & is.null(ggList$Deletion) & is.null(ggList$Mismatch)){
  ggList$Insertion <- ggList$Insertion + 
                        xlab("Position relative to nicking site") +
                        ggtitle(paste0("Editing Outcomes for ", gene.name, 
                                " - Edit +",des.pos,":",des.edit, " - ", cl,
                                " - ", day, " - ", repl)) 
}
# Deletion only
if(is.null(ggList$Insertion) & !is.null(ggList$Deletion) & is.null(ggList$Mismatch)){
  ggList$Deletion <- ggList$Deletion + 
                      ylab("Edit Frequency") + 
                      ggtitle(paste0("Editing Outcomes for ", gene.name, 
                      " - Edit +",des.pos,":",des.edit, " - ", cl,
                      " - ", day, " - ", repl)) 
}

# Combine plots
edit.plt <- cowplot::plot_grid(plotlist = ggList, ncol = 1, align = 'v')

print("Saving outcome characteristics...")

# Save characterizations
ggsave(paste0(out.path, "/outcome_characterizations_",pegID,".pdf"),
       edit.plt,
       device = "pdf",
       height = 7,
       width = 10)

print("Step 3 - Characterizing outcomes - is complete...")
print(paste0(pegID, " successfully analyzed. Script complete."))