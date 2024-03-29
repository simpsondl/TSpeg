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
                               help = "Maximum cycle number where a barcode can end  [default = %default]"),
                   make_option(c("-b","--barcode"), type = "character",
                               default = NULL, action = "store",
                               help = "Expected barcode for pegRNA"),
                   make_option(c("-s","--bcstartpos"), type = "integer",
                               default = 54, action = "store",
                               help = "Expected cycle number where barcode begins in read  [default = %default]"),
                   make_option(c("-m","--maxdist"), type = "integer",
                               default = 2, action = "store",
                               help = "Max mismatches allowed in barcode before filtering [default = %default]"),
                   make_option(c("-w","--targetregion"), type = "character",
                               default = NULL, action = "store",
                               help = "Expected target region for pegRNA"),
                   make_option(c("-n","--trorientation"), type = "character",
                               default = "f", action = "store",
                               help = "Orientation of target site in read compared to file (f or r) [default = %default]"),
                   make_option(c("-t","--trstartpos"), type = "integer",
                               default = 7, action = "store",
                               help = "Expected cycle number where target region begins in the read  [default = %default]"),
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

# Check if output directory exists and create if not
if(!dir.exists(out.path)){
  dir.create(out.path)
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
bcs.sum$bc <- rc(str_sub(as.character(sread(reads)), bcs.sum$start,bcs.sum$end))
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
if(orient == "f"){
  wtr.matches <- as.data.frame(as(vmatchPattern2(rc(exp.wtr), 
                                                 sread(reads.filt),
                                                 with.indels = TRUE,
                                                 max.mismatch = floor(.4*nchar(exp.wtr))), #at least 60% seq identity 
                                                 "CompressedIRangesList"))  
} else {
  wtr.matches <- as.data.frame(as(vmatchPattern2(exp.wtr, 
                                                 sread(reads.filt),
                                                 with.indels = TRUE,
                                                 max.mismatch = floor(.4*nchar(exp.wtr))), #at least 60% seq identity 
                                                 "CompressedIRangesList"))
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
if(orient == "f"){
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

wtrfreqfp <- paste0(out.path, "/corrected_outcome_frequencies_", pegID, ".txt")
print(getwd())
print(wtrfreqfp)
corrected.wtr <- read.delim(wtrfreqfp)

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
      insertions$edit[j] <- paste("+",str_sub(as.character(pattern(tmp.align)), insertions$start[j], insertions$end[j]), 
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
      deletions$edit[j] <- paste("-",str_sub(as.character(subject(tmp.align)), deletions$start[j], deletions$end[j]), 
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
  uetotal <- 0
  wttotal <- 0
  eonly <- 0
  for(i in 1:nrow(corrected.wtr)){
    if(corrected.wtr$wtr[i] == exp.wtr){
      wttotal <- wttotal + corrected.wtr$freq[i]
    } else {
      uetotal <- uetotal + corrected.wtr$freq[i]
    }
  }
  
  print("ERROR: No editing detected...")
  print("Saving statistics to file...")
  
  write.table(data.frame(pegID = pegID,
                         CellLine = cl,
                         Day = day,
                         Replicate = repl,
                         TotalReads = sum(corrected.wtr$N),
                         Freq.TotalWildType = wttotal,
                         Freq.TotalOffTarget = uetotal,
                         Freq.EditOnly = eonly,
			 Freq.Sum = wttotal + uetotal + eonly),
              paste0(out.path,"/summarized_outcomes_",pegID,".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  
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

outcome.char$id <- as.character(outcome.char$id)

missing <- setdiff(corrected.wtr$id, outcome.char$id)
if("Outcome_1" %in% missing){
  if(corrected.wtr$wtr[corrected.wtr$id == "Outcome_1"] == exp.wtr){
    outcome.char <- rbind(data.frame(id = "Outcome_1", HasEdit = 0, HasOtherMM = 0,
                                     HasIns = FALSE, HasDel = FALSE, freq = corrected.wtr$freq[1],
                                     EditOnly = FALSE, EditOtherMM = FALSE, OtherOnly = FALSE,
                                     stringsAsFactors = FALSE),
                          as.data.frame(outcome.char))  
  } else {
    outcome.char <- rbind(data.frame(id = "Outcome_1", HasEdit = 0, HasOtherMM = 0,
                                     HasIns = FALSE, HasDel = TRUE, freq = corrected.wtr$freq[i],
                                     EditOnly = FALSE, EditOtherMM = FALSE, OtherOnly = TRUE),
                          as.data.frame(outcome.char))
  }
  
}

still.missing <- as.numeric(gsub("Outcome_", "", setdiff(missing, "Outcome_1")))

for(i in still.missing){
  outcome.char <- rbind(data.frame(id = paste0("Outcome_", i), HasEdit = 0, HasOtherMM = 0,
                                   HasIns = FALSE, HasDel = TRUE, freq = corrected.wtr$freq[i],
                                   EditOnly = FALSE, EditOtherMM = FALSE, OtherOnly = TRUE),
                        as.data.frame(outcome.char))
}

outcome.char <- outcome.char[order(as.numeric(gsub("Outcome_","", outcome.char$id))),]

# Calculate frequencies of different editing event combinations
### Edit only
eonly <- sum(outcome.char$freq[outcome.char$EditOnly])

### Frequency of unintended edits alone
uetotal <- sum(outcome.char$freq[(outcome.char$HasOtherMM > 0 |
                                    outcome.char$HasIns |
                                    outcome.char$HasDel)])
### Frequency of WT (no edit)
wttotal <- sum(outcome.char$freq[outcome.char$HasEdit == 0 & 
                                   outcome.char$HasOtherMM == 0 &
                                   !outcome.char$HasIns &
                                   !outcome.char$HasDel])

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
write.table(data.frame(pegID = pegID,
                       CellLine = cl,
                       Day = day,
                       Replicate = repl,
                       TotalReads = sum(corrected.wtr$N),
                       Freq.TotalWildType = wttotal,
                       Freq.TotalOffTarget = uetotal,
                       Freq.EditOnly = eonly,
		       Freq.Sum = wttotal + uetotal + eonly),
            paste0(out.path,"/summarized_outcomes_",pegID,".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
