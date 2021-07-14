# Feature calculation - per tRF

library(stringr)
library(stringi)

# Get tRF info
trf_info <- read.csv("trfdb_info.csv")
trf_info <- trf_info[,c(2,4,10)]
trf_info = data.frame(lapply(trf_info, as.character), stringsAsFactors=FALSE)
trf_info$Type <- trimws(trf_info$Type)
trf_info$tRF.Sequence <- trimws(trf_info$tRF.Sequence)

# Get 3' UTR info for all transcripts in file containing gene information
all_genes <- read.csv("grch37_all_genes_curated.csv", row.names = 1)
all_genes <- data.frame(lapply(all_genes, as.character), stringsAsFactors = FALSE)

# First pass with 7mer-m1 seed match to tRF
trf_idx = 1
trf_seed <- substr(trf_info$tRF.Sequence[trf_idx], 1, 7)
# Generate complementary seed
comp_seed <- ""
for(i in 1:nchar(trf_seed)) {
  if(substr(trf_seed, i, i) == "C") {
    comp_seed <- paste(comp_seed, "G", sep = "")
  }
  else if(substr(trf_seed, i, i) == "G"){
    comp_seed <- paste(comp_seed, "C", sep = "")
  }
  else if(substr(trf_seed, i, i) == "T"){
    comp_seed <- paste(comp_seed, "A", sep = "")
  }
  else if(substr(trf_seed, i, i) == "A"){
    comp_seed <- paste(comp_seed, "T", sep = "")
  }
}
# Reverse complementary seed and find which genes contain the string
reverse_comp_seed = stri_reverse(comp_seed)
seed_locations <- which(grepl(reverse_comp_seed, all_genes$sequence))
seed_match <- data.frame(matrix(nrow = length(seed_locations), ncol = 23))
colnames(seed_match) <- c("trfdb_id", "trf_sequence", "tran_id", "tran_ver", "name", "mrna_sequence","chr","start_loc","end_loc","length_utr","trf_binding_loc", "mrna_binding_loc", "binding_energy", "seed", "au_content", "num_paired_pos", "binding_region_length", "longest_consecutive","pos_longest_consecutive", "three_prime_pairs", "seed_end_diff", "phylop_stem", "phylop_flanking")


# Second pass with binding energy for tRF
trf_idx = 1
# Create a FASTA file with tRF and all transcripts (see RNAduplex documentation for input format)
sink("duplex.txt")
for(k in 1:length(all_genes[,1])) {
  cat(paste(">",trf_info$tRF.ID[trf_idx], sep=""))
  cat("\n")
  cat(trf_info$tRF.Sequence[trf_idx])
  cat("\n")
  cat(paste(">", all_genes$tran_id[k], sep=""))
  cat("\n")
  cat(all_genes$sequence[k])
  cat("\n")
}
sink()

# Run RNAduplex via terminal 
terminal_command <- paste("RNAduplex < ", "duplex.txt", " > ", "rnaduplex_out.txt", sep = "")
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}
while(!file.exists("rnaduplex_out.txt")) {
  Sys.sleep(0.1)
}

# Get binding energies for each tRF-transcript duplex
rna_duplex_raw <- read.csv("rnaduplex_out.txt", header = FALSE, sep = "")
rna_duplex_raw <- data.frame(lapply(rna_duplex_raw, as.character), stringsAsFactors=FALSE)
rna_duplex <- data.frame(matrix(nrow = length(rna_duplex_raw[,1])/3, ncol = 5))
colnames(rna_duplex) = c("trf", "tran_id", "trf_binding_loc", "mrna_binding_loc", "binding_energy")
rna_duplex$trf = rna_duplex_raw$V1[seq(1, length(rna_duplex_raw[,1]), 3)]
rna_duplex$trf = substr(rna_duplex$trf, 2, nchar(rna_duplex$trf))
rna_duplex$tran_id = rna_duplex_raw$V1[seq(2, length(rna_duplex_raw[,1]), 3)]
rna_duplex$tran_id = substr(rna_duplex$tran_id, 2, nchar(rna_duplex$tran_id))
rna_duplex$trf_binding_loc = rna_duplex_raw$V2[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$mrna_binding_loc = rna_duplex_raw$V4[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$binding_energy = rna_duplex_raw$V5[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$binding_energy = as.numeric(substr(rna_duplex$binding_energy, 2, nchar(rna_duplex$binding_energy) - 1))

# Select transcripts with appropriate binding energy - see initial processing code for information on obtaining this cutoff
binding_locations <- which(rna_duplex$binding_energy <= trf_info$binding_cutoff[trf_idx])
binding_cutoff <- data.frame(matrix(nrow = length(binding_locations), ncol = 23))
colnames(binding_cutoff) <- c("trfdb_id", "trf_sequence", "tran_id", "tran_ver", "name", "mrna_sequence","chr","start_loc","end_loc","length_utr","trf_binding_loc", "mrna_binding_loc", "binding_energy", "seed", "au_content", "num_paired_pos", "binding_region_length", "longest_consecutive","pos_longest_consecutive", "three_prime_pairs", "seed_end_diff", "phylop_stem", "phylop_flanking")

# Get intersection of seed match and binding energy locations - use it to populate data frame for feature calculation
locations <- seed_locations[which(!is.na(match(seed_locations, binding_locations)))]
feature_profiles <- data.frame(matrix(nrow = length(locations), ncol = 23))
colnames(feature_profiles) <- c("trfdb_id", "trf_sequence", "tran_id", "tran_ver", "name", "mrna_sequence","chr","start_loc","end_loc","length_utr","trf_binding_loc", "mrna_binding_loc", "binding_energy", "seed", "au_content", "num_paired_pos", "binding_region_length", "longest_consecutive","pos_longest_consecutive", "three_prime_pairs", "seed_end_diff", "phylop_stem", "phylop_flanking")

# Get tRF ID
feature_profiles$trfdb_id = trf_info$tRF.ID[trf_idx]

# Get tRF Sequence
feature_profiles$trf_sequence = trf_info$tRF.Sequence[trf_idx]

# Get transcript ID with and without version, gene name, and 3' UTR mRNA sequence
feature_profiles$tran_id = all_genes$tran_id[locations]
feature_profiles$tran_ver = all_genes$tran_ver[locations]
feature_profiles$name = all_genes$name[locations]
feature_profiles$mrna_sequence = all_genes$sequence[locations]

# Get chromosomal coordinates and length of 3' UTR

feature_profiles$chr = all_genes$chr[locations]
feature_profiles$start_loc = all_genes$utr3_start_bp[locations]
feature_profiles$end_loc = all_genes$utr3_end_bp[locations]
feature_profiles$length_utr = nchar(all_genes$sequence[locations])

# Get binding locations on tRF and mRNA and binding energy
sink("duplex.txt")
for(k in 1:length(feature_profiles[,1])) {
  cat(paste(">",feature_profiles$trfdb_id[k], sep=""))
  cat("\n")
  cat(feature_profiles$trf_sequence[k])
  cat("\n")
  cat(paste(">", feature_profiles$tran_id[k], sep=""))
  cat("\n")
  cat(feature_profiles$mrna_sequence[k])
  cat("\n")
}
sink()

terminal_command <- paste("RNAduplex < ", "duplex.txt", " > ", "rnaduplex_out.txt", sep = "")
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}
while(!file.exists("rnaduplex_out.txt")) {
  Sys.sleep(0.1)
}

rna_duplex_raw <- read.csv("rnaduplex_out.txt", header = FALSE, sep = "")
rna_duplex_raw <- data.frame(lapply(rna_duplex_raw, as.character), stringsAsFactors=FALSE)
rna_duplex <- data.frame(matrix(nrow = length(rna_duplex_raw[,1])/3, ncol = 5))
colnames(rna_duplex) = c("trf", "tran_id", "trf_binding_loc", "mrna_binding_loc", "binding_energy")
rna_duplex$trf = rna_duplex_raw$V1[seq(1, length(rna_duplex_raw[,1]), 3)]
rna_duplex$trf = substr(rna_duplex$trf, 2, nchar(rna_duplex$trf))
rna_duplex$tran_id = rna_duplex_raw$V1[seq(2, length(rna_duplex_raw[,1]), 3)]
rna_duplex$tran_id = substr(rna_duplex$tran_id, 2, nchar(rna_duplex$tran_id))
rna_duplex$trf_binding_loc = rna_duplex_raw$V2[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$mrna_binding_loc = rna_duplex_raw$V4[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$binding_energy = rna_duplex_raw$V5[seq(3, length(rna_duplex_raw[,1]), 3)]
rna_duplex$binding_energy = as.numeric(substr(rna_duplex$binding_energy, 2, nchar(rna_duplex$binding_energy) - 1))

locations = match(rna_duplex$tran_id, feature_profiles$tran_id)
feature_profiles$trf_binding_loc = rna_duplex$trf_binding_loc[locations]
feature_profiles$mrna_binding_loc = rna_duplex$mrna_binding_loc[locations]
feature_profiles$binding_energy = rna_duplex$binding_energy[locations]

# Determine if seed exists on binding site
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  site_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  site_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  binding_site = substr(feature_profiles$mrna_sequence[k], site_start, site_end)
  if(grepl(reverse_comp_seed, binding_site, fixed = TRUE)) {
    feature_profiles$seed[k] <- 1
  }
  else {
    feature_profiles$seed[k] <- 0
  }
}

# Calculate AU content around binding site (30 nts upstream and downstream)
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  site_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  site_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  if(site_start <= 30) {
    downstream <- substr(feature_profiles$mrna_sequence[k], 1, site_start - 1)
  } else {
    downstream <- substr(feature_profiles$mrna_sequence[k], site_start - 30, site_start - 1)
  }
  if(site_end >= (nchar(feature_profiles$mrna_sequence[k]) - 30)) {
    upstream <- substr(feature_profiles$mrna_sequence[k], site_end + 1, nchar(feature_profiles$mrna_sequence[k]))
  } else {
    upstream <- substr(feature_profiles$mrna_sequence[k], site_end + 1, site_end + 30)
  }
  downstream <- stri_reverse(downstream)
  upscore <- 0
  downscore <- 0
  # Formula is a weighted sum with distance from binding site
  for(i in 1:30) {
    if(substr(upstream, i, i) == "A" || substr(upstream, i, i) == "T") {
      upscore <- upscore + (1/(i + 1))
    }
    if(substr(downstream, i, i) == "A" || substr(downstream, i, i) == "T") {
      downscore <- downscore + (1/(i + 1))
    }
  }
  feature_profiles$au_content[k] <- mean(c(upscore, downscore))
}

# Calculate number of paired positions on binding site

# Function to quickly return canonical nucleotide pair
isPaired <- function(nt) {
  if(nt == "A") {
    return("T")
  }
  else if(nt == "T") {
    return("A")
  }
  else if(nt == "G") {
    return("C")
  }
  else if(nt == "C") {
    return("G")
  }
  else {
    return("")
  }
}

for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  trf_start <- as.integer(substr(feature_profiles$trf_binding_loc[k], 1, str_locate(feature_profiles$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(feature_profiles$trf_binding_loc[k], str_locate(feature_profiles$trf_binding_loc[k], ",") + 1, nchar(feature_profiles$trf_binding_loc[k])))
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  count_paired <- 0
  for(i in 1:min(mrna_end-mrna_start + 1, trf_end-trf_start + 1)){
    if(isPaired(substr(feature_profiles$trf_sequence[k], trf_start + i - 1, trf_start + i - 1)) == substr(feature_profiles$mrna_sequence[k], mrna_start + i - 1, mrna_start + i - 1)) {
      count_paired <- count_paired + 1
    }
  }
  feature_profiles$num_paired_pos[k] <- count_paired
}


# Calculate length of binding region
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  feature_profiles$binding_region_length[k] <- mrna_end - mrna_start + 1
}

# Calculate the longest consecutive pairing and position of it on binding site
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  trf_start <- as.integer(substr(feature_profiles$trf_binding_loc[k], 1, str_locate(feature_profiles$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(feature_profiles$trf_binding_loc[k], str_locate(feature_profiles$trf_binding_loc[k], ",") + 1, nchar(feature_profiles$trf_binding_loc[k])))
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  paired_info <- c()
  for(i in 1:min(mrna_end-mrna_start + 1, trf_end-trf_start + 1)){
    if(isPaired(substr(feature_profiles$trf_sequence[k], trf_start + i - 1, trf_start + i - 1)) == substr(feature_profiles$mrna_sequence[k], mrna_start + i - 1, mrna_start + i - 1)) {
      paired_info <- c(paired_info, 1)
    }
    else {
      paired_info <- c(paired_info, 0)
    }
  }
  count_consec <- 0
  pos_consec <- 0
  for(i in 1:length(rle(paired_info)$values)) {
    if(rle(paired_info)$values[i] == 1) {
      if(rle(paired_info)$lengths[i] > count_consec) {
        count_consec <- rle(paired_info)$lengths[i]
        pos_consec <- sum(rle(paired_info)$lengths[1:i])
      }
    }
  }
  feature_profiles$longest_consecutive[k] <- count_consec
  feature_profiles$pos_longest_consecutive[k] <- pos_consec + trf_start - 1
}

# Calculate number of pairs in 3' region of tRF (last 7 nts)
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  trf_start <- as.integer(substr(feature_profiles$trf_binding_loc[k], 1, str_locate(feature_profiles$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(feature_profiles$trf_binding_loc[k], str_locate(feature_profiles$trf_binding_loc[k], ",") + 1, nchar(feature_profiles$trf_binding_loc[k])))
  three_prime_start <- trf_end - 7
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1)) + ((trf_end - trf_start) - 7)
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  count_paired <- 0
  for(i in 1:min(mrna_end-mrna_start + 1, 7)){
    if(isPaired(substr(feature_profiles$trf_sequence[k], trf_start + i - 1, trf_start + i - 1)) == substr(feature_profiles$mrna_sequence[k], mrna_start + i - 1, mrna_start + i - 1)) {
      count_paired <- count_paired + 1
    }
  }
  feature_profiles$three_prime_pairs[k] <- count_paired
}

# Difference in pairs between the seed region and 3' region of tRF
for(k in 1:length(feature_profiles[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  if(feature_profiles$seed[k] == 1) {
    feature_profiles$seed_end_diff[k] <- abs(feature_profiles$three_prime_pairs[k] - 6)
  } else {
    trf_seed <- substr(feature_profiles$trf_sequence[k], 2, 7)
    comp_seed <- ""
    for(i in 1:nchar(trf_seed)) {
      comp_seed <- paste(comp_seed, isPaired(substr(trf_seed, i, i)), sep = "")
    }
    mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
    mrna_comp_seed <- substr(feature_profiles$mrna_sequence[k], mrna_start + 1, mrna_start + 6)
    count_seed <- 0
    for(i in 1:nchar(trf_seed)) {
      if(substr(comp_seed, i, i) == substr(mrna_comp_seed, i, i)) {
        count_seed <- count_seed + 1
      }
    }
    feature_profiles$seed_end_diff[k] <- abs(feature_profiles$three_prime_pairs[k] - count_seed)
  }
}

# Calculate phyloP scores - stem region (binding site) and flanking region ((40 nts upstream and downstream of stem)

# Stem region - first generate unsorted bed files
#   Need chromosome, 3' UTR coordinates, and transcript ID

sink("phylop_stem_unsorted.bed")
for(k in 1:length(feature_profiles[,1])) {
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  cat(paste("chr", as.character(feature_profiles$chr[k]), sep = ""))
  cat("\t")
  if(!is.na(suppressWarnings(as.integer(feature_profiles$start_loc[k])))) {
    cat(as.character(as.integer(feature_profiles$start_loc[k]) + mrna_start - 1))
    cat("\t")
    cat(as.character(as.integer(feature_profiles$end_loc[k]) - (feature_profiles$length_utr[k] - mrna_end)))
    cat("\t")
    cat(feature_profiles$tran_id[k])
    cat("\n")
  } else {
    # Get piece-wise 3' UTR chromosomal coordinates
    x = unlist(str_locate_all(feature_profiles$start_loc[k], ";"))
    x = x[1:(length(x)/2)]
    y = unlist(str_locate_all(feature_profiles$end_loc[k], ";"))
    y = y[1:(length(y)/2)]
    all_start_locs = c(as.integer(substr(feature_profiles$start_loc[k], 1, x[1] - 1)))
    if(length(x)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[length(x)] + 1, nchar(feature_profiles$start_loc[k]))))
    all_end_locs = c(as.integer(substr(feature_profiles$end_loc[k], 1, x[1] - 1)))
    if(length(y)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[length(x)] + 1, nchar(feature_profiles$end_loc[k]))))
    all_start_locs = sort(all_start_locs)
    all_end_locs = sort(all_end_locs)
    differences = all_end_locs - all_start_locs + 1
    for(i in 2:length(differences)) {
      differences[i] = sum(differences[(i-1):i])
    }
    if(sum(mrna_start > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1)))
      cat("\t")
    } else{ 
      which_start_loc = max(which(mrna_start > differences))
      cat(as.character(all_start_locs[which_start_loc + 1] + (mrna_start-differences[which_start_loc] - 1)))
      cat("\t")
    }
    if(sum(mrna_end > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1) + (mrna_end - mrna_start)))
      cat("\t")
    } else {
      which_end_loc = max(which(mrna_end > differences))
      cat(as.character(all_start_locs[which_end_loc + 1] + (mrna_end - differences[which_end_loc] - 1)))
      cat("\t")  
    }
    cat(feature_profiles$tran_id[k])
    cat("\n")
  }
}
sink()

# Sort bed file for mapping
terminal_command <- "sort-bed phylop_stem_unsorted.bed > phylop_stem_sorted.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

# Break up sorted bed file by 1-base sections to get score for each individual nucleotide in region
terminal_command <- "awk \' { regionChromosome = $1; regionStart = $2; regionStop = $3; regionID = $4; baseIdx = 0; for (baseStart = regionStart; baseStart < regionStop; baseStart++) { baseStop = baseStart + 1; print regionChromosome\"\\t\"baseStart\"\\t\"baseStop\"\\t\"regionID\"-\"baseIdx; baseIdx++; } }\' phylop_stem_sorted.bed > phylop_stem_perBase.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

# Flanking region - same process

sink("phylop_downstream_unsorted.bed")
for(k in 1:length(feature_profiles[,1])) {
  mrna_start <- as.integer(substr(feature_profiles$mrna_binding_loc[k], 1, str_locate(feature_profiles$mrna_binding_loc[k], ",") - 1))
  cat(paste("chr", as.character(feature_profiles$chr[k]), sep = ""))
  cat("\t")
  if(!is.na(suppressWarnings(as.integer(feature_profiles$start_loc[k])))) {
    cat(as.character(as.integer(feature_profiles$start_loc[k]) + mrna_start - 1 - 40))
    cat("\t")
    cat(as.character(as.integer(feature_profiles$start_loc[k]) + mrna_start - 1 - 1))
    cat("\t")
    cat(feature_profiles$tran_id[k])
    cat("\n")
  } else {
    x = unlist(str_locate_all(feature_profiles$start_loc[k], ";"))
    x = x[1:(length(x)/2)]
    y = unlist(str_locate_all(feature_profiles$end_loc[k], ";"))
    y = y[1:(length(y)/2)]
    all_start_locs = c(as.integer(substr(feature_profiles$start_loc[k], 1, x[1] - 1)))
    if(length(x)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[length(x)] + 1, nchar(feature_profiles$start_loc[k]))))
    all_end_locs = c(as.integer(substr(feature_profiles$end_loc[k], 1, x[1] - 1)))
    if(length(y)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[length(x)] + 1, nchar(feature_profiles$end_loc[k]))))
    all_start_locs = sort(all_start_locs)
    all_end_locs = sort(all_end_locs)
    differences = all_end_locs - all_start_locs + 1
    for(i in 2:length(differences)) {
      differences[i] = sum(differences[(i-1):i])
    }
    if(sum(mrna_start > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1) - 40))
      cat("\t")
      cat(as.character(all_start_locs[1] + (mrna_start - 1) - 1))
      cat("\t")
    } else{ 
      which_start_loc = max(which(mrna_start > differences))
      cat(as.character(all_start_locs[which_start_loc + 1] + (mrna_start-differences[which_start_loc] - 1) - 40))
      cat("\t")
      cat(as.character(all_start_locs[which_start_loc + 1] + (mrna_start-differences[which_start_loc] - 1) - 1))
      cat("\t")
    }
    cat(feature_profiles$tran_id[k])
    cat("\n")
  }
}
sink()

sink("phylop_upstream_unsorted.bed")
for(k in 1:length(feature_profiles[,1])) {
  mrna_end <- as.integer(substr(feature_profiles$mrna_binding_loc[k], str_locate(feature_profiles$mrna_binding_loc[k], ",") + 1, nchar(feature_profiles$mrna_binding_loc[k])))
  cat(paste("chr", as.character(feature_profiles$chr[k]), sep = ""))
  cat("\t")
  if(!is.na(suppressWarnings(as.integer(feature_profiles$start_loc[k])))) {
    cat(as.character(as.integer(feature_profiles$end_loc[k]) - (feature_profiles$length_utr[k] - mrna_end) + 1))
    cat("\t")
    cat(as.character(as.integer(feature_profiles$end_loc[k]) - (feature_profiles$length_utr[k] - mrna_end) + 40))
    cat("\t")
    cat(feature_profiles$tran_id[k])
    cat("\n")
  } else {
    x = unlist(str_locate_all(feature_profiles$start_loc[k], ";"))
    x = x[1:(length(x)/2)]
    y = unlist(str_locate_all(feature_profiles$end_loc[k], ";"))
    y = y[1:(length(y)/2)]
    all_start_locs = c(as.integer(substr(feature_profiles$start_loc[k], 1, x[1] - 1)))
    if(length(x)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_start_locs = c(all_start_locs, as.integer(substr(feature_profiles$start_loc[k], x[length(x)] + 1, nchar(feature_profiles$start_loc[k]))))
    all_end_locs = c(as.integer(substr(feature_profiles$end_loc[k], 1, x[1] - 1)))
    if(length(y)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_end_locs = c(all_end_locs, as.integer(substr(feature_profiles$end_loc[k], x[length(x)] + 1, nchar(feature_profiles$end_loc[k]))))
    all_start_locs = sort(all_start_locs)
    all_end_locs = sort(all_end_locs)
    differences = all_end_locs - all_start_locs + 1
    for(i in 2:length(differences)) {
      differences[i] = sum(differences[(i-1):i])
    }
    if(sum(mrna_end > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1) + (mrna_end - mrna_start) + 1))
      cat("\t")
      cat(as.character(all_start_locs[1] + (mrna_start - 1) + (mrna_end - mrna_start) + 40))
      cat("\t")
    } else {
      which_end_loc = max(which(mrna_end > differences))
      cat(as.character(all_start_locs[which_end_loc + 1] + (mrna_end - differences[which_end_loc] - 1) + 1))
      cat("\t")  
      cat(as.character(all_start_locs[which_end_loc + 1] + (mrna_end - differences[which_end_loc] - 1) + 40))
      cat("\t")
    }
    cat(feature_profiles$tran_id[k])
    cat("\n")
  }
}
sink()

terminal_command <- "sort-bed phylop_downstream_unsorted.bed > phylop_downstream_sorted.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

terminal_command <- "sort-bed phylop_upstream_unsorted.bed > phylop_upstream_sorted.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

terminal_command <- "awk \' { regionChromosome = $1; regionStart = $2; regionStop = $3; regionID = $4; baseIdx = 0; for (baseStart = regionStart; baseStart < regionStop; baseStart++) { baseStop = baseStart + 1; print regionChromosome\"\\t\"baseStart\"\\t\"baseStop\"\\t\"regionID\"-\"baseIdx; baseIdx++; } }\' phylop_downstream_sorted.bed > phylop_downstream_perBase.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

terminal_command <- "awk \' { regionChromosome = $1; regionStart = $2; regionStop = $3; regionID = $4; baseIdx = 0; for (baseStart = regionStart; baseStart < regionStop; baseStart++) { baseStop = baseStart + 1; print regionChromosome\"\\t\"baseStart\"\\t\"baseStop\"\\t\"regionID\"-\"baseIdx; baseIdx++; } }\' phylop_upstream_sorted.bed > phylop_upstream_perBase.bed"
termID <- rstudioapi::terminalExecute(terminal_command)
while(is.null(rstudioapi::terminalExitCode(termID))) {
  Sys.sleep(0.1)
}
while(rstudioapi::terminalExitCode(termID) != 0) {
  rstudioapi::terminalKill(termID)
  termID <- rstudioapi::terminalExecute(terminal_command)
  while(is.null(rstudioapi::terminalExitCode(termID))) {
    Sys.sleep(0.1)
  }
}

# Use bedmap to get scores using reference files (see bedops documentation)

# Get results from files (broken up by chromosome)

bed_files = list.files()
for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], header = FALSE, sep = "\t")
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  chr[,c(2,3)] <- data.frame(lapply(chr[,c(2,3)], as.integer), stringsAsFactors=FALSE)
  chr$score <- 0
  colnames(chr) <- c("chr", "start_loc", "end_loc", "gene_base", "score")
  chr$score <- as.double(substr(chr$gene_base, str_locate(chr$gene_base, "\\|") + 1, nchar(chr$gene_base)))
  chr$gene_base <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "\\|") - 1)
  write.csv(chr, bed_files[k])
}
for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], row.names = 1)
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  genes <- unique(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  # Fill in scores for any duplicate coordinates
  chrloc_pairs <- paste(chr$chr, chr$start_loc)
  duplicate_matches <- match(chrloc_pairs, unique(chrloc_pairs))
  for(i in 1:length(duplicate_matches)) {
    if(i%%1000==0){print(i)}
    x <- chr$score[which(duplicate_matches == duplicate_matches[i])]
    best_value <- x[!is.na(x)][1]
    x <- rep(best_value, length(x))
    chr$score[which(duplicate_matches == duplicate_matches[i])] = x
  }
  # Get score for relevant region
  chr$score[which(is.na(chr$score))] = 0
  rle <- rle(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  condensed_chr <- data.frame(matrix(nrow = length(rle$values), ncol = 0))
  condensed_chr$gene_id <- genes
  condensed_chr$avg_score <- NA
  gene_nobp <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1)
  for(i in 1:length(genes)) {
    condensed_chr$avg_score[i] <- mean(chr$score[which(genes[i] == gene_nobp)])
  }
  
  loc_genes <- match(genes, feature_profiles$tran_id)
  for(i in 1:length(loc_genes)) {
    if(i%%1000==0){print(i)}
    if(is.na(loc_genes[i])) {
      next
    } else {
      feature_profiles$phylop_stem[loc_genes[i]] <- condensed_chr$avg_score[i]
    }
  }
}

feature_profiles$phylop_upstream = NA
bed_files = list.files()
for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], header = FALSE, sep = "\t")
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  chr[,c(2,3)] <- data.frame(lapply(chr[,c(2,3)], as.integer), stringsAsFactors=FALSE)
  chr$score <- 0
  colnames(chr) <- c("chr", "start_loc", "end_loc", "gene_base", "score")
  chr$score <- as.double(substr(chr$gene_base, str_locate(chr$gene_base, "\\|") + 1, nchar(chr$gene_base)))
  chr$gene_base <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "\\|") - 1)
  write.csv(chr, bed_files[k])
}

for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], row.names = 1)
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  genes <- unique(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  
  chrloc_pairs <- paste(chr$chr, chr$start_loc)
  duplicate_matches <- match(chrloc_pairs, unique(chrloc_pairs))
  for(i in 1:length(duplicate_matches)) {
    if(i%%1000==0){print(i)}
    x <- chr$score[which(duplicate_matches == duplicate_matches[i])]
    best_value <- x[!is.na(x)][1]
    x <- rep(best_value, length(x))
    chr$score[which(duplicate_matches == duplicate_matches[i])] = x
  }
  chr$score[which(is.na(chr$score))] = 0
  rle <- rle(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  condensed_chr <- data.frame(matrix(nrow = length(rle$values), ncol = 0))
  condensed_chr$gene_id <- genes
  condensed_chr$avg_score <- NA
  gene_nobp <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1)
  for(i in 1:length(genes)) {
    condensed_chr$avg_score[i] <- mean(chr$score[which(genes[i] == gene_nobp)])
  }
  
  loc_genes <- match(genes, feature_profiles$tran_id)
  for(i in 1:length(loc_genes)) {
    if(i%%1000==0){print(i)}
    if(is.na(loc_genes[i])) {
      next
    } else {
      feature_profiles$phylop_upstream[loc_genes[i]] <- condensed_chr$avg_score[i]
    }
  }
}

feature_profiles$phylop_downstream = NA
bed_files = list.files()
for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], header = FALSE, sep = "\t")
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  chr[,c(2,3)] <- data.frame(lapply(chr[,c(2,3)], as.integer), stringsAsFactors=FALSE)
  chr$score <- 0
  colnames(chr) <- c("chr", "start_loc", "end_loc", "gene_base", "score")
  chr$score <- as.double(substr(chr$gene_base, str_locate(chr$gene_base, "\\|") + 1, nchar(chr$gene_base)))
  chr$gene_base <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "\\|") - 1)
  write.csv(chr, bed_files[k])
}
for(k in 1:length(bed_files)) {
  print(k)
  chr <- read.csv(bed_files[k], row.names = 1)
  chr[,c(1,4)] <- data.frame(lapply(chr[,c(1,4)], as.character), stringsAsFactors=FALSE)
  genes <- unique(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  
  chrloc_pairs <- paste(chr$chr, chr$start_loc)
  duplicate_matches <- match(chrloc_pairs, unique(chrloc_pairs))
  for(i in 1:length(duplicate_matches)) {
    if(i%%1000==0){print(i)}
    x <- chr$score[which(duplicate_matches == duplicate_matches[i])]
    best_value <- x[!is.na(x)][1]
    x <- rep(best_value, length(x))
    chr$score[which(duplicate_matches == duplicate_matches[i])] = x
  }
  chr$score[which(is.na(chr$score))] = 0
  rle <- rle(substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1))
  condensed_chr <- data.frame(matrix(nrow = length(rle$values), ncol = 0))
  condensed_chr$gene_id <- genes
  condensed_chr$avg_score <- NA
  gene_nobp <- substr(chr$gene_base, 1, str_locate(chr$gene_base, "-")[,1] - 1)
  for(i in 1:length(genes)) {
    condensed_chr$avg_score[i] <- mean(chr$score[which(genes[i] == gene_nobp)])
  }
  
  loc_genes <- match(genes, feature_profiles$tran_id)
  for(i in 1:length(loc_genes)) {
    if(i%%1000==0){print(i)}
    if(is.na(loc_genes[i])) {
      next
    } else {
      feature_profiles$phylop_downstream[loc_genes[i]] <- condensed_chr$avg_score[i]
    }
  }
}

# Take mean of upstream and downstream phyloP scores to get flanking score
for(k in 1:length(feature_profiles[,1])) {
  if(k%%1000==0){print(k)}
  feature_profiles$phylop_flanking[k] = mean(c(feature_profiles$phylop_upstream[k], feature_profiles$phylop_downstream[k]))
}

# Remove individual upstream and downstream scores
feature_profiles = feature_profiles[,c(-24,-25)]