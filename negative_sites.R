# Generate negative corresponding sites for random forest feature set for training

# Retrieve binding sites
#   "rf_features" is a data frame containing all ground-truth positive sites with features 
all_binding_sites <- data.frame(matrix(, nrow = length(rf_features[,1]), ncol = 1))
for(k in 1:length(rf_features[,1])) {
  mrna_start <- as.integer(substr(rf_features$mrna_binding_loc[k], 1, str_locate(rf_features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(rf_features$mrna_binding_loc[k], str_locate(rf_features$mrna_binding_loc[k], ",") + 1, nchar(rf_features$mrna_binding_loc[k])))
  all_binding_sites[k, 1] <- substr(rf_features$fixed_mrna_transcript[k], mrna_start, mrna_end)
}

# Generate corresponding negative site for each ground-truth positive site.
for(k in 1:length(rf_features[,1])) {
  # Obtain location of binding sites for tRF-gene pair 
  trf_start <- as.integer(substr(rf_features$trf_binding_loc[k], 1, str_locate(rf_features$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(rf_features$trf_binding_loc[k], str_locate(rf_features$trf_binding_loc[k], ",") + 1, nchar(rf_features$trf_binding_loc[k])))
  mrna_start <- as.integer(substr(rf_features$mrna_binding_loc[k], 1, str_locate(rf_features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(rf_features$mrna_binding_loc[k], str_locate(rf_features$mrna_binding_loc[k], ",") + 1, nchar(rf_features$mrna_binding_loc[k])))
  # Obtain mRNA transcript split up into binding site and areas before and after it
  before_site <- substr(rf_features$fixed_mrna_transcript[k], 1, mrna_start - 1)
  after_site <- substr(rf_features$fixed_mrna_transcript[k], mrna_end + 1, nchar(rf_features$fixed_mrna_transcript[k]))
  binding_site <- substr(rf_features$fixed_mrna_transcript[k], mrna_start, mrna_end)
  # Set sliding window size
  window_size <- mrna_end - mrna_start
  i <- 1
  # Obtain possible candidate sites for a given trf-gene pair
  candidate_sites <- data.frame(matrix(, nrow = 1, ncol = 4))
  colnames(candidate_sites) <- c("trf_id", "fixed_trf_sequence", "mrna_id", "mrna_neg_sites")
  candidate_sites <- data.frame(lapply(candidate_sites, as.character), stringsAsFactors=FALSE)
  cand_site_idx <- 1
  
  # Look in the 3' UTR before the binding site to find candidate negative sites and filter them
  while((i + window_size) < mrna_start) {
    candidate_neg <- substr(before_site, i, i + window_size)
    # Filter by CG dinucleotide frequency
    candidate_cg <- length(as.data.frame(str_locate_all(candidate_neg, "CG"))[,1])
    binding_cg <- length(as.data.frame(str_locate_all(binding_site, "CG"))[,1])
    if(candidate_cg == binding_cg) {
      # Filter by G nucleotide frequency
      candidate_g <- c()
      binding_g <- c()
      for(j in 1:window_size) {
        if(substr(candidate_neg, j, j) == "G") {
          candidate_g <- c(candidate_g, 1)
        }
        else {
          candidate_g <- c(candidate_g, 0)
        }
        if(substr(binding_site, j, j) == "G") {
          binding_g <- c(binding_g, 1)
        }
        else {
          binding_g <- c(binding_g, 0)
        }
      }
      # Perform a statistical test to determine if the candidate and binding site have similar G nucleotide frequency
      if(suppressWarnings(wilcox.test(candidate_g, binding_g, paired = FALSE)$p.value) >= 0.05 || is.na(suppressWarnings(wilcox.test(candidate_g, binding_g, paired = FALSE)$p.value))) {
        check_binding <- 0
        # If found, add site to list of candidate sites
        if(check_binding == 0) {
          candidate_sites[nrow(candidate_sites)+1,] <- NA
          candidate_sites$trf_id[cand_site_idx] <- rf_features$trf_id[k]
          candidate_sites$fixed_trf_sequence[cand_site_idx] <- rf_features$fixed_trf_sequence[k]
          candidate_sites$mrna_id[cand_site_idx] <- rf_features$mrna_id[k]
          candidate_sites$mrna_neg_sites[cand_site_idx] <- candidate_neg
          cand_site_idx <- cand_site_idx + 1
        }
      }
    }
    i <- i + 1
  }
  i <- 1
  # Look in the 3' UTR after the binding site to find candidate negative sites and filter them
  while((i + window_size) <= nchar(rf_features$fixed_mrna_transcript[k]) - mrna_end) {
    candidate_neg <- substr(after_site, i, i + window_size)
    # Filter by CG dinucleotide frequency
    candidate_cg <- length(as.data.frame(str_locate_all(candidate_neg, "CG"))[,1])
    binding_cg <- length(as.data.frame(str_locate_all(binding_site, "CG"))[,1])
    if(candidate_cg == binding_cg) {
      # Filter by G nucleotide frequency
      candidate_g <- c()
      binding_g <- c()
      for(j in 1:window_size) {
        if(substr(candidate_neg, j, j) == "G") {
          candidate_g <- c(candidate_g, 1)
        }
        else {
          candidate_g <- c(candidate_g, 0)
        }
        if(substr(binding_site, j, j) == "G") {
          binding_g <- c(binding_g, 1)
        }
        else {
          binding_g <- c(binding_g, 0)
        }
      }
      # Perform a statistical test to determine if the candidate and binding site have similar G nucleotide frequency
      if(suppressWarnings(wilcox.test(candidate_g, binding_g, paired = FALSE)$p.value) >= 0.05 || is.na(suppressWarnings(wilcox.test(candidate_g, binding_g, paired = FALSE)$p.value))) {
        check_binding <- 0
        # If found, add site to list of candidate sites
        if(check_binding == 0) {
          candidate_sites[nrow(candidate_sites)+1,] <- NA
          candidate_sites$trf_id[cand_site_idx] <- rf_features$trf_id[k]
          candidate_sites$fixed_trf_sequence[cand_site_idx] <- rf_features$fixed_trf_sequence[k]
          candidate_sites$mrna_id[cand_site_idx] <- rf_features$mrna_id[k]
          candidate_sites$mrna_neg_sites[cand_site_idx] <- candidate_neg
          cand_site_idx <- cand_site_idx + 1
        }
      }
    }
    i <- i + 1
  }
  # Filter all candidate sites by binding energy
  sink(duplex_filename)
  for(k in 1:length(candidate_sites[,1])) {
    cat(paste(">",candidate_sites[k,1], sep=""))
    cat("\n")
    cat(candidate_sites[k,2])
    cat("\n")
    cat(paste(">",candidate_sites[k,3], sep=""))
    cat("\n")
    cat(candidate_sites[k,4])
    cat("\n")
  }
  sink()
  # Run code on terminal and prevent R script from running until terminal script is complete
  terminal_command <- paste("RNAduplex < ", duplex_filename, " > ", output_filename, sep = "")
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
  while(!file.exists(output_filename1)) {
    Sys.sleep(0.1)
  }
  # Get binding energies
  rnadup <- read.csv(output_filename1, header = FALSE, sep = "\t")
  rnadup <- data.frame(lapply(rnadup, as.character), stringsAsFactors=FALSE)
  rnadup <- as.data.frame(rnadup[!grepl(">", rnadup$V1),])
  colnames(rnadup) <- c("mfe")
  rnadup <- data.frame(lapply(rnadup, as.character), stringsAsFactors=FALSE)
  if(length(rnadup[,1]) == 0) {
    next
  }
  for(p in 1:length(rnadup[,1])) {
    x <- as.data.frame(str_locate_all(rnadup$mfe[p], ("\\(")))
    y <- as.data.frame(str_locate_all(rnadup$mfe[p], ("\\)")))
    rnadup$mfe[p] <- as.numeric(substr(rnadup$mfe[p], x[length(x[,1]),1] + 1, y[length(y[,1]),1] - 1))
  }
  # Choose candidate site with lowest binding energy and add to rf_features as a negative site
  rnadup$mfe <- as.numeric(rnadup$mfe)
  rf_features[nrow(rf_features)+1,] <- NA
  rf_features$istarget[nrow(rf_features)] <- 0
  rf_features$trf_id[nrow(rf_features)] <- candidate_sites$trf_id[1]
  rf_features$fixed_trf_sequence[nrow(rf_features)] <- candidate_sites$fixed_trf_sequence[1]
  rf_features$mrna_id[nrow(rf_features)] <- candidate_sites$mrna_id[1]
  rf_features$fixed_mrna_transcript[nrow(rf_features)] <- candidate_sites$mrna_neg_sites[which(grepl(min(rnadup$mfe), rnadup$mfe))[1]]
  rf_features$binding_energy[nrow(rf_features)] <- min(rnadup$mfe)
}