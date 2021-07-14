library(stringr)
library(stringi)

# "Species" is one of the species found in tRFdb
species_trf_info <- read.csv("species_trf_info.csv", row.names = 1)
species_trf_info = data.frame(lapply(species_trf_info, as.character), stringsAsFactors=FALSE)

# Collapse and process info to remove repeated tRFs
species_trf_info = species_trf_info[match(unique(species_trf_info$tRF.ID), species_trf_info$tRF.ID),]
species_trf_info$tRF.ID = trimws(species_trf_info$tRF.ID)
species_trf_info$Type = trimws(species_trf_info$Type)
species_trf_info$tRF.Sequence = trimws(species_trf_info$tRF.Sequence)

# Calculate binding cutoff for tRFs - 60% of energy of perfect tRF-mRNA duplex ligated by "AAAA" sequence
sink("species_trf_rnafold.txt")
for(k in 1:length(species_trf_info[,1])) {
  trf_seq = species_trf_info$tRF.Sequence[k]
  comp_seq = ""
  for(i in 1:nchar(trf_seq)) {
    if(substr(trf_seq, i, i) == "C") {
      comp_seq <- paste(comp_seq, "G", sep = "")
    }
    else if(substr(trf_seq, i, i) == "G"){
      comp_seq <- paste(comp_seq, "C", sep = "")
    }
    else if(substr(trf_seq, i, i) == "T"){
      comp_seq <- paste(comp_seq, "A", sep = "")
    }
    else if(substr(trf_seq, i, i) == "A"){
      comp_seq <- paste(comp_seq, "T", sep = "")
    }
  }
  reverse_comp_seq = stri_reverse(comp_seq)
  cat(paste(trf_seq,"AAAA",reverse_comp_seq, sep = ""))
  cat("\n")
}
sink()

# Uses RNAfold (ViennaRNA)
terminal_command <- paste("RNAfold --noPS < ", "species_trf_rnafold.txt", " > ", "species_trf_rnafold_out.txt", sep = "")
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

species_rnafold = read.csv("species_trf_rnafold_out.txt", header = FALSE)
species_rnafold = as.character(species_rnafold[seq(2,length(species_rnafold[,1]), 2), ])
species_trf_info$binding_cutoff = ceiling(0.6*as.numeric(substr(species_rnafold, str_locate(species_rnafold, "-"), nchar(species_rnafold)-1)))
species_trf_info$tRF.Sequence = trimws(species_trf_info$tRF.Sequence)
rownames(species_trf_info) = seq(1,length(species_trf_info[,1]))
write.csv(species_trf_info, "species_trf_info.csv")

# Process all_genes file obtained from Ensembl
biomart_import = read.csv("all_species_genes.txt", header = FALSE)
biomart_import = data.frame(lapply(biomart_import, as.character), stringsAsFactors=FALSE)

all_species_genes = data.frame(matrix(nrow = length(biomart_import[,1])/2, ncol = 2))
all_species_genes$X1 = biomart_import[c(seq(1,length(biomart_import[,1]),2)),]
all_species_genes$X2 = biomart_import[c(seq(2,length(biomart_import[,1]),2)),]

all_species_genes$X1 = substr(all_species_genes$X1, 2, nchar(all_species_genes$X1))
x = str_locate_all(all_species_genes$X1, "\\|")

all_species_genes$tran_id = ""
all_species_genes$tran_ver = ""
all_species_genes$name = ""
all_species_genes$chr = ""
all_species_genes$gene_start_bp = ""
all_species_genes$gene_end_bp = ""
all_species_genes$transcript_start_bp = ""
all_species_genes$transcript_end_bp = ""
all_species_genes$utr3_start_bp = ""
all_species_genes$utr3_end_bp = ""
all_species_genes$length_seq = ""

for(k in 1:length(all_species_genes[,1])) {
  if(k%%1000==0){print(k)}
  seq = all_species_genes$X1[k]
  all_species_genes$tran_id[k] = substr(seq, 1, x[[k]][1] - 1)
  all_species_genes$tran_ver[k] = substr(seq, x[[k]][1] + 1, x[[k]][2] - 1)
  all_species_genes$name[k] = substr(seq, x[[k]][2] + 1, x[[k]][3] - 1)
  all_species_genes$chr[k] = substr(seq, x[[k]][3] + 1, x[[k]][4] - 1)
  all_species_genes$gene_start_bp[k] = substr(seq, x[[k]][4] + 1, x[[k]][5] - 1)
  all_species_genes$gene_end_bp[k] = substr(seq, x[[k]][5] + 1, x[[k]][6] - 1)
  all_species_genes$transcript_start_bp[k] = substr(seq, x[[k]][6] + 1, x[[k]][7] - 1)
  all_species_genes$transcript_end_bp[k] = substr(seq, x[[k]][7] + 1, x[[k]][8] - 1)
  all_species_genes$utr3_start_bp[k] = substr(seq, x[[k]][8] + 1, x[[k]][9] - 1)
  all_species_genes$utr3_end_bp[k] = substr(seq, x[[k]][9] + 1, nchar(seq))
}
colnames(all_species_genes) = c("all", "sequence", colnames(all_species_genes)[3:13])
all_species_genes = all_species_genes[,c(3,4,5,2,6,7,8,9,10,11,12,13)]
all_species_genes$length_seq = nchar(all_species_genes$sequence)

# Remove transcripts with unavailable sequence data
all_species_genes = all_species_genes[which(all_species_genes$sequence != "Sequence unavailable"),]

# Remove genes shorter than longest tRF for the species
max_trf_length = max(nchar(as.character(species_trf_info$tRF.Sequence)))
all_species_genes = all_species_genes[which(all_species_genes$length_seq >= max_trf_length),]

write.csv(all_species_genes, "all_species_genes_curated.csv")


