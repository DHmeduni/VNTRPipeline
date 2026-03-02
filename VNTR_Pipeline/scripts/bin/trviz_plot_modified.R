# ---- R Script Header: Install CRAN and Bioconductor Packages ----

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("plotrix", "ggplot2", "openxlsx")
bioc_packages <- c("ShortRead", "Biostrings")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(plotrix)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(openxlsx)

GENETIC_CODE <- Biostrings::GENETIC_CODE

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
motifs_char_sheme <- args[2]
motifs_char_scheme_protein <- args[3]
non_coding <- if (length(args) >= 4 && !is.na(args[4])) args[4] else "N"

count_sequence_lengths <- function(input_file_path) {
  sequences <- readBStringSet(input_file_path)
  seq_names <- gsub("^>", "", names(sequences))
  seq_lengths <- width(sequences)
  data.frame(Sequence = seq_names, `Repeat Number` = seq_lengths, stringsAsFactors = FALSE)
}

translate_dna_to_protein <- function(dna_vec) {
  dna <- DNAStringSet(toupper(dna_vec))
  aa  <- suppressWarnings(translate(dna, if.fuzzy.codon = "X"))
  as.character(aa)
}

count.stop.codons <- function(sequence, reading.frame = 1, genetic.code = GENETIC_CODE) {
  if (!reading.frame %in% c(1,2,3)) stop("reading.frame must be 1, 2, or 3")
  if (class(sequence) != "DNAString") stop("sequence must be a DNAString object")
  if (!("*" %in% genetic.code)) stop("Your genetic code does not specify any stop codons")
  
  l <- length(sequence) + 1 - reading.frame
  if (l < 3) {
    warning(sprintf("Cannot calculate stop codons on sequence of length %d in reading frame %d", 
                    length(sequence), reading.frame))
    return(NULL)
  }
  
  tri <- trinucleotideFrequency(sequence[reading.frame:length(sequence)], step = 3)
  names(tri) <- genetic.code[names(tri)]
  freqs <- sapply(split(tri, names(tri)), sum)
  stops <- freqs["*"]
  return(as.numeric(stops))
}

create_plot <- function(align_output, col_code_out, motifs, mode = c("dna","protein")) {
  mode <- match.arg(mode)
  
  align_output <- rev(align_output)
  colormatrix <- c(as.character(motifs$V3), "#FFA500", "red", "white")
  names(colormatrix) <- c(as.character(motifs$V2), "N", "L", "-")
  
  motifs$V1 <- toupper(as.character(motifs$V1))
  align_output <- sort(align_output, decreasing = TRUE)
  align_motifs <- as.character(align_output)
  names(align_motifs) <- NULL
  align_names <- names(align_output)
  align_names <- gsub(">", "", align_names)
  align_names <- gsub("./", "", align_names)
  
  if (mode == "dna") {
    not_known <- which(is.na(col_code_out$V4))
    if (length(not_known) > 0) {
      for (elem in not_known) { 
        if (non_coding != "Y") {
          if (nchar(col_code_out$V1[elem]) %% 3 != 0 ||
              count.stop.codons(DNAString(col_code_out$V1[elem]), reading.frame = 1) > 0)
            col_code_out$V4[elem] <- "L"
          else
            col_code_out$V4[elem] <- "N"
        }
      }
    }
  }
  
  final_color_code <- character()
  n_seq <- list()
  l_seq <- list()
  i <- 1
  
  for (element in align_motifs) {
    char_vector <- unlist(strsplit(element, split = ""))
    
    pos_match <- match(char_vector, as.character(col_code_out$V2))
    col_code_seq <- col_code_out$V4[as.numeric(pos_match)]
    col_code_seq[is.na(col_code_seq)] <- "-"
    
    aa_label <- col_code_seq
    
    x_vec <- seq(1, length(col_code_seq))
    y_vec <- rep(i, length(col_code_seq))
    final_color_entry <- cbind(col_code_seq, x_vec, y_vec, aa_label)
    final_color_code <- rbind(final_color_code, final_color_entry)
    
    if (mode == "dna") {
      pos_n <- which(col_code_out$V4 == "N")
      pos_l <- which(col_code_out$V4 == "L")
      m_n <- which(!is.na(match(pos_match, pos_n)))
      m_l <- which(!is.na(match(pos_match, pos_l)))
      mc_n <- pos_match[m_n]
      mc_l <- pos_match[m_l]
      
      if (length(m_n) > 0) {
        n_seq[[i]] <- cbind(align_names[i], "N", as.numeric(na.omit(m_n)), 
                            as.character(col_code_out$V1)[as.numeric(mc_n)])
      }
      if (length(m_l) > 0) {
        l_seq[[i]] <- cbind(align_names[i], "L", as.numeric(na.omit(m_l)), 
                            as.character(col_code_out$V1)[as.numeric(mc_l)])
      }
    }
    i <- i + 1
  }
  
  final_color_code <- as.data.frame(final_color_code)
  colnames(final_color_code) <- c("col_code_seq", "x_vec", "y_vec", "aa_label")
  final_color_code$x_vec <- as.numeric(as.character(final_color_code$x_vec))
  final_color_code$y_vec <- as.numeric(as.character(final_color_code$y_vec))
  
  if ((length(n_seq) > 0) & (length(l_seq) > 0)) {
    del <- which(lapply(n_seq, length) < 4)
    if (length(del) > 0) n_seq <- n_seq[-del]
    n_seq <- do.call(rbind.data.frame, n_seq)
    
    del <- which(lapply(l_seq, length) < 4)
    if (length(del) > 0) l_seq <- l_seq[-del]
    l_seq <- do.call(rbind.data.frame, l_seq)
    t_seq <- rbind(n_seq, l_seq)
  } else if (length(n_seq) > 0 & (length(l_seq) == 0)) {
    del <- which(lapply(n_seq, length) < 4)
    if (length(del) > 0) n_seq <- n_seq[-del]
    n_seq <- do.call(rbind.data.frame, n_seq)
    t_seq <- n_seq
  } else if (length(l_seq) > 0 & (length(n_seq) == 0)) {
    del <- which(lapply(l_seq, length) < 4)
    if (length(del) > 0) l_seq <- l_seq[-del]
    l_seq <- do.call(rbind.data.frame, l_seq)
    t_seq <- l_seq
  } else {
    t_seq <- data.frame(Haplotype = NA, `N/L` = NA, `Motif position` = NA, Sequence = NA,
                        stringsAsFactors = FALSE)
  }
  
  n_rows <- length(align_names)
  n_cols <- max(final_color_code$x_vec)
  base_size <- 0.02
  max_size <- 3
  text_size <- min(max_size, max_size - (base_size * n_rows))
  
  p <- ggplot(final_color_code, aes(x_vec, y_vec, fill = col_code_seq)) + 
    geom_tile(color = "white", show.legend = FALSE) +
    coord_fixed(ratio = 1) +
    geom_text(aes(label = aa_label), size = text_size) +
    scale_fill_manual(values = colormatrix) +
    theme(
      axis.text.y = element_text(size = text_size * 3),
      axis.text.x = element_text(size = text_size * 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(1, length(final_color_code[, 1]) / length(align_names)),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(1, length(align_names)), labels = align_names,
                       expand = c(0.0, 0.0))
  
  colnames(t_seq) <- c("Haplotype", "N/L", "Motif position", "Sequence")
  return(list(data = t_seq, plot = p))
}

# ---- Protein-specific preprocessing: classify motifs BEFORE mapping ----
preprocess_protein_motifs <- function(col_code_out, motifs_protein) {
  # Step 1: LOF check (length not div by 3 OR has stop codon) -> mark as L
  dna_motifs <- toupper(as.character(col_code_out$V1))
  non_div3 <- nchar(dna_motifs) %% 3 != 0
  aa_motifs <- translate_dna_to_protein(dna_motifs)
  has_stop <- grepl("\\*", aa_motifs)
  
  # Create classification column
  col_code_out$V4_PROT_CLASS <- "-"
  lof_idx <- which(non_div3 | has_stop)
  col_code_out$V4_PROT_CLASS[lof_idx] <- "L"
  
  # Step 2: For non-LOF motifs, check if protein seq in motifs_protein$V1 -> N if not found
  protein_motifs <- toupper(as.character(motifs_protein$V1))
  non_lof_idx <- which(col_code_out$V4_PROT_CLASS == "-")
  
  for (i in non_lof_idx) {
    prot_seq <- toupper(col_code_out$V1_PROT[i])
    if (!is.na(prot_seq) && prot_seq != "" && !(prot_seq %in% protein_motifs)) {
      col_code_out$V4_PROT_CLASS[i] <- "N"
    }
  }
  
  return(col_code_out)
}

# ---- Main execution ----

input_alignment_file <- paste0(work_dir, "/output_TRViz/best_trviz_alignment_input.fa")
vntr_lengths <- count_sequence_lengths(input_alignment_file)
write.xlsx(vntr_lengths, paste0(work_dir, "/output_TRViz/vntr_length.xlsx"))

align_output <- readBStringSet(paste0(work_dir, "/output_TRViz/best_trviz_alignment_output.fa"))
col_code_out <- read.table(paste0(work_dir, "/output_TRViz/best_trviz_motif_map.txt"), stringsAsFactors = FALSE)
motifs <- read.table(motifs_char_sheme, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
if (file.exists(motifs_char_scheme_protein)) {
  motifs_protein <- read.table(motifs_char_scheme_protein, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
}

motifs$V1 <- toupper(as.character(motifs$V1))
if (exists("motifs_protein")) {
  motifs_protein$V1 <- toupper(as.character(motifs_protein$V1))
}

# DNA mapping: V1 (DNA) -> V4_DNA via motifs
col_code_out$V1 <- toupper(as.character(col_code_out$V1))
col_code_out$V4_DNA <- as.character(motifs$V2[match(col_code_out$V1, motifs$V1)])

# DNA plot: use V4_DNA  
col_code_out$V4 <- col_code_out$V4_DNA
res_dna <- create_plot(align_output, col_code_out, motifs, mode = "dna")
p_dna <- res_dna$plot

# New/LOF sequences from DNA plot data
newlof_output <- readBStringSet(paste0(work_dir, "/output_TRViz/best_trviz_alignment_input.fa"))
col_code_out$V4 <- col_code_out$V4_DNA
newlof <- create_plot(newlof_output, col_code_out, motifs, mode = "dna")
seq_dat <- newlof$data
write.xlsx(seq_dat, paste0(work_dir, "/output_TRViz/new_and_lof_seqs.xlsx"))

# ---- Write DNA motif-symbol FASTA (parallel to protein) ----

dna_haplotype_codes <- list()

align_names <- gsub(">", "", names(align_output))
align_names <- gsub("./", "", align_names)
align_output_sorted <- sort(align_output, decreasing = TRUE)

for (i in seq_along(align_names)) {
  haplotype_seq <- as.character(align_output_sorted[i])
  char_vector <- unlist(strsplit(haplotype_seq, split = ""))

  pos_match <- match(char_vector, as.character(col_code_out$V2))
  dna_codes <- col_code_out$V4_DNA[as.numeric(pos_match)]
  dna_codes[is.na(dna_codes)] <- "-"

  clean_seq <- paste(dna_codes[dna_codes != "-"], collapse = "")
  dna_haplotype_codes[[align_names[i]]] <- clean_seq
}

dna_fasta_file <- paste0(work_dir, "/output_TRViz/best_trviz_dna_motif_sequence.fa")
sink(dna_fasta_file)
for (haplotype in names(dna_haplotype_codes)) {
  cat(">", haplotype, "\n")
  cat(dna_haplotype_codes[[haplotype]], "\n")
}
sink()

# Create Plot

ggsave(filename = paste0(work_dir, "/output_TRViz/best_trviz_fig.png"),
       plot = p_dna, width = 15, height = 6, dpi = 300)

# Create Protein Sequence unless non_coding

if (non_coding != "Y" && motifs_char_scheme_protein != "N") {

  # Protein translation
  col_code_out$V1_PROT <- translate_dna_to_protein(col_code_out$V1)

  # PROTEIN-SPECIFIC: Classify motifs first (L/N logic)
  col_code_out <- preprocess_protein_motifs(col_code_out, motifs_protein)

  # Protein motif mapping: V1_PROT -> V4_PROT
  col_code_out$V4_PROT <- as.character(
    motifs_protein$V2[match(col_code_out$V1_PROT, motifs_protein$V1)]
  )

  # Override motif codes with L/N classification where applicable
  col_code_out$V4_PROT[col_code_out$V4_PROT_CLASS == "L"] <- "L"
  col_code_out$V4_PROT[col_code_out$V4_PROT_CLASS == "N"] <- "N"
  col_code_out$V4_PROT[is.na(col_code_out$V4_PROT)] <- "-"

  # --- Diagnostics ---
  cat("=== PROTEIN MAPPING DIAGNOSTICS ===\n")
  print(head(col_code_out[, c("V1", "V1_PROT", "V4_PROT_CLASS", "V4_PROT")]))
  print(table(col_code_out$V4_PROT_CLASS, useNA = "always"))
  print(table(col_code_out$V4_PROT, useNA = "always"))

  # Protein plot
  col_code_out$V4 <- col_code_out$V4_PROT
  res_prot <- create_plot(align_output, col_code_out, motifs_protein, mode = "protein")
  p_prot <- res_prot$plot

  # Reconstruct protein haplotypes
  protein_haplotype_seqs <- list()
  align_names <- gsub(">", "", names(align_output))
  align_names <- gsub("./", "", align_names)
  align_output_sorted <- sort(align_output, decreasing = TRUE)

  for (i in seq_along(align_names)) {
    haplotype_seq <- as.character(align_output_sorted[i])
    char_vector <- unlist(strsplit(haplotype_seq, split = ""))
    pos_match <- match(char_vector, as.character(col_code_out$V2))
    prot_codes <- col_code_out$V4[as.numeric(pos_match)]
    prot_codes[is.na(prot_codes)] <- "-"
    clean_seq <- paste(prot_codes[prot_codes != "-"], collapse = "")
    protein_haplotype_seqs[[align_names[i]]] <- clean_seq
  }

  # Write FASTA
  fasta_file <- paste0(work_dir, "/output_TRViz/best_trviz_protein_motif_sequence.fa")
  sink(fasta_file)
  for (haplotype in names(protein_haplotype_seqs)) {
    cat(">",haplotype,"\n")
    cat(protein_haplotype_seqs[[haplotype]], "\n")
  }
  sink()

  cat("Wrote", length(protein_haplotype_seqs),
      "protein haplotype sequences to:\n", fasta_file, "\n")

  # Save protein figure
  ggsave(filename = paste0(work_dir, "/output_TRViz/best_trviz_protein_fig.png"),
         plot = p_prot, width = 15, height = 6, dpi = 300)
}







