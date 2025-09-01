# Set CRAN mirror to avoid errors in non-interactive sessions
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of required packages and their sources
cran_packages <- c("stringi")
bioc_packages <- c("Biostrings", "ShortRead")

# Install CRAN packages if missing
for(pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Install Bioconductor packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for(pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

library(Biostrings)
library(ShortRead)
library(stringi)
args <- commandArgs(trailingOnly = TRUE)
work_dir=args[1]
setwd(work_dir)
all_sample_names=list.dirs(recursive=FALSE)

if (!dir.exists("output_TRViz")) dir.create("output_TRViz")

for (sample_name in all_sample_names)
{
  
  
  if (file.exists(paste0(sample_name,"/",sample_name,".trimmedReads_rev.fasta"))&&file.exists(paste0(sample_name,"/",sample_name,".trimmedReads_for.fasta")))
  {
    dna_set=readDNAStringSet(paste0(sample_name,"/",sample_name,".trimmedReads_for.fasta"))
    dna_set_rev=readDNAStringSet(paste0(sample_name,"/",sample_name,".trimmedReads_rev.fasta"))
   
  
    dna_set_rev=reverseComplement(dna_set_rev)
    all_dna_set=c(dna_set,dna_set_rev)
  }
  else if (file.exists(paste0(sample_name,"/",sample_name,".trimmedReads_for.fasta")))
      all_dna_set=readDNAStringSet(paste0(sample_name,"/",sample_name,".trimmedReads_for.fasta"))
  else if (file.exists(paste0(sample_name,"/",sample_name,".trimmedReads_rev.fasta")))
  {
    dna_set_rev=readDNAStringSet(paste0(sample_name,"/",sample_name,".trimmedReads_rev.fasta"))
    dna_set_rev=reverseComplement(dna_set_rev)
    all_dna_set=dna_set_rev
  }
  else
  {
    next
  }
  all_dna_table=sort(table(all_dna_set),decreasing = TRUE)
  dna_table_char=as.character(names(all_dna_table))
  sample_name_clean=sub("^\\./", "", sample_name)
  seq_names=paste(sample_name_clean, "_seq",seq(1:length(all_dna_table)), " (",as.character(all_dna_table),")",sep="")
  names(dna_table_char)=seq_names
  
  id_first=unlist(lapply(strsplit(sample_name,"_"),"[[",1))
  id_last=stri_sub(sample_name,-1,-1)
  
  
  writeFasta(DNAStringSet(dna_table_char),paste0("output_TRViz/",sample_name,"_result.fasta"))
  
  best_hit=dna_table_char[1]
  names(best_hit)=sample_name_clean
  writeFasta(DNAStringSet(best_hit),paste0("output_TRViz/",sample_name,"_best_hit.fasta"))
}

