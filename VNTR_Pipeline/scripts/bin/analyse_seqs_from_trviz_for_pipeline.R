
# Set CRAN mirror to avoid errors in non-interactive sessions
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Define packages and their sources
cran_packages <- c("openxlsx")
bioc_packages <- c("ShortRead")

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
library(openxlsx)
library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)
motifs_path=args[1]
work_dir=args[2]


#FUNCTION, count stop codons ----------------------------------------------
count.stop.codons <- function(sequence, reading.frame = 1, genetic.code = GENETIC_CODE){
  
  if(!reading.frame %in% c(1,2,3)){ stop("reading.frame must be 1, 2, or 3")}
  if(class(sequence)!='DNAString'){ stop("sequence must be a DNAString object")}
  if(!("*" %in% genetic.code)) { stop("Your genetic code does not specify any stop codons")}
  
  l = length(sequence) + 1 - reading.frame
  
  if(l < 3){
    warning(sprintf("Cannot calculate stop codons on sequence of length %d in reading frame %d", 
                    length(sequence), reading.frame))
    return(NULL)
  }
  
  # this comes almost straight from the BioStrings manual
  tri = trinucleotideFrequency(sequence[reading.frame:length(sequence)], step=3)
  
  names(tri) <- genetic.code[names(tri)]
  
  freqs = sapply(split(tri, names(tri)), sum)
  
  stops = freqs["*"]
  
  return(as.numeric(stops))
}
#--------------------------------------

#motiv evaluation
motiv_map=paste0(work_dir,"/output_TRViz/best_trviz_motif_map.txt")
input_seqs=read.table(motiv_map,sep="\t")

#all_motivs=paste0(motifs_path,"/motifs_muc1.txt")
all_motivs_in_seqs=read.table(paste0(work_dir,"/output_TRViz/best_trviz_alignment_input.fa"))
all_motivs=read.table(motifs_path)
all_motivs=toupper(all_motivs$V1)

matched_motifs=match(input_seqs$V1,all_motivs)
dat_for_xls=cbind(input_seqs,matched_motifs)

colnames(dat_for_xls)=c("Sequence","Seq ID","# Ocurrences","ID of sequence of known motifs")
write.xlsx(dat_for_xls,paste0(work_dir,"/motifs_evaluation.xlsx"))


#check for frameshifts and stops.
all_sample_names=as.data.frame(all_motivs_in_seqs$V1[seq(1, length(all_motivs_in_seqs$V1), 2)])
all_seqs=all_motivs_in_seqs$V1[seq(2, length(all_motivs_in_seqs$V1), 2)]
all_seqs=as.character(all_seqs)

all_compl_seqs=character()
stop_counts=character()
for (seqs in all_seqs)
{
  complete_seq=character()
  for (i in 1:nchar(seqs))
  {
    current_seq=as.character(input_seqs[match(substring(seqs,i,i),input_seqs$V2),1])
    complete_seq=paste0(complete_seq,current_seq)
  }
  all_compl_seqs=c(all_compl_seqs,complete_seq)
  stop_counts=c(stop_counts,count.stop.codons(DNAString(complete_seq),reading.frame=1)) #change here if seq does not start with first base of amino acid
}


in_frame_seq=nchar(all_compl_seqs)%%3==0 # change here if seq is not dividable by 3


#lof_check=cbind(all_sample_names,in_frame_seq,stop_counts)
#colnames(lof_check)=c("File Name","In Frame?","Stop Codon Count")

#write.xlsx(lof_check, paste0(work_dir,"/LOF_check.xlsx"))

