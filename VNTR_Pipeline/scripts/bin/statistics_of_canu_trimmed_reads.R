
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
work_dir=args[1]

setwd(paste0(work_dir,"/output_TRViz/"))
all_seqs=list.files("./","result.fasta")

seqs_dist=as.character()
frac_all=as.character()
frac_second_all=as.character()
i=1
for (seqs in all_seqs)
{
  sread(readFasta(seqs))
  names=ShortRead::id(readFasta(seqs))
  names=gsub("\\(","z",names)
  names=gsub("\\)","z",names)
  dist=as.numeric(unlist(lapply(strsplit(as.vector(names),"z"),"[[",2)))
  frac=sum(dist[-1])/dist[1]
  frac_second=dist[2]/dist[1]
  seqs_dist[i]=paste(unlist(lapply(strsplit(as.vector(names),"z"),"[[",2)), collapse=", ")
  frac_all[i]=frac
  frac_second_all[i]=frac_second
  i=i+1
}

seqs_dat=cbind(seqs_dist,frac_all,frac_second_all)
rownames(seqs_dat)=all_seqs
seqs_dat=as.data.frame(seqs_dat)
seqs_dat[,2]=round(as.numeric(as.character(seqs_dat[,2])),digits=2)
seqs_dat[,3]=round(as.numeric(as.character(seqs_dat[,3])),digits=2)

colnames(seqs_dat)=c("sequences distribution","sum of other seqs/# of most freq seq","# of second most freq seq/# of most freq seq")
write.xlsx(seqs_dat, "seqs_distribution.xlsx", rowNames=TRUE)

