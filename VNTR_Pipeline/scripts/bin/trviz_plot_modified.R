# ---- R Script Header: Install CRAN and Bioconductor Packages ----

# Install BiocManager if needed (for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("plotrix", "ggplot2", "openxlsx")
bioc_packages <- c("ShortRead")

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load all libraries
library(plotrix)
library(ShortRead)
library(ggplot2)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
work_dir=args[1]
motifs_char_sheme=args[2]
non_coding=args[3]

#function
create_plot=function(align_output,col_code_out)
{
  align_output=rev(align_output)
  
  colormatrix=c(as.character(motifs$V3), "#FFA500","red","white")
  
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
  
  
  
  motifs$V1=toupper(as.character(motifs$V1))
  align_output = sort(align_output, decreasing = TRUE)
  align_motifs=as.character(align_output)
  #align_motifs=sort(align_motifs)
  names(align_motifs)=NULL
  align_names=names(align_output)
  align_names=gsub(">","",align_names)
  align_names=gsub("./","",align_names)
  
  #--------------------
  col_code_out$V4=as.character(motifs$V2[match(col_code_out$V1, motifs$V1)])
  not_known=which(is.na(col_code_out$V4))

  #print(non_coding)
  #readline(prompt = "Press [Enter] to continue...")

  if (length(not_known)> 0)
  {
    for (elem in not_known)
    { 
        if (non_coding!="Y")
	    {
      if(nchar(as.character(col_code_out$V1)[elem])%%3!=0 | count.stop.codons(DNAString(col_code_out$V1[elem]),reading.frame=1)>0)
        col_code_out$V4[elem]="L"
      else
        col_code_out$V4[elem]="N"
      }
       else
	    col_code_out$V4[elem]="N"
    }
  }
  
  final_color_entry=character()
  final_color_code=character()
  
  i=1
  n_seq=list()
  l_seq=list()
  element=align_motifs[2]
  for (element in align_motifs)
  {
    char_vector=unlist(strsplit(element, split=""))
    pos_match=character()
    for (char in char_vector)
    {
      pos_match=c(pos_match,match(char, as.character(col_code_out$V2)))
    }
    col_code_seq=col_code_out$V4[as.numeric(pos_match)]
    col_code_seq[is.na(col_code_seq)]="-"
    x_vec=seq(1,length(col_code_seq))
    y_vec=rep(i,length(col_code_seq))
    final_color_entry=cbind(col_code_seq,x_vec,y_vec)
    final_color_code=rbind(final_color_code,final_color_entry)
    
     #for new/lof motif seqs
     pos_n=which(col_code_out$V4=="N")
     pos_l=which(col_code_out$V4=="L")
     m_n=which(!is.na(match(pos_match,pos_n)))
     m_l=which(!is.na(match(pos_match,pos_l)))
     mc_n=pos_match[m_n]
     mc_l=pos_match[m_l]

     if (length(m_n)>0)
     {
    n_seq[[i]]=cbind(align_names[i],"N",as.numeric(na.omit(m_n)),as.character(col_code_out$V1)[as.numeric(mc_n)])
     }
     if (length(m_l)>0)
     {
    #l_seq=rbind(l_seq,cbind(align_names[i],"L,",as.numeric(na.omit(m_l)),as.character(col_code_out$V1)[pos_l[which(!is.na(m_l))]]))
    l_seq[[i]]=cbind(align_names[i],"L",as.numeric(na.omit(m_l)),as.character(col_code_out$V1)[as.numeric(mc_l)])
     }
     i=i+1
   }

  final_color_code=as.data.frame(final_color_code)
  final_color_code$x_vec=as.numeric(as.character(final_color_code$x_vec))
  final_color_code$y_vec=as.numeric(as.character(final_color_code$y_vec))
  
  
  names(colormatrix)=c(as.character(motifs$V2),"N","L","-")
  
  #4 possibilities
#1 Ns and Ls present
if ((length(n_seq)>0) & (length(l_seq)>0))
{
  del=which(lapply(n_seq, length)<4)
  if (length(del) > 0)
    n_seq <- n_seq[-del]
  
  n_seq <- do.call(rbind.data.frame, n_seq)
  
  del=which(lapply(l_seq, length)<4)
  if (length(del) > 0)
    l_seq <- l_seq[-del]
  
  l_seq <- do.call(rbind.data.frame, l_seq)
  t_seq <- rbind(n_seq, l_seq)
}
#Only Ns present
else if (length(n_seq)>0 & (length(l_seq)==0))
{
  del=which(lapply(n_seq, length)<4)
  if (length(del)>0)
    n_seq <- n_seq[-del]
  
  n_seq <- do.call(rbind.data.frame, n_seq)
  t_seq <- n_seq
}
#Only Ls are present
else if (length(l_seq)>0 & (length(n_seq)==0))
{
  del=which(lapply(l_seq, length)<4)
  if (length(del)>0)
    l_seq <- l_seq[-del]
  
  l_seq <- do.call(rbind.data.frame, l_seq)
  t_seq <- l_seq
}
#none is present
else
{
  t_seq <- data.frame(
    Haplotype = NA,
    `N/L` = NA,
    `Motif position` = NA,
    Sequence = NA,
    stringsAsFactors = FALSE
  )
}

  
  n_rows <- length(align_names)
  n_cols <- max(final_color_code$x_vec)
  base_size <- 0.02
  max_size <- 3
  text_size <- min(max_size, max_size - (base_size * n_rows) )

  p <- ggplot(final_color_code, aes(x_vec, y_vec, fill = col_code_seq)) + 
    geom_tile(color = 'white', show.legend = FALSE) +
    coord_fixed(ratio = 1) +    # maintain a fixed aspect ratio so rectangles are wider than tall coord_fixed(ratio = 1 / 1) +
    geom_text(aes(label = col_code_seq), size = text_size) +
    #geom_text(aes(label = col_code_seq), cex = 2.5) + scale_fill_manual(values = colormatrix) +
    scale_fill_manual(values = colormatrix) +
    theme(
      axis.text.y = element_text(size = text_size * 3),
      axis.text.x = element_text(size = text_size * 2 ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    scale_x_continuous(breaks = seq(1:length(final_color_code[, 1]) / length(align_names)),
                       expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(1:length(align_names)),
      labels = align_names,
      expand = c(0.0, 0.0)
    )
    
  colnames(t_seq)=c("Haplotype","N/L","Motif position","Sequence")
  return(list(data = t_seq, plot = p))
}



align_output=readBStringSet(paste0(work_dir,"/output_TRViz/","best_trviz_alignment_output.fa"))
col_code_out=read.table(paste0(work_dir,"/output_TRViz/","best_trviz_motif_map.txt"))
motifs=read.table(motifs_char_sheme)

res <- create_plot(align_output, col_code_out)
seq_dat <- res$data
p <- res$plot
#seq_dat=create_plot(align_output,col_code_out)
write.xlsx(seq_dat,paste0(work_dir,"/output_TRViz/","new_and_lof_seqs.xlsx"))

ggsave(
  filename = paste0(work_dir, "/output_TRViz/best_trviz_fig.png"),
  plot = p,
  width = 15,
  height = 6,
  dpi = 300
)
#ggsave(paste0(work_dir,"/output_TRViz/","best_trviz_fig.png"))


