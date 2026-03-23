create_fimo <- function(Eshijixuliequchong,TSSshijixuliequchong,PWM_E,PWM_TSS,E_motif_index,TSS_motif_index,output_folder_path){
  
  Eshijixulie <- Eshijixuliequchong
  TSSshijixulie <- TSSshijixuliequchong
  
  Eshijixulie <- unlist(Eshijixulie)
  TSSshijixulie <- unlist(TSSshijixulie)
  
  E_pro_A <- length(which(Eshijixulie==0))/length(Eshijixulie)
  E_pro_T <- length(which(Eshijixulie==1))/length(Eshijixulie)
  E_pro_C <- length(which(Eshijixulie==2))/length(Eshijixulie)
  E_pro_G <- length(which(Eshijixulie==3))/length(Eshijixulie)
  cat("E_shijigailv:","E_pro_A:",E_pro_A,"E_pro_T:",E_pro_T,"E_pro_C:",E_pro_C,"E_pro_G:",E_pro_G,"\n")
  
  TSS_pro_A <- length(which(TSSshijixulie==0))/length(TSSshijixulie)
  TSS_pro_T <- length(which(TSSshijixulie==1))/length(TSSshijixulie)
  TSS_pro_C <- length(which(TSSshijixulie==2))/length(TSSshijixulie)
  TSS_pro_G <- length(which(TSSshijixulie==3))/length(TSSshijixulie)
  cat("TSS_shijigailv:","TSS_pro_A:",TSS_pro_A,"TSS_pro_T:",TSS_pro_T,"TSS_pro_C:",TSS_pro_C,"TSS_pro_G:",TSS_pro_G,"\n")
  
  quchong_index_E <- E_motif_index
  
  
  for(i in 1:length(E_motif_index)){
    quchong_index_E[[i]] <- E_motif_index[[i]][1]
  }
  
  quchong_index_TSS <- TSS_motif_index
  for(i in 1:length(TSS_motif_index)){
    quchong_index_TSS[[i]] <- TSS_motif_index[[i]][1]
  }
  
  quchong_index_E <- unlist(quchong_index_E)
  quchong_index_TSS <- unlist(quchong_index_TSS)
  
  
  quchong_motif_E <- PWM_E[quchong_index_E]
  quchong_motif_TSS <- PWM_TSS[quchong_index_TSS]
  
  cat("MEME version 5.5.0\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"))
  cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("ALPHABET = ACGT\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("strands : + -\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("Background letter frequencies\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("A" ,E_pro_A ,"C" ,E_pro_C ,"G", E_pro_G, "T" ,E_pro_T,"\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
  
  index <- "MOTIF MA000_"
  for(i in 1:length(quchong_motif_E)){
    cat("i:",i,"\n")
    a <- paste(index,i,sep="")
    cat(a, file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    w <- nrow(quchong_motif_E[[i]])
    cat("letter-probability matrix: alength= 4", "w=", w ,"nsites= 20 E= 0", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    for(j in 1:w){
      cat(quchong_motif_E[[i]][j,], file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
      cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
      
    }
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_E_FIMO.txt",sep="/"),append=TRUE )
    
  }
  
  ##TSS
  cat( "MEME version 5.5.0\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"))
  cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("ALPHABET = ACGT\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("strands : + -\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
  cat("Background letter frequencies\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("A" ,TSS_pro_A ,"C" ,TSS_pro_C ,"G", TSS_pro_G, "T" ,TSS_pro_T,"\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE)
  cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
  
  for(i in 1:length(quchong_motif_TSS)){
    cat("i:",i,"\n")
    a <- paste(index,i,sep="")
    cat(a, file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    w <- nrow(quchong_motif_TSS[[i]])
    cat("letter-probability matrix: alength= 4", "w=", w ,"nsites= 20 E= 0", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    for(j in 1:w){
      cat(quchong_motif_TSS[[i]][j,], file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
      cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
      
    }
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    cat("\n", file=paste(output_folder_path,"Motif_TSS_FIMO.txt",sep="/"),append=TRUE )
    
  }
  
  
  return(invisible(list(
    motif_E = file.path(output_folder_path, "Motif_E_FIMO.txt"),
    motif_TSS = file.path(output_folder_path, "Motif_TSS_FIMO.txt")
  )))
  
  
}