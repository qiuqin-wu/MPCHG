pair_distance  <- function(PWM_E,PWM_TSS,E_motif_index,TSS_motif_index,fimo_E,peaks_E,fimo_TSS,peaks_TSS){
  
  Eshijisitequan <- peaks_E
  TSSshijisitequan <- peaks_TSS
  
  fimo_E$chr_start <- 0
  fimo_E$chr_end <- 0
  fimo_TSS$chr_start <- 0
  fimo_TSS$chr_end <- 0
  
  
  for(i in 1:length(fimo_E$motif_id)){
    index <- which((fimo_E$start[i]  >= Eshijisitequan$V2) & (fimo_E$stop[i] <= Eshijisitequan$V3))
    fimo_E$chr_start[i] <- Eshijisitequan[index,2]
    fimo_E$chr_end[i] <- Eshijisitequan[index,3]
    
  }
  colnames(fimo_E)[3] <- "chr"  
  
  fimo_E$sequence_name <- paste(fimo_E$chr,fimo_E$chr_start,sep=":")
  fimo_E$sequence_name <- paste(fimo_E$sequence_name,fimo_E$chr_end,sep="-")
  fimo2_E <- fimo_E[,c(1,13,4,5,3,11,12)]
  fimo_E <- fimo2_E
  data_E <- fimo_E
  
  for(i in 1:length(fimo_TSS$motif_id)){
    index <- which((fimo_TSS$start[i]  >= TSSshijisitequan$V2) & (fimo_TSS$stop[i] <= TSSshijisitequan$V3))
    fimo_TSS$chr_start[i] <- TSSshijisitequan[index,2]
    fimo_TSS$chr_end[i] <- TSSshijisitequan[index,3]
    
  }
  colnames(fimo_TSS)[3] <- "chr"  
  
  fimo_TSS$sequence_name <- paste(fimo_TSS$chr,fimo_TSS$chr_start,sep=":")
  fimo_TSS$sequence_name <- paste(fimo_TSS$sequence_name,fimo_TSS$chr_end,sep="-")
  fimo2_TSS <- fimo_TSS[,c(1,13,4,5,3,11,12)]
  fimo_TSS <- fimo2_TSS
  data_TSS <- fimo_TSS
  
  library(Matrix)
  
  paste1 <- paste(Eshijisitequan$V1,Eshijisitequan$V2,sep=":")
  paste2 <- paste(paste1,Eshijisitequan$V3,sep="-")
  Eshijisitequan$V4 <- paste2
  Eshijisitequan$V2 <- as.numeric(Eshijisitequan$V2)
  Eshijisitequan$V3 <- as.numeric(Eshijisitequan$V3)
  colnames(Eshijisitequan) <- c("E_chr","E_start","E_end","E_sequence_name")
  
  paste3 <- paste(TSSshijisitequan$V1,TSSshijisitequan$V2,sep=":")
  paste4 <- paste(paste3,TSSshijisitequan$V3,sep="-")
  TSSshijisitequan$V4 <- paste4
  TSSshijisitequan$V2 <- as.numeric(TSSshijisitequan$V2)
  TSSshijisitequan$V3 <- as.numeric(TSSshijisitequan$V3)
  colnames(TSSshijisitequan) <- c("TSS_chr","TSS_start","TSS_end","TSS_sequence_name")
  
  library(stringr)
  fimo_E <- data_E
  
  fimo_TSS <- data_TSS
  
  motif_name_E <- vector(length=length(PWM_E))
  name <- "MA000_"
  for(i in 1:length(E_motif_index)){
    index <- E_motif_index[[i]][1]
    name2 <- paste(name,index,sep="")
    motif_name_E[E_motif_index[[i]]] <- name2
  }
  
  
  motif_name_TSS <- vector(length=length(PWM_TSS))
  for(i in 1:length(TSS_motif_index)){
    index <- TSS_motif_index[[i]][1]
    name2 <- paste(name,index,sep="")
    motif_name_TSS[TSS_motif_index[[i]]] <- name2
  }
  
  
  distance_pair <- function(x){
    
    name_E <- motif_name_E[x]
    name_TSS <- motif_name_TSS[x]                  
    #index_E_8000 <- which(fimo_E_8000$motif_id==name_E)
    index_E <- which(fimo_E$motif_id==name_E)
    
    index_TSS <- which(fimo_TSS$motif_id==name_TSS)
    
    fimo_E_new <- fimo_E[index_E,]
    
    
    fimo_TSS_new <- fimo_TSS[index_TSS,]  
    
    
    sequence_name_E <- unique(fimo_E_new$sequence_name)
    
    
    sequence_name_TSS <- unique(fimo_TSS_new$sequence_name)
    
    
    
    
    index_E_inshiji <- which(Eshijisitequan$E_sequence_name %in% sequence_name_E)
    index_TSS_inshiji <- which(TSSshijisitequan$TSS_sequence_name %in% sequence_name_TSS)
    index_common <- which(index_E_inshiji %in% index_TSS_inshiji)
    index_common2 <- index_E_inshiji[index_common]
    common_E <-Eshijisitequan[index_common2,] 
    common_TSS <-TSSshijisitequan[index_common2,]
    
    if(nrow(common_E) > 0){
      
      index_site_E <- match(common_E$E_sequence_name,fimo_E_new$sequence_name)
      index_E <- which(!is.na(index_site_E))
      
      
      
      index_site_TSS <- match(common_TSS$TSS_sequence_name,fimo_TSS_new$sequence_name)
      index_TSS <- which(!is.na(index_site_TSS))
      
      
      common_E$E_kaishi <- NA
      common_E$E_kaishi[index_E] <- fimo_E_new$start[index_site_E[index_E]]
      
      
      
      common_E$E_jieshu <- NA
      common_E$E_jieshu[index_E] <- fimo_E_new$stop[index_site_E[index_E]]
      
      
      common_TSS$TSS_kaishi <- NA
      common_TSS$TSS_kaishi[index_TSS] <- fimo_TSS_new$start[index_site_TSS[index_TSS]]
      
      
      common_TSS$TSS_jieshu <- NA
      common_TSS$TSS_jieshu[index_TSS] <- fimo_TSS_new$stop[index_site_TSS[index_TSS]]
      
      
      
      pair_site <- cbind(common_E$E_chr,common_E$E_kaishi,common_E$E_jieshu,common_TSS$TSS_chr,common_TSS$TSS_kaishi,common_TSS$TSS_jieshu)
      pair_site <- as.data.frame(pair_site)
      colnames(pair_site) <- c("E_chr","E_start","E_end","TSS_chr","TSS_start","TSS_end")
      pair_site$E_start <- as.numeric(pair_site$E_start)
      pair_site$E_end <- as.numeric(pair_site$E_end)
      pair_site$TSS_start <- as.numeric(pair_site$TSS_start)
      pair_site$TSS_end <- as.numeric(pair_site$TSS_end)
      return(pair_site)
    }
    
    
    
  }
  
  library(foreach)
  library(doParallel)
  cl <- makeCluster(38)
  registerDoParallel(cl)
  
  distance_pair2 <- foreach(x=1:length(motif_name_E)) %dopar% distance_pair(x)
  
  stopCluster(cl)
  
  distance_pair_index <- vector()
  
  for(i in 1:length(distance_pair2)){
   
    num <- nrow(distance_pair2[[i]])
   
    if(is.null(num)){
      distance_pair_index <- c(distance_pair_index,i)
    }
    
  }
  
  distance_pair_index2 <- c(1:length(distance_pair2))
  distance_pair_index2 <- distance_pair_index2[-distance_pair_index]
  
  distance_pair3 <- distance_pair2[distance_pair_index2]
  
  
  motif_name_E_new <- motif_name_E[distance_pair_index2]
  motif_name_TSS_new <- motif_name_TSS[distance_pair_index2]
  
  
  pair_name <- vector()
  for(i in 1:length(motif_name_E_new)){
    pair_name2 <- paste(motif_name_E_new[i],motif_name_TSS_new[i],sep="-")
    pair_name <- c(pair_name,pair_name2)
  }
  
  names(distance_pair3) <- pair_name
  
  distance_num <- vector()
  for(i in 1:length(distance_pair3)){
    distance_num <- c(distance_num,nrow(distance_pair3[[i]]))
  }
  
  names(PWM_E) <- motif_name_E
  names(PWM_TSS) <- motif_name_TSS
  
  index_pair <- vector()
  for(i in 1:length(motif_name_E_new)){
    index <- which(names(PWM_E)==motif_name_E_new[i] & names(PWM_TSS)==motif_name_TSS_new[i])
    index_pair <- c(index_pair,index)
  }
  
  PWM_E_zuizhong <- PWM_E[index_pair]
  PWM_TSS_zuizhong <- PWM_TSS[index_pair]
  
  return(list(PWM_E_zuizhong,PWM_TSS_zuizhong))
  
}