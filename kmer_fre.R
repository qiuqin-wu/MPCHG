kmer_fre <- function(Eshijixulie_num_sum,kmer_num_sum,kmer_rev_complementary_num_sum,kmer_rev_freq,Ebeijingxulie_num_sum,
                     TSSshijixulie_num_sum,TSSbeijingxulie_num_sum){
  
  #jiang lie biao zhuan hua TSSi yi ge xiang liang
  list_union <- unlist(Eshijixulie_num_sum)
  
  kmer_freq <- lapply(1:nrow(kmer_num_sum),function(x){
    sum <- length(which(list_union == kmer_num_sum[x,1]))
    
  })
  
  kmer_freq <- as.data.frame(unlist(kmer_freq))
  Eshiji_kmer_num_freq <- cbind(kmer_num_sum,kmer_freq)
  colnames(Eshiji_kmer_num_freq) <- c("kmer","kmer_freq")
  
  # compute rev_complementary freq
  kmer_rev_freq <- lapply(1:nrow(kmer_rev_complementary_num_sum),function(x){
    sum <- length(which(list_union == kmer_rev_complementary_num_sum[x,1]))
    
  })
  
  kmer_rev_freq <- as.data.frame(unlist(kmer_rev_freq))
  Eshiji_kmer_rev_num_freq <- cbind(kmer_rev_complementary_num_sum,kmer_rev_freq)
  colnames(Eshiji_kmer_rev_num_freq) <- c("kmer","kmer_rev_freq")
  
  
  
  k3 <- Eshiji_kmer_num_freq$kmer_freq + Eshiji_kmer_rev_num_freq$kmer_rev_freq
  Eshiji_kmer_num_freq  <- cbind(Eshiji_kmer_num_freq,k3)
  Eshiji_kmer_rev_num_freq <- cbind(Eshiji_kmer_rev_num_freq,k3)
  colnames(Eshiji_kmer_num_freq)[3] <- "zheng_fu_chain_sum"
  colnames(Eshiji_kmer_rev_num_freq)[3] <- "zheng_fu_chain_sum"
  Eshiji_kmer_freq <- as.data.frame(matrix(nrow=nrow(Eshiji_kmer_num_freq),ncol=5))
  colnames(Eshiji_kmer_freq) <- c("kmer","kmer_freq","kmer_rev","kmer_rev_freq","sum")
  Eshiji_kmer_freq$kmer <- Eshiji_kmer_num_freq$kmer
  Eshiji_kmer_freq$kmer_freq <- Eshiji_kmer_num_freq$kmer_freq
  Eshiji_kmer_freq$sum <- Eshiji_kmer_num_freq$zheng_fu_chain_sum
  Eshiji_kmer_freq$kmer_rev <- Eshiji_kmer_rev_num_freq$kmer
  Eshiji_kmer_freq$kmer_rev_freq <- Eshiji_kmer_rev_num_freq$kmer_rev_freq
  output_file_path_name <- paste(output_folder_path,"Eshiji_kmer_freq.RDS",sep="/")
  #saveRDS(Eshiji_kmer_freq,output_file_path_name)
  
  
  
  
  list_union2 <- unlist(Ebeijingxulie_num_sum)
  
  kmer_freq2 <- lapply(1:nrow(kmer_num_sum),function(x){
    sum <- length(which(list_union2 == kmer_num_sum[x,1]))
    
  })
  
  
  
  
  
  kmer_freq2 <- as.data.frame(unlist(kmer_freq2))
  Ebeijing_kmer_num_freq <- cbind(kmer_num_sum,kmer_freq2)
  colnames(Ebeijing_kmer_num_freq) <- c("kmer","kmer_freq")
  
  # compute rev_complementary freq
  kmer_rev_freq2 <- lapply(1:nrow(kmer_rev_complementary_num_sum),function(x){
    sum <- length(which(list_union2 == kmer_rev_complementary_num_sum[x,1]))
    
  })
  
  kmer_rev_freq2 <- as.data.frame(unlist(kmer_rev_freq2))
  Ebeijing_kmer_rev_num_freq <- cbind(kmer_rev_complementary_num_sum,kmer_rev_freq2)
  colnames(Ebeijing_kmer_rev_num_freq) <- c("kmer","kmer_rev_freq")
  
  
  k3 <- Ebeijing_kmer_num_freq$kmer_freq + Ebeijing_kmer_rev_num_freq$kmer_rev_freq
  Ebeijing_kmer_num_freq  <- cbind(Ebeijing_kmer_num_freq,k3)
  Ebeijing_kmer_rev_num_freq <- cbind(Ebeijing_kmer_rev_num_freq,k3)
  colnames(Ebeijing_kmer_num_freq)[3] <- "zheng_fu_chain_sum"
  colnames(Ebeijing_kmer_rev_num_freq)[3] <- "zheng_fu_chain_sum"
  Ebeijing_kmer_freq <- as.data.frame(matrix(nrow=nrow(Ebeijing_kmer_num_freq),ncol=5))
  colnames(Ebeijing_kmer_freq) <- c("kmer","kmer_freq","kmer_rev","kmer_rev_freq","sum")
  Ebeijing_kmer_freq$kmer <- Ebeijing_kmer_num_freq$kmer
  Ebeijing_kmer_freq$kmer_freq <- Ebeijing_kmer_num_freq$kmer_freq
  Ebeijing_kmer_freq$sum <- Ebeijing_kmer_num_freq$zheng_fu_chain_sum
  Ebeijing_kmer_freq$kmer_rev <- Ebeijing_kmer_rev_num_freq$kmer
  Ebeijing_kmer_freq$kmer_rev_freq <- Ebeijing_kmer_rev_num_freq$kmer_rev_freq
  output_file_path_name <- paste(output_folder_path,"Ebeijing_kmer_freq.RDS",sep="/")
  #saveRDS(Ebeijing_kmer_freq,output_file_path_name)
  
  
  
  list_union3 <- unlist(TSSshijixulie_num_sum)
  
  kmer_freq3 <- lapply(1:nrow(kmer_num_sum),function(x){
    sum <- length(which(list_union3 == kmer_num_sum[x,1]))
    
  })
  
  kmer_freq3 <- as.data.frame(unlist(kmer_freq3))
  TSSshiji_kmer_num_freq <- cbind(kmer_num_sum,kmer_freq3)
  colnames(TSSshiji_kmer_num_freq) <- c("kmer","kmer_freq")
  
  
  # compute rev_complementary freq
  kmer_rev_freq3 <- lapply(1:nrow(kmer_rev_complementary_num_sum),function(x){
    sum <- length(which(list_union3 == kmer_rev_complementary_num_sum[x,1]))
    
  })
  
  kmer_rev_freq3 <- as.data.frame(unlist(kmer_rev_freq3))
  TSSshiji_kmer_rev_num_freq <- cbind(kmer_rev_complementary_num_sum,kmer_rev_freq3)
  colnames(TSSshiji_kmer_rev_num_freq) <- c("kmer","kmer_rev_freq")
  
  
  
  
  k3 <- TSSshiji_kmer_num_freq$kmer_freq + TSSshiji_kmer_rev_num_freq$kmer_rev_freq
  TSSshiji_kmer_num_freq  <- cbind(TSSshiji_kmer_num_freq,k3)
  TSSshiji_kmer_rev_num_freq <- cbind(TSSshiji_kmer_rev_num_freq,k3)
  colnames(TSSshiji_kmer_num_freq)[3] <- "zheng_fu_chain_sum"
  colnames(TSSshiji_kmer_rev_num_freq)[3] <- "zheng_fu_chain_sum"
  TSSshiji_kmer_freq <- as.data.frame(matrix(nrow=nrow(TSSshiji_kmer_num_freq),ncol=5))
  colnames(TSSshiji_kmer_freq) <- c("kmer","kmer_freq","kmer_rev","kmer_rev_freq","sum")
  TSSshiji_kmer_freq$kmer <- TSSshiji_kmer_num_freq$kmer
  TSSshiji_kmer_freq$kmer_freq <- TSSshiji_kmer_num_freq$kmer_freq
  TSSshiji_kmer_freq$sum <- TSSshiji_kmer_num_freq$zheng_fu_chain_sum
  TSSshiji_kmer_freq$kmer_rev <- TSSshiji_kmer_rev_num_freq$kmer
  TSSshiji_kmer_freq$kmer_rev_freq <- TSSshiji_kmer_rev_num_freq$kmer_rev_freq
  output_file_path_name <- paste(output_folder_path,"TSSshiji_kmer_freq.RDS",sep="/")
  #saveRDS(TSSshiji_kmer_freq,output_file_path_name)
  
  
  list_union4 <- unlist(TSSbeijingxulie_num_sum)
  
  kmer_freq4 <- lapply(1:nrow(kmer_num_sum),function(x){
    sum <- length(which(list_union4 == kmer_num_sum[x,1]))
    
  })
  
  kmer_freq4 <- as.data.frame(unlist(kmer_freq4))
  TSSbeijing_kmer_num_freq <- cbind(kmer_num_sum,kmer_freq4)
  colnames(TSSbeijing_kmer_num_freq) <- c("kmer","kmer_freq")
  
  
  
  # compute rev_complementary freq
  kmer_rev_freq4 <- lapply(1:nrow(kmer_rev_complementary_num_sum),function(x){
    sum <- length(which(list_union4 == kmer_rev_complementary_num_sum[x,1]))
    
  })
  
  kmer_rev_freq4 <- as.data.frame(unlist(kmer_rev_freq4))
  TSSbeijing_kmer_rev_num_freq <- cbind(kmer_rev_complementary_num_sum,kmer_rev_freq4)
  colnames(TSSbeijing_kmer_rev_num_freq) <- c("kmer","kmer_rev_freq")
  
  k3 <- TSSbeijing_kmer_num_freq$kmer_freq + TSSbeijing_kmer_rev_num_freq$kmer_rev_freq
  TSSbeijing_kmer_num_freq  <- cbind(TSSbeijing_kmer_num_freq,k3)
  TSSbeijing_kmer_rev_num_freq <- cbind(TSSbeijing_kmer_rev_num_freq,k3)
  colnames(TSSbeijing_kmer_num_freq)[3] <- "zheng_fu_chain_sum"
  colnames(TSSbeijing_kmer_rev_num_freq)[3] <- "zheng_fu_chain_sum"
  TSSbeijing_kmer_freq <- as.data.frame(matrix(nrow=nrow(TSSbeijing_kmer_num_freq),ncol=5))
  colnames(TSSbeijing_kmer_freq) <- c("kmer","kmer_freq","kmer_rev","kmer_rev_freq","sum")
  TSSbeijing_kmer_freq$kmer <- TSSbeijing_kmer_num_freq$kmer
  TSSbeijing_kmer_freq$kmer_freq <- TSSbeijing_kmer_num_freq$kmer_freq
  TSSbeijing_kmer_freq$sum <- TSSbeijing_kmer_num_freq$zheng_fu_chain_sum
  TSSbeijing_kmer_freq$kmer_rev <- TSSbeijing_kmer_rev_num_freq$kmer
  TSSbeijing_kmer_freq$kmer_rev_freq <- TSSbeijing_kmer_rev_num_freq$kmer_rev_freq
  output_file_path_name <- paste(output_folder_path,"TSSbeijing_kmer_freq.RDS",sep="/")
  #saveRDS(TSSbeijing_kmer_freq,output_file_path_name)
  
 return(list(Eshiji_kmer_freq,Ebeijing_kmer_freq,TSSshiji_kmer_freq,TSSbeijing_kmer_freq)) 
}