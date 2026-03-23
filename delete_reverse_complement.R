delete_reverse_complement <- function(TSSbeijing_kmer_freq,Eshiji_kmer_freq,TSSshiji_kmer_freq,Ebeijing_kmer_freq,kmer){
  
  # delete fan xaing hu bu TSSi zi shen de kmer
  index <- which(TSSbeijing_kmer_freq$kmer==TSSbeijing_kmer_freq$kmer_rev)
  TSSshiji_kmer_freq[index,5] <- TSSshiji_kmer_freq[index,5]/2
  TSSbeijing_kmer_freq[index,5] <- TSSbeijing_kmer_freq[index,5]/2
  Eshiji_kmer_freq[index,5] <- Eshiji_kmer_freq[index,5]/2
  Ebeijing_kmer_freq[index,5] <- Ebeijing_kmer_freq[index,5]/2
  output_file_path_name <- paste(output_folder_path,"TSSbeijing_kmer_freq.RDS",sep="/")
  saveRDS(TSSbeijing_kmer_freq,output_file_path_name)
  output_file_path_name <- paste(output_folder_path,"TSSshiji_kmer_freq.RDS",sep="/")
  saveRDS(TSSshiji_kmer_freq,output_file_path_name)
  output_file_path_name <- paste(output_folder_path,"Ebeijing_kmer_freq.RDS",sep="/")
  saveRDS(Ebeijing_kmer_freq,output_file_path_name)
  output_file_path_name <- paste(output_folder_path,"Eshiji_kmer_freq.RDS",sep="/")
  saveRDS(Eshiji_kmer_freq,output_file_path_name)
  
  
  index_hubu <- rep(0,time=1,each=nrow(Eshiji_kmer_freq))
  
  for(i in 2 :nrow(Eshiji_kmer_freq)){
    index <- Eshiji_kmer_freq$kmer_rev[i] %in% Eshiji_kmer_freq$kmer[1:(i-1)]
    if(index==TRUE){
      index_hubu[i] <- 1
    }
  }
  index <- which(index_hubu==0)  
  
  Eshiji_kmer2078_freq <- Eshiji_kmer_freq[index,]
  Ebeijing_kmer2078_freq <- Ebeijing_kmer_freq[index,]
  TSSshiji_kmer2078_freq <- TSSshiji_kmer_freq[index,]
  TSSbeijing_kmer2078_freq <- TSSbeijing_kmer_freq[index,]
  
  
  kmer_2078 <- kmer[index,]
  Eshiji_kmer2078_freq_mean <- mean(Eshiji_kmer2078_freq$sum)
  Ebeijing_kmer2078_freq_mean <- mean(Ebeijing_kmer2078_freq$sum)
  TSSshiji_kmer2078_freq_mean <- mean(TSSshiji_kmer2078_freq$sum)
  TSSbeijing_kmer2078_freq_mean <- mean(TSSbeijing_kmer2078_freq$sum)
  
  
  Eshiji_kmer2078_freq_sd <- sd(Eshiji_kmer2078_freq$sum)
  Ebeijing_kmer2078_freq_sd <- sd(Ebeijing_kmer2078_freq$sum)
  TSSshiji_kmer2078_freq_sd <- sd(TSSshiji_kmer2078_freq$sum)
  TSSbeijing_kmer2078_freq_sd <- sd(TSSbeijing_kmer2078_freq$sum)
  
  index1 <- which(Eshiji_kmer2078_freq$sum > Eshiji_kmer2078_freq_mean - Eshiji_kmer2078_freq_sd)
  index2 <- which(TSSshiji_kmer2078_freq$sum > TSSshiji_kmer2078_freq_mean - TSSshiji_kmer2078_freq_sd)
  
  Eshiji_kmer_freq_remain <- Eshiji_kmer2078_freq[index1,]
  TSSshiji_kmer_freq_remain <- TSSshiji_kmer2078_freq[index2,]
  
  index_Ebeijing <- match(Eshiji_kmer_freq_remain$kmer,Ebeijing_kmer2078_freq$kmer)
  index_TSSbeijing <- match(TSSshiji_kmer_freq_remain$kmer,TSSbeijing_kmer2078_freq$kmer)
  
  Ebeijing_kmer_freq_remain <- Ebeijing_kmer2078_freq[index_Ebeijing,]
  TSSbeijing_kmer_freq_remain <- TSSbeijing_kmer2078_freq[index_TSSbeijing,]
  
  pro_E_shiji <- as.data.frame(Eshiji_kmer_freq_remain$sum / sum(Eshiji_kmer_freq_remain$sum))
  pro_E_beijing <-  as.data.frame(Ebeijing_kmer_freq_remain$sum / sum(Ebeijing_kmer_freq_remain$sum))
  pro_TSS_shiji <-  as.data.frame(TSSshiji_kmer_freq_remain$sum / sum(TSSshiji_kmer_freq_remain$sum))
  pro_TSS_beijing <-  as.data.frame(TSSbeijing_kmer_freq_remain$sum / sum(TSSbeijing_kmer_freq_remain$sum))
  
  Eshiji_kmer_freq_remain <- cbind(Eshiji_kmer_freq_remain,pro_E_shiji)
  colnames(Eshiji_kmer_freq_remain)[6] <- "pro_E_shiji"
  Ebeijing_kmer_freq_remain <- cbind(Ebeijing_kmer_freq_remain,pro_E_beijing)
  colnames(Ebeijing_kmer_freq_remain)[6] <- "pro_E_beijing"
  TSSshiji_kmer_freq_remain <- cbind(TSSshiji_kmer_freq_remain,pro_TSS_shiji)
  colnames(TSSshiji_kmer_freq_remain)[6] <- "pro_TSS_shiji"
  TSSbeijing_kmer_freq_remain <- cbind(TSSbeijing_kmer_freq_remain,pro_TSS_beijing)
  colnames(TSSbeijing_kmer_freq_remain)[6] <- "pro_TSS_beijing"
  
  p_E_shiji_beijing <- (Eshiji_kmer_freq_remain$sum + Ebeijing_kmer_freq_remain$sum) / (sum(Eshiji_kmer_freq_remain$sum) + sum(Ebeijing_kmer_freq_remain$sum))
  p_TSS_shiji_beijing <- (TSSshiji_kmer_freq_remain$sum + TSSbeijing_kmer_freq_remain$sum) / (sum(TSSshiji_kmer_freq_remain$sum) + sum(TSSbeijing_kmer_freq_remain$sum))
  
  
  
  yijian_p_E_shiji_beijing <- rep(1,time=1,each=length(p_E_shiji_beijing))
  yijian_p_TSS_shiji_beijing <- rep(1,time=1,each=length(p_TSS_shiji_beijing))
  
  z_E_shiji_beijing <- (Eshiji_kmer_freq_remain$pro_E_shiji - Ebeijing_kmer_freq_remain$pro_E_beijing) / sqrt(((p_E_shiji_beijing * (yijian_p_E_shiji_beijing - p_E_shiji_beijing)) *((yijian_p_E_shiji_beijing / sum(Eshiji_kmer_freq_remain$sum)) + ((yijian_p_E_shiji_beijing / sum(Ebeijing_kmer_freq_remain$sum))))))
  z_TSS_shiji_beijing <- (TSSshiji_kmer_freq_remain$pro_TSS_shiji - TSSbeijing_kmer_freq_remain$pro_TSS_beijing) / sqrt(((p_TSS_shiji_beijing * (yijian_p_TSS_shiji_beijing - p_TSS_shiji_beijing)) *((yijian_p_TSS_shiji_beijing / sum(TSSshiji_kmer_freq_remain$sum)) + ((yijian_p_TSS_shiji_beijing / sum(TSSbeijing_kmer_freq_remain$sum))))))
  
  
  p_E_whole <- as.data.frame(matrix(nrow=nrow(Eshiji_kmer_freq_remain),ncol=4))
  p_E_whole[,1] <- Eshiji_kmer_freq_remain$pro_E_shiji
  p_E_whole[,2] <- Ebeijing_kmer_freq_remain$pro_E_beijing
  p_E_whole[,3] <- p_E_shiji_beijing
  p_E_whole[,4] <- z_E_shiji_beijing
  colnames(p_E_whole) <- c("pro_E_shiji","pro_E_beijing","p_E_shiji_beijing","z_E_shiji_beijing")
  
  p_TSS_whole <- as.data.frame(matrix(nrow=nrow(TSSshiji_kmer_freq_remain),ncol=4))
  p_TSS_whole[,1] <- TSSshiji_kmer_freq_remain$pro_TSS_shiji
  p_TSS_whole[,2] <- TSSbeijing_kmer_freq_remain$pro_TSS_beijing
  p_TSS_whole[,3] <- p_TSS_shiji_beijing
  p_TSS_whole[,4] <- z_TSS_shiji_beijing
  colnames(p_TSS_whole) <- c("pro_TSS_shiji","pro_TSS_beijing","p_TSS_shiji_beijing","z_TSS_shiji_beijing")
  
  kmer_E_index <- Eshiji_kmer_freq_remain[,c(1,3)]
  kmer_TSS_index <- TSSshiji_kmer_freq_remain[,c(1,3)]
  
  p_E_whole <- cbind(kmer_E_index,p_E_whole)
  p_E_whole <- cbind(p_E_whole,Eshiji_kmer_freq_remain$sum)
  colnames(p_E_whole)[7] <- "sum"
  p_TSS_whole <- cbind(kmer_TSS_index,p_TSS_whole)
  p_TSS_whole <- cbind(p_TSS_whole,TSSshiji_kmer_freq_remain$sum)
  colnames(p_TSS_whole)[7] <- "sum"
  
  #yuzhi z > 1.96zuoTSSi he xin ji he
  z_E_core_index <- which(p_E_whole$z_E_shiji_beijing > 1.96)
  z_TSS_core_index <- which(p_TSS_whole$z_TSS_shiji_beijing > 1.96)
  
  #重新排序kmer以及对应的频数和zscore
  E_core <- p_E_whole[z_E_core_index,]
  E_core <- E_core[order(-E_core$sum),]
  
  E_nocore <- p_E_whole[-(z_E_core_index),]
  E_nocore <- E_nocore[order(-E_nocore$sum),]
  
  E_core_nocore <- rbind(E_core,E_nocore)
  E_core_nocore <- E_core_nocore[order(-E_core_nocore[,6]),]
  E_core <- E_core_nocore[which(E_core_nocore$z_E_shiji_beijing > 1.96),]
  E_nocore <- E_core_nocore[which(E_core_nocore$z_E_shiji_beijing <= 1.96),]
  
  E_core_kmer_num <- nrow(E_core)
  E_nocore_kmer_num <- nrow(E_nocore)
  
  TSS_core <- p_TSS_whole[z_TSS_core_index,]
  TSS_core <- TSS_core[order(-TSS_core$sum),]
  TSS_nocore <- p_TSS_whole[-(z_TSS_core_index),]
  TSS_nocore <- TSS_nocore[order(-TSS_nocore$sum),]
  TSS_core_nocore <- rbind(TSS_core,TSS_nocore)
  
  TSS_core_nocore <- TSS_core_nocore[order(-TSS_core_nocore[,6]),]
  
  TSS_core <- TSS_core_nocore[which(TSS_core_nocore$z_TSS_shiji_beijing > 1.96),]
  TSS_nocore <- TSS_core_nocore[which(TSS_core_nocore$z_TSS_shiji_beijing <= 1.96),]
  
  TSS_core_kmer_num <- nrow(TSS_core)
  TSS_nocore_kmer_num <- nrow(TSS_nocore)
  
  
  
  return(list(Eshiji_kmer_freq_remain,Ebeijing_kmer_freq_remain,TSSshiji_kmer_freq_remain,TSSbeijing_kmer_freq_remain,p_E_whole,p_TSS_whole,
              E_core,E_nocore,E_core_nocore,TSS_core,TSS_nocore,TSS_core_nocore))
  
  
  
  
  
}