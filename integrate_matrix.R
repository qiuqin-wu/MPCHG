integrate_matrix <- function(E_union2_0.4,E_core_nocore,TSS_union2_0.4,TSS_core_nocore,E_linyu,E_core_nocore_kmer,TSS_core_nocore_kmer){
  
  E_union3_0.4 <- E_union2_0.4  #E_union3_0.5 为待合并的kmer，每个元素表示Kmer的四进制数
  
  for(i in 1:length(E_union3_0.4)){
    E_union3_0.4[[i]] <- E_core_nocore$kmer[E_union2_0.4[[i]]]
  }
  
  TSS_union3_0.4 <- TSS_union2_0.4  #TSS_union3_0.5 为待合并的kmer，每个元素表示Kmer的四进制数
  for(i in 1:length(TSS_union3_0.4)){
    TSS_union3_0.4[[i]] <- TSS_core_nocore$kmer[TSS_union2_0.4[[i]]]
  }
  #待合并的kmer的数量
  E_union2_0.4_kmer_sum <- lapply(1:length(E_union2_0.4),function(x){
    index <- E_core_nocore$sum[E_union2_0.4[[x]]]
  })
  TSS_union2_0.4_kmer_sum <- lapply(1:length(TSS_union2_0.4),function(x){
    index <- TSS_core_nocore$sum[TSS_union2_0.4[[x]]]
  })
  
  #先判断Kmer与第一个是什么关系，若没关系比较与第二个什么关系,但结果显示都与第一个kmer有关系
  hhh1 <- vector()
  hhh2 <- vector()
  hhh3 <- vector()
  for(i in 1:length(E_union2_0.4)){
    guanxi <- as.vector(matrix(0,nrow=1,ncol=length(E_union2_0.4[[i]])-1))
    for(j in 1:length(E_union2_0.4[[i]])){
      diyi <- E_union2_0.4[[i]][1]
      index <- which(E_linyu[[diyi]] > 0)
      #index_site <- which(E_union2_0.5[[i]] %in% index)
      jiaoji <- intersect(index,E_union2_0.4[[i]][1:length(E_union2_0.4[[i]])])
      guanxi <- E_linyu[[diyi]][jiaoji]
      
      #cat("length(guanxi):",length(guanxi),"\n")
      #cat("length(E_union2_0.4[[i]]):",length(E_union2_0.4[[i]]),"\n")
      hhh1 <- c(hhh1,length(guanxi))
      hhh2 <- c(hhh2,length(E_union2_0.4[[i]]))
      hhh7 <- hhh2-hhh1
      hhh4 <- length(which(hhh7 > 1))
    }
    hhh3 <- c(hhh3,hhh4)
  }
  
  for(i in 1:length(TSS_union2_0.4)){
    #guanxi <- as.vector(matrix(0,nrow=1,ncol=length(E_union2_0.5[[i]])-1))
    for(j in 1:length(TSS_union2_0.4[[i]])){
      diyi <- TSS_union2_0.4[[i]][1]
      index <- which(TSS_linyu[[diyi]] > 0)
      #index_site <- which(E_union2_0.5[[i]] %in% index)
      jiaoji <- intersect(index,TSS_union2_0.4[[i]][1:length(TSS_union2_0.4[[i]])])
      guanxi <- TSS_linyu[[diyi]][jiaoji]
      
      #  cat("length(guanxi):",length(guanxi),"\n")
      #  cat("length(TSS_union2_0.5[[i]]):",length(TSS_union2_0.5[[i]]),"\n")
      #  hhh1 <- c(hhh1,length(guanxi))
      #  hhh2 <- c(hhh2,length(TSS_union2_0.5[[i]]))
    }
  }
  
  #先全部扩充为长度为10的向量，扩充的字母全部用B代替
  
  juzhen_E <- lapply(1:length(E_union2_0.4),function(x){
    kmer_dataframe <- E_core_nocore_kmer[E_union2_0.4[[x]],]
    matrix_10 <- as.data.frame(matrix("B",nrow=length(E_union2_0.4[[x]]),ncol=10))
    matrix_10[1,3:8] <-  kmer_dataframe[1,]
    index <- E_union2_0.4[[x]]
    index <- index[-1]
    linyu_index <- E_linyu[[E_union2_0.4[[x]][1]]][index]
    for(j in 1:length(linyu_index)){
      if(linyu_index[j] == 1){
        matrix_10[j+1,3:8] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 2){
        matrix_10[j+1,1:6] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 3){
        matrix_10[j+1,2:7] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 4){
        matrix_10[j+1,5:10] <- kmer_dataframe[j+1,]
      }
      else {
        matrix_10[j+1,4:9] <- kmer_dataframe[j+1,]
      }
    }
    
    matrix2_10 <- matrix_10
  })
  
  
  juzhen_TSS <- lapply(1:length(TSS_union2_0.4),function(x){
    kmer_dataframe <- TSS_core_nocore_kmer[TSS_union2_0.4[[x]],]
    matrix_10 <- as.data.frame(matrix("B",nrow=length(TSS_union2_0.4[[x]]),ncol=10))
    matrix_10[1,3:8] <-  kmer_dataframe[1,]
    index <- TSS_union2_0.4[[x]]
    index <- index[-1]
    linyu_index <- TSS_linyu[[TSS_union2_0.4[[x]][1]]][index]
    for(j in 1:length(linyu_index)){
      if(linyu_index[j] == 1){
        matrix_10[j+1,3:8] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 2){
        matrix_10[j+1,1:6] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 3){
        matrix_10[j+1,2:7] <- kmer_dataframe[j+1,]
      }
      else if(linyu_index[j] == 4){
        matrix_10[j+1,5:10] <- kmer_dataframe[j+1,]
      }
      else {
        matrix_10[j+1,4:9] <- kmer_dataframe[j+1,]
      }
    }
    
    matrix2_10 <- matrix_10
  })
  
  PWM_E <- lapply(1:length(juzhen_E),function(x){
    matrix_E <- matrix(0,nrow=4,ncol=10)
    rownames(matrix_E) <- c("A","C","G","T")
    for(i in 1:10){
      index_A <- which(juzhen_E[[x]][,i] == 'A')
      index_T <- which(juzhen_E[[x]][,i] == 'T')
      index_C <- which(juzhen_E[[x]][,i] == 'C')
      index_G <- which(juzhen_E[[x]][,i] == 'G')
      sum_A <- sum(E_union2_0.4_kmer_sum[[x]][index_A])
      sum_T <- sum(E_union2_0.4_kmer_sum[[x]][index_T])
      sum_C <- sum(E_union2_0.4_kmer_sum[[x]][index_C])
      sum_G <- sum(E_union2_0.4_kmer_sum[[x]][index_G])
      if(sum_A > 0 | sum_C > 0 | sum_G > 0 | sum_T > 0){
        sum_whole <- sum(sum_A,sum_T,sum_C,sum_G)
        avg_A <- sum_A/sum_whole
        avg_T <- sum_T/sum_whole
        avg_C <- sum_C/sum_whole
        avg_G <- sum_G/sum_whole
        matrix_E[1,i] <- as.numeric(sprintf("%0.4f", avg_A))
        matrix_E[2,i] <- as.numeric(sprintf("%0.4f", avg_C))
        matrix_E[3,i] <- as.numeric(sprintf("%0.4f", avg_G))
        matrix_E[4,i] <- as.numeric(sprintf("%0.4f", avg_T))
      }
      else{
        matrix_E[,i] <- 0
      }
    }
    matrix2_E <- t(matrix_E)
    
  })
  
  PWM_TSS <- lapply(1:length(juzhen_TSS),function(x){
    matrix_TSS <- matrix(0,nrow=4,ncol=10)
    rownames(matrix_TSS) <- c("A","C","G","T")
    for(i in 1:10){
      index_A <- which(juzhen_TSS[[x]][,i] == 'A')
      index_T <- which(juzhen_TSS[[x]][,i] == 'T')
      index_C <- which(juzhen_TSS[[x]][,i] == 'C')
      index_G <- which(juzhen_TSS[[x]][,i] == 'G')
      sum_A <- sum(TSS_union2_0.4_kmer_sum[[x]][index_A])
      sum_T <- sum(TSS_union2_0.4_kmer_sum[[x]][index_T])
      sum_C <- sum(TSS_union2_0.4_kmer_sum[[x]][index_C])
      sum_G <- sum(TSS_union2_0.4_kmer_sum[[x]][index_G])
      if(sum_A > 0 | sum_C > 0 | sum_G > 0 | sum_T > 0){
        sum_whole <- sum(sum_A,sum_T,sum_C,sum_G)
        avg_A <- sum_A/sum_whole
        avg_T <- sum_T/sum_whole
        avg_C <- sum_C/sum_whole
        avg_G <- sum_G/sum_whole
        matrix_TSS[1,i] <- as.numeric(sprintf("%0.4f", avg_A))
        matrix_TSS[2,i] <- as.numeric(sprintf("%0.4f", avg_C))
        matrix_TSS[3,i] <- as.numeric(sprintf("%0.4f", avg_G))
        matrix_TSS[4,i] <- as.numeric(sprintf("%0.4f", avg_T))
      }
      else{
        matrix_TSS[,i] <- 0
      }
    }
    matrix2_TSS <- t(matrix_TSS)
    
  })
  
  PWM_E2 <- PWM_E
  for(i in 1:length(PWM_E)){
    index <- rowSums(PWM_E[[i]][,1:4])
    index2 <- which(index == 0)
    if(length(index2) > 0){
      PWM_E2[[i]] <- PWM_E2[[i]][-index2,]
    }
  }
  
  PWM_TSS2 <- PWM_TSS
  for(i in 1:length(PWM_TSS)){
    index <- rowSums(PWM_TSS[[i]][,1:4])
    index2 <- which(index == 0)
    if(length(index2) > 0){
      PWM_TSS2[[i]] <- PWM_TSS2[[i]][-index2,]
    }
  }
  
  PWM_E <- PWM_E2
  PWM_TSS <- PWM_TSS2
  
 return(list(E_union3_0.4,TSS_union3_0.4,E_union2_0.4_kmer_sum,TSS_union2_0.4_kmer_sum,PWM_E,PWM_TSS)) 
}