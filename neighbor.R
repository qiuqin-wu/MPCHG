neighbor <- function(E_core,E_nocore,TSS_core,TSS_nocore,Eshijixulie_num_sum,TSSshijixulie_num_sum,E_core_nocore,kmer_num_sum,kmer,kmer_num,kmer_rev_complementary_num_sum,
                     kmer_rev_complementary, kmer_rev_complementary_num,TSS_core_nocore,E_kmer_in_xuliehao,TSS_kmer_in_xuliehao){
  
  
  E_core_kmer_num <- nrow(E_core)
  E_nocore_kmer_num <- nrow(E_nocore)
  
  TSS_core_kmer_num <- nrow(TSS_core)
  TSS_nocore_kmer_num <- nrow(TSS_nocore)
  
  Eshijixulie_num_sum_pair <- Eshijixulie_num_sum
  TSSshijixulie_num_sum_pair <- TSSshijixulie_num_sum
  
  index_site_E <- match(E_core_nocore$kmer, kmer_num_sum$kmer_num_sum)
  E_core_nocore_kmer <- kmer[index_site_E,]
  E_core_nocore_kmer_num <- kmer_num[index_site_E,]
  
  index_site_rev_E <- match(E_core_nocore$kmer_rev, kmer_rev_complementary_num_sum$kmer_rev_complementary_num_sum)
  E_core_nocore_kmer_rev_complementary <- kmer_rev_complementary[index_site_rev_E,]
  E_core_nocore_kmer_rev_complementary_num <- kmer_rev_complementary_num[index_site_rev_E,]
  
  index_site_TSS <- match(TSS_core_nocore$kmer, kmer_num_sum$kmer_num_sum)
  TSS_core_nocore_kmer <- kmer[index_site_TSS,]
  TSS_core_nocore_kmer_num <- kmer_num[index_site_TSS,]
  index_site_rev_TSS <- match(TSS_core_nocore$kmer_rev, kmer_rev_complementary_num_sum$kmer_rev_complementary_num_sum)
  TSS_core_nocore_kmer_rev_complementary <- kmer_rev_complementary[index_site_rev_TSS,]
  TSS_core_nocore_kmer_rev_complementary_num <- kmer_rev_complementary_num[index_site_rev_TSS,]
  
  
  
  # compute E TSS single side neighbor mismatch overlap
  
  library(foreach)
  library(doParallel)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  k <- 6
  
  fun_Elinyu <- function(x){
    linyu_site_E <- rep(0,nrow(E_core_nocore_kmer_num))
    for(j in 1:nrow(E_core_nocore_kmer_num)){
      s <- 0
      qian4 <- 0
      qian5 <- 0
      hou4 <- 0
      hou5 <- 0
      if(j != x){
        for(z in 1:k){
          if(E_core_nocore_kmer_num[x,z]== E_core_nocore_kmer_num[j,z]){
            s=s+1
          }
          
        }
        if(s==5){
          linyu_site_E[j] <- 1
        }
        if(E_core_nocore_kmer_num[j,3]==E_core_nocore_kmer_num[x,1] & E_core_nocore_kmer_num[j,4]==E_core_nocore_kmer_num[x,2] & E_core_nocore_kmer_num[j,5]==E_core_nocore_kmer_num[x,3] & E_core_nocore_kmer_num[j,6]==E_core_nocore_kmer_num[x,4]){
          linyu_site_E[j] <- 2
        }
        if(E_core_nocore_kmer_num[j,2]==E_core_nocore_kmer_num[x,1] & E_core_nocore_kmer_num[j,3]==E_core_nocore_kmer_num[x,2] & E_core_nocore_kmer_num[j,4]==E_core_nocore_kmer_num[x,3] & E_core_nocore_kmer_num[j,5]==E_core_nocore_kmer_num[x,4] & E_core_nocore_kmer_num[j,6]==E_core_nocore_kmer_num[x,5]){
          linyu_site_E[j] <- 3
        }
        if(E_core_nocore_kmer_num[j,1]==E_core_nocore_kmer_num[x,3] & E_core_nocore_kmer_num[j,2]==E_core_nocore_kmer_num[x,4] & E_core_nocore_kmer_num[j,3]==E_core_nocore_kmer_num[x,5] & E_core_nocore_kmer_num[j,4]==E_core_nocore_kmer_num[x,6]){
          linyu_site_E[j] <- 4
        }
        if(E_core_nocore_kmer_num[j,1]==E_core_nocore_kmer_num[x,2] & E_core_nocore_kmer_num[j,2]==E_core_nocore_kmer_num[x,3] & E_core_nocore_kmer_num[j,3]==E_core_nocore_kmer_num[x,4] & E_core_nocore_kmer_num[j,4]==E_core_nocore_kmer_num[x,5] & E_core_nocore_kmer_num[j,5]==E_core_nocore_kmer_num[x,6]){
          linyu_site_E[j] <- 5
        }
        
      }
      
      
    }
    return(linyu_site_E)
  }
    
    fun_TSSlinyu <- function(x){
      linyu_site_TSS <- rep(0,nrow(TSS_core_nocore_kmer_num))
      for(j in 1:nrow(TSS_core_nocore_kmer_num)){
        s <- 0
        qian4 <- 0
        qian5 <- 0
        hou4 <- 0
        hou5 <- 0
        if(j != x){
          for(z in 1:k){
            if(TSS_core_nocore_kmer_num[x,z]== TSS_core_nocore_kmer_num[j,z]){
              s=s+1
            }
            
          }
          if(s==5){
            linyu_site_TSS[j] <- 1
          }
          if(TSS_core_nocore_kmer_num[j,3]==TSS_core_nocore_kmer_num[x,1] & TSS_core_nocore_kmer_num[j,4]==TSS_core_nocore_kmer_num[x,2] & TSS_core_nocore_kmer_num[j,5]==TSS_core_nocore_kmer_num[x,3] & TSS_core_nocore_kmer_num[j,6]==TSS_core_nocore_kmer_num[x,4]){
            linyu_site_TSS[j] <- 2
          }
          if(TSS_core_nocore_kmer_num[j,2]==TSS_core_nocore_kmer_num[x,1] & TSS_core_nocore_kmer_num[j,3]==TSS_core_nocore_kmer_num[x,2] & TSS_core_nocore_kmer_num[j,4]==TSS_core_nocore_kmer_num[x,3] & TSS_core_nocore_kmer_num[j,5]==TSS_core_nocore_kmer_num[x,4] & TSS_core_nocore_kmer_num[j,6]==TSS_core_nocore_kmer_num[x,5]){
            linyu_site_TSS[j] <- 3
          }
          if(TSS_core_nocore_kmer_num[j,1]==TSS_core_nocore_kmer_num[x,3] & TSS_core_nocore_kmer_num[j,2]==TSS_core_nocore_kmer_num[x,4] & TSS_core_nocore_kmer_num[j,3]==TSS_core_nocore_kmer_num[x,5] & TSS_core_nocore_kmer_num[j,4]==TSS_core_nocore_kmer_num[x,6]){
            linyu_site_TSS[j] <- 4
          }
          if(TSS_core_nocore_kmer_num[j,1]==TSS_core_nocore_kmer_num[x,2] & TSS_core_nocore_kmer_num[j,2]==TSS_core_nocore_kmer_num[x,3] & TSS_core_nocore_kmer_num[j,3]==TSS_core_nocore_kmer_num[x,4] & TSS_core_nocore_kmer_num[j,4]==TSS_core_nocore_kmer_num[x,5] & TSS_core_nocore_kmer_num[j,5]==TSS_core_nocore_kmer_num[x,6]){
            linyu_site_TSS[j] <- 5
          }
          
        }
        
        
      }
      return(linyu_site_TSS)
    }
    
    E_linyu <- foreach(x=1:nrow(E_core_nocore_kmer_num)) %dopar% fun_Elinyu(x)
    TSS_linyu <- foreach(x=1:nrow(TSS_core_nocore_kmer_num)) %dopar% fun_TSSlinyu(x)
    
    #别忘了结束并行
    stopCluster(cl)
    
    c_E <- rep(0,nrow(E_core_nocore_kmer_num))
    for(i in 1:nrow(E_core_nocore_kmer_num)){
      c_E[i] <- length(which(E_linyu[[i]] != 0))
    }
    
    c_TSS <- rep(0,nrow(TSS_core_nocore_kmer_num))
    for(i in 1:nrow(TSS_core_nocore_kmer_num)){
      c_TSS[i] <- length(which(TSS_linyu[[i]] != 0))
    }
    
    # compute pair num between E and TSS, and site in sequence
    
    E_TSS_pairnum <- as.data.frame(matrix(0,nrow=(nrow(E_core_nocore_kmer_num) * nrow(TSS_core_nocore_kmer_num)),ncol=3))
    colnames(E_TSS_pairnum) <- c("E","TSS","num")
    
    E <- rep(1:nrow(E_core_nocore_kmer_num),each=nrow(TSS_core_nocore_kmer_num))
    TSS <- rep(1:nrow(TSS_core_nocore_kmer_num),nrow(E_core_nocore_kmer_num))
    
    E_TSS_pairnum$E <- E
    E_TSS_pairnum$TSS <- TSS
    
    numm <- as.numeric(vector()) 
    pairnum <- lapply(1:nrow(E_core_nocore_kmer_num),function(x){
      for(j in 1:nrow(TSS_core_nocore_kmer_num)){
        intersection <- intersect(E_kmer_in_xuliehao[[x]],TSS_kmer_in_xuliehao[[j]])
        numm <- c(numm,length(intersection)) 
      }
      numa <- numm
      
    })
    
    pairnum2 <- unlist(pairnum)
    
    E_TSS_pairnum$num <- pairnum2
    
    E_TSS_pairnum$new_weight <- rep(0,nrow(E_TSS_pairnum))
    
    diff  <- max(E_TSS_pairnum$num)-min(E_TSS_pairnum$num)
    
    library(foreach)
    library(doParallel)
    cl <- makeCluster(4)
    registerDoParallel(cl)
    
    weight <- function(x){
      
      lin_e <- which(E_linyu[[x]] != 0)
      index <- which(E_TSS_pairnum$E==1)
      lin_TSS <- as.numeric(vector())
      for(h in 1:length(index)){
        lin <- which(TSS_linyu[[h]] != 0)
        lin_TSS <- c(lin_TSS,length(lin))
        
      }
      lin_TSS <- rep(length(lin_e),nrow(TSS_core_nocore))
      weight_topo <- (lin_TSS / nrow(E_core)) + (lin_TSS / nrow(TSS_core))
      
      index2 <- which(E_TSS_pairnum$E==x)
      
      
      weight_true <- (E_TSS_pairnum$num[index2] - min(E_TSS_pairnum$num))/diff
      weight_new <- weight_topo + weight_true
      
      
      return(weight_new)
      
    }
  
    weight_new <- foreach(x=1:nrow(E_core_nocore)) %dopar% weight(x)
    stopCluster(cl)
    
    weight_new2 <- unlist(weight_new)
    E_TSS_pairnum$new_weight <- weight_new2
    E_TSS_pairnum_newweight <- E_TSS_pairnum
    
    E_TSS_pairnum_core_index <- which(E_TSS_pairnum_newweight$E <= E_core_kmer_num & E_TSS_pairnum_newweight$TSS <= TSS_core_kmer_num)
    E_TSS_pairnum_core <- E_TSS_pairnum_newweight[E_TSS_pairnum_core_index,]
    E_TSS_pairnum_nocore <- E_TSS_pairnum_newweight[-E_TSS_pairnum_core_index,]
    
    E_TSS_pairnum_core <- E_TSS_pairnum_core[order(-E_TSS_pairnum_core$new_weight),]
    E_TSS_pairnum_nocore <- E_TSS_pairnum_nocore[order(-E_TSS_pairnum_nocore$new_weight),]
    
    E_TSS_pairnum_core_nocore <- rbind(E_TSS_pairnum_core,E_TSS_pairnum_nocore)
    
    E_TSS_pairnum_core_num <- nrow(E_TSS_pairnum_core)
    E_TSS_pairnum_nocore_num <- nrow(E_TSS_pairnum_nocore)
    E_TSS_pairnum_core_nocore$standweight <- (E_TSS_pairnum_core_nocore$new_weight - min(E_TSS_pairnum_core_nocore$new_weight)) / (max(E_TSS_pairnum_core_nocore$new_weight - min(E_TSS_pairnum_core_nocore$new_weight)))
    whole_edg_num <- length(which(E_TSS_pairnum_core_nocore$standweight != 0))
    whole_weight <- sum(E_TSS_pairnum_core_nocore$standweight)
    whole_avg <- whole_weight / whole_edg_num
    
  
  
  return(list(E_core_nocore_kmer,E_core_nocore_kmer_num,E_core_nocore_kmer_rev_complementary,E_core_nocore_kmer_rev_complementary_num,
              TSS_core_nocore_kmer,TSS_core_nocore_kmer_num,TSS_core_nocore_kmer_rev_complementary,TSS_core_nocore_kmer_rev_complementary_num,
              E_linyu,TSS_linyu,weight_new,E_TSS_pairnum,E_TSS_pairnum_core,E_TSS_pairnum_nocore,E_TSS_pairnum_core_nocore))
  
  
}