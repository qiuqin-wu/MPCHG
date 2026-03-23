dense_subgraph <- function(E_TSS_pairnum_core_nocore,E_linyu,TSS_linyu){
  
  whole_edg_num <- length(which(E_TSS_pairnum_core_nocore$standweight != 0))
  whole_weight <- sum(E_TSS_pairnum_core_nocore$standweight)
  whole_avg <- whole_weight / whole_edg_num
  E_TSS_pairnum_core$standweight  <- E_TSS_pairnum_core_nocore$standweight[1:nrow(E_TSS_pairnum_core)]
  E_TSS_pairnum_thre <- E_TSS_pairnum_core[which(E_TSS_pairnum_core$standweight > 0.4),]
  
  uni_E <- unique(E_TSS_pairnum_thre$E)
  len_E <- vector()
  for(i in 1:length(uni_E)){
    index <- which(E_linyu[[uni_E[i]]] > 0)
    len_E <- c(len_E,length(index))
  }
  uni_TSS <- unique(E_TSS_pairnum_thre$TSS)
  len_TSS <- vector()
  for(i in 1:length(uni_TSS)){
    index <- which(TSS_linyu[[uni_TSS[i]]] > 0)
    len_TSS <- c(len_TSS,length(index))
  }
  
  
  hehe <- function(q){
    #对于每一行开始找邻域
    out <- list()
    #cat("q:",q,"\n")
    TSS_point <- i
    m0 <- 1  # edg num
    E_point <- q
    index <- which(E_TSS_pairnum_thre$TSS == i & E_TSS_pairnum_thre$E == q)
    f_weight <- E_TSS_pairnum_core_nocore$standweight[index] # 上一步中也是原始第一个Pair的边的权重
    f_shiying <- f_weight
    TSS_neigh <- which(TSS_linyu[[E_TSS_pairnum_core_nocore$TSS[index]]] != 0) #找出第一个pair的E中点的邻域
    E_neigh <-  which(E_linyu[[E_TSS_pairnum_core_nocore$E[index]]] != 0) #找出第一个pair中TSS中的点的邻域
    
    
    
    if(length(E_neigh) == 0 & length(TSS_neigh) == 0){
      out2 <- list(E=E_point,TSS=TSS_point)
      #out <- c(out,out2)
      return(out2)
      break
    }
    
    max_index_E <- 1
    max_index_TSS <- 1
    while(length(E_neigh) > 0 | length(TSS_neigh) > 0){
     
      if(length(E_neigh) > 0){
        a <- lapply(1:length(E_neigh),function(w){
          index_index <- E_TSS_pairnum_core_nocore[which(E_TSS_pairnum_core_nocore$E == E_neigh[w]),]
          add_edg_index <- which(index_index$TSS %in% TSS_point)
          add_edg_weight2 <- index_index$standweight[add_edg_index]
          add_edg_weight_E <- sum(add_edg_weight2)
          m_new_add_E <-  length(which(add_edg_weight2 > 0))
          b <- c(add_edg_weight_E,m_new_add_E)
        }) 
        a <- unlist(a)
        f_new_E <- a[seq(1,length(a),2)]
        m_new_E <- a[seq(2,length(a),2)]
        f_new_E_weight <- f_new_E + f_weight
        f_new_E <- f_new_E_weight - 0.5 * whole_avg * ((length(E_point) + length(TSS_point) + 1) * (length(E_point) + length(TSS_point)) - 2*(m0 + m_new_E))
        delta_E <- f_new_E - f_shiying
      }
      
      if(length(TSS_neigh) > 0){
        c <- lapply(1:length(TSS_neigh),function(t){
          index_index2 <- E_TSS_pairnum_core_nocore[which(E_TSS_pairnum_core_nocore$TSS == TSS_neigh[t]),]
          add_edg_index <- which(index_index2$E %in% E_point)
          add_edg_weight2 <- index_index2$standweight[add_edg_index]
          add_edg_weight_TSS <- sum(add_edg_weight2)
          m_new_add_TSS <- length(which(add_edg_weight2 > 0))
          d <- c(add_edg_weight_TSS,m_new_add_TSS)
          
        })
        c <- unlist(c)
        f_new_TSS <- c[seq(1,length(c),2)]
        m_new_TSS <- c[seq(2,length(c),2)]
        f_new_TSS_weight <- f_new_TSS + f_weight
        f_new_TSS <- f_new_TSS_weight - 0.5 * whole_avg * ((length(E_point) + length(TSS_point) + 1) * (length(E_point) + length(TSS_point)) - 2*(m0 + m_new_TSS))
        delta_TSS <- f_new_TSS - f_shiying
      }
    
      if(length(delta_E) > 0){
        max_index_E <- max(delta_E)
      }
      else{
        max_index_E <- -10000
      }
      
      if(length(delta_TSS) > 0){
        max_index_TSS <- max(delta_TSS)
      }
      else{
        max_index_TSS <- -10000
      }
      
      if(max_index_E <=0 & max_index_TSS <=0){
        out2 <- list(E=E_point,TSS=TSS_point)
        #out <- c(out,out2)
        #return(out)
        break
      }
      cat("hh1","\n")
      
      if(length(TSS_neigh) == 0){
        max_index_TSS <- -10000
      }
      
      if(length(E_neigh) == 0){
        max_index_E <- -10000
      }
      
      if(max_index_E >= max_index_TSS & max_index_E > 0){
        index_E <- which(delta_E == max_index_E)
        if(length(index_E) > 1){
          index_E <- index_E[1]
        }
        E_point <- c(E_point,E_neigh[index_E])
        E_neigh <- E_neigh[-index_E]
        m0 <- m0 + m_new_E[index_E]
        f_weight <- f_new_E_weight[index_E]
        f_shiying <- f_new_E[index_E]
        delta_E <- delta_E[-index_E]
        #cat("length(E_neigh):",length(E_neigh),"\n")
      }
      
      if(max_index_E < max_index_TSS & max_index_TSS > 0){
        index_TSS <- which(delta_TSS == max_index_TSS)
        if(length(index_TSS) > 1){
          index_TSS <- index_TSS[1]
        }
        TSS_point <- c(TSS_point,TSS_neigh[index_TSS])
        TSS_neigh <- TSS_neigh[-index_TSS]
        m0 <- m0 + m_new_TSS[index_TSS]
        f_weight <- f_new_TSS_weight[index_TSS]
        f_shiying <- f_new_TSS[index_TSS]
        delta_TSS <- delta_TSS[-index_TSS]
        #cat("length(TSS_neigh):",length(TSS_neigh),"\n")
      }
      
      if(length(E_neigh) == 0 & length(TSS_neigh) == 0){
        break
      }
      
    }
    out2 <- list(E=E_point,TSS=TSS_point)
    return(out2)
  }
  
  library(foreach)
  library(doParallel)
  cl <- makeCluster(38)
  registerDoParallel(cl)
  
  TSS  <- unique(E_TSS_pairnum_thre$TSS)
  inner <- foreach(i=TSS[1:length(TSS)], .packages = 'foreach') %dopar% {
    E_indexx <- which(E_TSS_pairnum_thre$TSS == i)
    E <- unique(E_TSS_pairnum_thre$E[E_indexx])
    innerr2 <- foreach(x=E) %dopar% hehe(x)
    innerr <- innerr2
    
  }
  hebing_fina_0.4 <- inner
  
  stopCluster(cl)
  
  return(hebing_fina_0.4)
  
  
}