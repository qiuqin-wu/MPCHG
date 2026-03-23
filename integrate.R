integrate <- function(hebing_fina_0.4){
  
  a_0.4 <- list()
  b_0.4 <- list()
  for(i in 1:length(hebing_fina_0.4)){
    a_0.4 <- c(a_0.4,hebing_fina_0.4[[i]])
  }
  
  for(i in 1:length(a_0.4)){
    b_0.4 <- c(b_0.4,a_0.4[[i]])
  }
  
  E_union_0.4 <- b_0.4[seq(1,length(b_0.4),2)]
  TSS_union_0.4 <- b_0.4[seq(2,length(b_0.4),2)]
  
  num_E_0.4 <- lengths(E_union_0.4,use.names = F)
  num_TSS_0.4 <- lengths(TSS_union_0.4,use.names = F)
  index_E_0.4 <- which(num_E_0.4  <= 5)
  index_TSS_0.4 <- which(num_TSS_0.4  <= 5)
  index_0.4 <- unique(c(index_E_0.4,index_TSS_0.4))
  if(length(index_0.4) > 0){
    E_union_0.4 <- E_union_0.4[-index_0.4]
    TSS_union_0.4 <- TSS_union_0.4[-index_0.4]
  }
  
  inter_E <- function(x){
    num <- as.numeric(vector())
    pro <- as.numeric(vector())
    for(i in 1:length(E_union_0.4)){
      num2 <- length(intersect(E_union_0.4[[x]],E_union_0.4[[i]]))
      num <- c(num,num2)
      pro2 <- num2/length(E_union_0.4[[x]])
      pro <- c(pro,pro2)
    }
    pro[x] <- 0
    return(pro)
  }
  
  library(foreach)   
  library(doParallel)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  inters_E <- foreach(x=1:length(E_union_0.4)) %dopar% inter_E(x)
  stopCluster(cl)
  
  k_nn <- 50  #k最近邻
  b_E <- vector("list", length(inters_E))
  for(i in 1:length(inters_E)){
    #index <- length(which(inters[[i]] > 0.5))
    b_E[[i]] <- order(inters_E[[i]],decreasing=TRUE)[1:(k_nn-1)]
    
  }
  
  nn_map_E <- as.data.frame(do.call(rbind,b_E)) # 每个pair_E边的邻居，每一列都表示邻居
  nn_map_E <- cbind(nn_map_E, seq_len(nrow(nn_map_E)))
  
  r <- 50  #zui jin lin
  good_choices_E <- seq_len(length(E_union_0.4)) # 生成一个向量1:length(E_union_0.5)
  choice_E <- sample(seq_len(length(good_choices_E)), size = 1, replace = FALSE)
  chosen_E <- good_choices_E[choice_E] # 被抽样的选项
  good_choices_E <- good_choices_E[good_choices_E != good_choices_E[choice_E]] #备选项
  it <- 0
  r2 <- r * 2
  # function for sapply
  get_shared_E <- function(other, this_choice) {
    r2 - length(union(cell_sample_E[other,], this_choice))
  }
  
  while (length(good_choices_E) > 0 & it < 5000) { # slow
    
    it <- it + 1
    choice_E <- sample(seq_len(length(good_choices_E)), size = 1, replace = FALSE)
    new_chosen_E <- c(chosen_E, good_choices_E[choice_E])
    good_choices_E <- good_choices_E[good_choices_E != good_choices_E[choice_E]]
    cell_sample_E <- nn_map_E[new_chosen_E,]
    
    others_E <- seq_len(nrow(cell_sample_E) - 1)
    this_choice_E <- cell_sample_E[nrow(cell_sample_E),]
    shared_E <- sapply(others_E, get_shared_E, this_choice = this_choice_E)
    cat("shared:",max(shared_E),"\n")
    if (max(shared_E) < .9 * r) {
      chosen_E <- new_chosen_E
    }
  }
  
  inter_TSS <- function(x){
    num <- as.numeric(vector())
    pro <- as.numeric(vector())
    for(i in 1:length(TSS_union_0.4)){
      num2 <- length(intersect(TSS_union_0.4[[x]],TSS_union_0.4[[i]]))
      num <- c(num,num2)
      pro2 <- num2/length(TSS_union_0.4[[x]])
      pro <- c(pro,pro2)
    }
    pro[x] <- 0
    return(pro)
  }
  
  library(foreach)   
  library(doParallel)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  inters_TSS <- foreach(x=1:length(TSS_union_0.4)) %dopar% inter_TSS(x)
  stopCluster(cl)
  
  k_nn <- 50  #k最近邻
  b_TSS <- vector("list", length(inters_TSS))
  for(i in 1:length(inters_TSS)){
    #index <- length(which(inters[[i]] > 0.5))
    b_TSS[[i]] <- order(inters_TSS[[i]],decreasing=TRUE)[1:(k_nn-1)]
    
  }
  nn_map_TSS <- as.data.frame(do.call(rbind,b_TSS)) # 每个pair_E边的邻居，每一列都表示邻居
  nn_map_TSS <- cbind(nn_map_TSS, seq_len(nrow(nn_map_TSS)))
  
  r <- 50  #zui jin lin
  good_choices_TSS <- seq_len(length(TSS_union_0.4)) # 生成一个向量1:length(E_union_0.5)
  choice_TSS <- sample(seq_len(length(good_choices_TSS)), size = 1, replace = FALSE)
  chosen_TSS <- good_choices_TSS[choice_TSS] # 被抽样的选项
  good_choices_TSS <- good_choices_TSS[good_choices_TSS != good_choices_TSS[choice_TSS]] #备选项
  it <- 0
  r2 <- r * 2
  # function for sapply
  get_shared_TSS <- function(other, this_choice) {
    r2 - length(union(cell_sample_TSS[other,], this_choice))
  }
  
  while (length(good_choices_TSS) > 0 & it < 5000) { # slow
    
    it <- it + 1
    choice_TSS <- sample(seq_len(length(good_choices_TSS)), size = 1, replace = FALSE)
    new_chosen_TSS <- c(chosen_TSS, good_choices_TSS[choice_TSS])
    good_choices_TSS <- good_choices_TSS[good_choices_TSS != good_choices_TSS[choice_TSS]]
    cell_sample_TSS <- nn_map_TSS[new_chosen_TSS,]
    
    others_TSS <- seq_len(nrow(cell_sample_TSS) - 1)
    this_choice_TSS <- cell_sample_TSS[nrow(cell_sample_TSS),]
    shared_TSS <- sapply(others_TSS, get_shared_TSS, this_choice = this_choice_TSS)
    cat("shared:",max(shared_TSS),"\n")
    if (max(shared_TSS) < .9 * r) {
      chosen_TSS <- new_chosen_TSS
    }
  }
  
  #first 找出KNN 得出的cluster
  
  E_union2_0.4 <- E_union_0.4
  TSS_union2_0.4 <- TSS_union_0.4
  
  nn_map_E2 <- nn_map_E[chosen_E,]
  for(i in 1:nrow(nn_map_E2)){
    for(j in 1:r){
      a <- E_union2_0.4[[nn_map_E2[i,j]]]
      b <- E_union2_0.4[[nn_map_E2[i,r]]]
      bili <- length(which(a %in% b))/length(b)
      if(bili < 0.9){
        nn_map_E2[i,j] <- nn_map_E2[i,r]
      }
    }
  }
  nn_map_E3 <- nn_map_E
  nn_map_E3[chosen_E,] <- nn_map_E2
  
  for(i in 1:length(chosen_E)){
    index <- chosen_E[i]
    index_vector <- unique(as.numeric(nn_map_E3[index,]))
    union <- unique(unlist(E_union2_0.4[index_vector]))
    for(j in 1:length(index_vector)){
      E_union2_0.4[[index_vector[j]]] <- union
    }
  }
  
  # 计算E motif 的个数
  list_E <- vector("list", length(E_union2_0.4))
  index_done_E <- matrix(0,nrow=length(E_union2_0.4),ncol=1)
  
  for(i in 1:length(E_union2_0.4)){
    vec <- vector()
    if(index_done_E[i] ==0){
      for(j in i:length(E_union2_0.4)){
        
        if(index_done_E[j] ==0){ 
          
          index <- identical(E_union2_0.4[[i]],E_union2_0.4[[j]])
          if(index==TRUE){
            vec <- c(vec,j)
            index_done_E[j] <- 1
          }
          
        }
      }
      if(length(vec) > 0){
        list_E[[i]] <- vec
      }
    }
    
  }
  E_motif_num <- length(which(lengths(list_E)>0))
  E_motif_index  <- list_E[which(lengths(list_E)>0)] 
  
  nn_map_TSS2 <- nn_map_TSS[chosen_TSS,]
  for(i in 1:nrow(nn_map_TSS2)){
    for(j in 1:r){
      a <- TSS_union2_0.4[[nn_map_TSS2[i,j]]]
      b <- TSS_union2_0.4[[nn_map_TSS2[i,r]]]
      bili <- length(which(a %in% b))/length(b)
      if(bili < 0.9){
        nn_map_TSS2[i,j] <- nn_map_TSS2[i,r]
      }
    }
  }
  nn_map_TSS3 <- nn_map_TSS
  nn_map_TSS3[chosen_TSS,] <- nn_map_TSS2
  
  for(i in 1:length(chosen_TSS)){
    index <- chosen_TSS[i]
    index_vector <- unique(as.numeric(nn_map_TSS3[index,]))
    union <- unique(unlist(TSS_union2_0.4[index_vector]))
    for(j in 1:length(index_vector)){
      TSS_union2_0.4[[index_vector[j]]] <- union
    }
  }
  
  # 计算TSS motif 的个数
  list_TSS <- vector("list", length(TSS_union2_0.4))
  index_done_TSS <- matrix(0,nrow=length(TSS_union2_0.4),ncol=1)
  
  for(i in 1:length(TSS_union2_0.4)){
    vec <- vector()
    if(index_done_TSS[i] ==0){
      for(j in i:length(TSS_union2_0.4)){
        
        if(index_done_TSS[j] ==0){ 
          
          index <- identical(TSS_union2_0.4[[i]],TSS_union2_0.4[[j]])
          if(index==TRUE){
            vec <- c(vec,j)
            index_done_TSS[j] <- 1
          }
          
        }
      }
      if(length(vec) > 0){
        list_TSS[[i]] <- vec
      }
    }
    
  }
  
  TSS_motif_num <- length(which(lengths(list_TSS)>0))
  TSS_motif_index  <- list_TSS[which(lengths(list_TSS)>0)] 
  
  
  
  
  return(list(E_union_0.4,TSS_union_0.4,inters_E,chosen_E,inters_TSS,chosen_TSS,E_union2_0.4,E_motif_index,TSS_union2_0.4,TSS_motif_index))
  
}