kmer_in_seq <- function(Eshijixulie_num_sum,E_core_nocore,TSSshijixulie_num_sum,TSS_core_nocore){
  
  fun_E <-function(x) {
    index <- as.numeric(vector())
    for(i in 1:length(Eshijixulie_num_sum)){
      index2 <- which(Eshijixulie_num_sum[[i]] %in% E_core_nocore$kmer[x])
      
      if(length(index2) > 0){
        index <- c(index,i)
      }
      
    }
    return(index)
  }
  
  
  fun_TSS <-function(x) {
    index <- as.numeric(vector())
    for(i in 1:length(TSSshijixulie_num_sum)){
      index2 <- which(TSSshijixulie_num_sum[[i]] %in% TSS_core_nocore$kmer[x])
      
      
      if(length(index2) > 0){
        index <- c(index,i)
      }
      
    }
    return(index)
  }
  
  library(foreach)
  library(doParallel)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  E_kmer_in_xuliehao <- foreach(x=1:nrow(E_core_nocore)) %dopar% fun_E(x)
  TSS_kmer_in_xuliehao <- foreach(x=1:nrow(TSS_core_nocore)) %dopar% fun_TSS(x)
  
  
  #别忘了结束并行
  stopCluster(cl)
  
  
  
  
  
  return(list(E_kmer_in_xuliehao,TSS_kmer_in_xuliehao))
  
}