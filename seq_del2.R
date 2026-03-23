seq_del2 <- function(Eshijixuliequan,TSSshijixuliequan){
  k=6
  Eshijixuliequan2 <- lapply(1:length(Eshijixuliequan),function(x){
    line <- unlist(strsplit(Eshijixuliequan[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
  })
  
  Eshijixuliequan <- Eshijixuliequan2
  
  TSSshijixuliequan2 <- lapply(1:length(TSSshijixuliequan),function(x){
    line <- unlist(strsplit(TSSshijixuliequan[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
  })
  
  TSSshijixuliequan <- TSSshijixuliequan2
  
  Eshijixuliequan_length <- lengths(Eshijixuliequan)
  TSSshijixuliequan_length <- lengths(TSSshijixuliequan)
  
  num <- as.numeric(vector())
  Eshijixulie_num_sum <- lapply(1:length(Eshijixuliequan),function(x){
    for(i in 1:(length(Eshijixuliequan[[x]])-k+1)){
      numm <- Eshijixuliequan[[x]][i]*(4**5)+Eshijixuliequan[[x]][i+1]*(4**4)+Eshijixuliequan[[x]][i+2]*(4**3)+Eshijixuliequan[[x]][i+3]*(4**2)+Eshijixuliequan[[x]][i+4]*4+Eshijixuliequan[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  num <- as.numeric(vector())
  TSSshijixulie_num_sum <- lapply(1:length(TSSshijixuliequan),function(x){
    for(i in 1:(length(TSSshijixuliequan[[x]])-k+1)){
      numm <- TSSshijixuliequan[[x]][i]*(4**5)+TSSshijixuliequan[[x]][i+1]*(4**4)+TSSshijixuliequan[[x]][i+2]*(4**3)+TSSshijixuliequan[[x]][i+3]*(4**2)+TSSshijixuliequan[[x]][i+4]*4+TSSshijixuliequan[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  
  
  return(list(Eshijixuliequan,TSSshijixuliequan,Eshijixulie_num_sum,TSSshijixulie_num_sum))
  
}