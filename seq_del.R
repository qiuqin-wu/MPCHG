seq_del <- function(Eshijixuliequchong,Ebeijingxulie,TSSshijixuliequchong,TSSbeijingxulie){
  k=6
  Eshijixuliequchong2 <- lapply(1:length(Eshijixuliequchong),function(x){
    line <- unlist(strsplit(Eshijixuliequchong[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
  })
  
  Eshijixuliequchong <- Eshijixuliequchong2
  EbeijingxuliTSS <- lapply(1:length(Ebeijingxulie),function(x){
    line <- unlist(strsplit(Ebeijingxulie[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
  })
  
  Ebeijingxulie <- EbeijingxuliTSS
  
  
  TSSshijixuliequchong2 <- lapply(1:length(TSSshijixuliequchong),function(x){
    line <- unlist(strsplit(TSSshijixuliequchong[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
  })
  
  TSSshijixuliequchong <- TSSshijixuliequchong2
  
  
  TSSbeijingxuliTSS <- lapply(1:length(TSSbeijingxulie),function(x){
    line <- unlist(strsplit(TSSbeijingxulie[[x]],split=''))
    line[line=="A" | line=="a"] <- 0
    line[line=="T" | line=="t"] <- 1
    line[line=="C" | line=="c"] <- 2
    line[line=="G" | line=="g"] <- 3
    linTSS <- line
    line3 <- as.numeric(linTSS)
    
  })
  
  TSSbeijingxulie <- TSSbeijingxuliTSS
  
  
  Eshijixuliequchong_length <- lengths(Eshijixuliequchong)
  TSSshijixuliequchong_length <- lengths(TSSshijixuliequchong)
  
  num <- as.numeric(vector())
  Eshijixulie_num_sum <- lapply(1:length(Eshijixuliequchong),function(x){
    for(i in 1:(length(Eshijixuliequchong[[x]])-k+1)){
      numm <- Eshijixuliequchong[[x]][i]*(4**5)+Eshijixuliequchong[[x]][i+1]*(4**4)+Eshijixuliequchong[[x]][i+2]*(4**3)+Eshijixuliequchong[[x]][i+3]*(4**2)+Eshijixuliequchong[[x]][i+4]*4+Eshijixuliequchong[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  num <- as.numeric(vector())
  Ebeijingxulie_num_sum <- lapply(1:length(Ebeijingxulie),function(x){
    for(i in 1:(length(Ebeijingxulie[[x]])-k+1)){
      numm <- Ebeijingxulie[[x]][i]*(4**5)+Ebeijingxulie[[x]][i+1]*(4**4)+Ebeijingxulie[[x]][i+2]*(4**3)+Ebeijingxulie[[x]][i+3]*(4**2)+Ebeijingxulie[[x]][i+4]*4+Ebeijingxulie[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  
  num <- as.numeric(vector())
  TSSshijixulie_num_sum <- lapply(1:length(TSSshijixuliequchong),function(x){
    for(i in 1:(length(TSSshijixuliequchong[[x]])-k+1)){
      numm <- TSSshijixuliequchong[[x]][i]*(4**5)+TSSshijixuliequchong[[x]][i+1]*(4**4)+TSSshijixuliequchong[[x]][i+2]*(4**3)+TSSshijixuliequchong[[x]][i+3]*(4**2)+TSSshijixuliequchong[[x]][i+4]*4+TSSshijixuliequchong[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  num <- as.numeric(vector())
  TSSbeijingxulie_num_sum <- lapply(1:length(TSSbeijingxulie),function(x){
    for(i in 1:(length(TSSbeijingxulie[[x]])-k+1)){
      numm <- TSSbeijingxulie[[x]][i]*(4**5)+TSSbeijingxulie[[x]][i+1]*(4**4)+TSSbeijingxulie[[x]][i+2]*(4**3)+TSSbeijingxulie[[x]][i+3]*(4**2)+TSSbeijingxulie[[x]][i+4]*4+TSSbeijingxulie[[x]][i+5]
      num <- c(num,numm) 
    }
    
    num2 <- num
  })
  
  
  return(list(Eshijixuliequchong,Ebeijingxulie,TSSshijixuliequchong,TSSbeijingxuliTSS,Eshijixulie_num_sum,Ebeijingxulie_num_sum,TSSshijixulie_num_sum,
              TSSbeijingxulie_num_sum))
  
}