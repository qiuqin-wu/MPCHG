create_kmer <- function(k = 6, output_folder_path) {
  
  mer <- c('A','T','C','G')
  k=6
  kmer <- as.data.frame(matrix(nrow=1,ncol=6))
  for(i in 1:4){
    for(j in 1:4){
      for(k in 1:4){
        for(p in 1:4){
          for(q in 1:4){
            for(l in 1:4){
              kmer2 <- as.data.frame(matrix(nrow=1,ncol=6))
              kmer2[1,1] <- i
              kmer2[1,2] <- j
              kmer2[1,3] <- k
              kmer2[1,4] <- p
              kmer2[1,5] <- q
              kmer2[1,6] <- l
              kmer <- rbind(kmer,kmer2)
            }
          }
        }
      }
    }
  } 
  kmer <- kmer[-1,]
  
  kmer3 <- as.data.frame(lapply(kmer,as.character))
  
  index1 <- which(kmer3[,1]=="1" & kmer3[,2]=="1" & kmer3[,3]=="1" & kmer3[,4]=="1" & kmer3[,5]=="1" & kmer3[,6]=="1")
  index2 <- which(kmer3[,1]=="2" & kmer3[,2]=="2" & kmer3[,3]=="2" & kmer3[,4]=="2" & kmer3[,5]=="2" & kmer3[,6]=="2")
  index3 <- which(kmer3[,1]=="3" & kmer3[,2]=="3" & kmer3[,3]=="3" & kmer3[,4]=="3" & kmer3[,5]=="3" & kmer3[,6]=="3")
  index4 <- which(kmer3[,1]=="4" & kmer3[,2]=="4" & kmer3[,3]=="4" & kmer3[,4]=="4" & kmer3[,5]=="4" & kmer3[,6]=="4")
  
  kmer3 <- kmer3[-c(index1,index2,index3,index4),]
  kmer3[kmer3 == "1"] <- "A"
  kmer3[kmer3 == "2"] <- "T"
  kmer3[kmer3 == "3"] <- "C"
  kmer3[kmer3 == "4"] <- "G"
  output_file_path_name <- paste(output_folder_path,"kmer.RDS",sep="/")
  saveRDS(kmer3,output_file_path_name)
  
  return(kmer3)
}