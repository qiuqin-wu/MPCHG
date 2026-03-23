kmer_num <- function(kmer){
  char <- as.data.frame(matrix(nrow=nrow(kmer),ncol=1))
  char[,1] <- paste(kmer[,1],kmer[,2],kmer[,3],kmer[,4],kmer[,5],kmer[,6],sep="")
  complementary_DNA=mgsub(char,   #原始序列
                          c("A","T","G","C","a","g","t","c","N","n"),  #原始碱基
                          c("T","A","C","G","t","c","a","g","N","n")   #互补碱基
  )
  #使用stri_reverse获取反向序列
  rev_complementary_DNA=stri_reverse(complementary_DNA[,1])
  char_rev_complementary <- as.data.frame(rev_complementary_DNA)
  
  kmer_rev_complementary <- as.data.frame(matrix(nrow=nrow(kmer),ncol=6))
  
  char_rev_complementary <- str_split_fixed(char_rev_complementary$rev_complementary_DNA,"",6)
  kmer_rev_complementary <- as.data.frame(char_rev_complementary)
  
  kmer_bei  <- kmer
  kmer_rev_complementary_bei <- kmer_rev_complementary
  
  kmer[kmer == "A"] <- 0
  kmer[kmer == "T"] <- 1
  kmer[kmer == "C"] <- 2
  kmer[kmer == "G"] <- 3
  kmer_num <- kmer
  kmer_num <- as.data.frame(lapply(kmer_num,as.numeric))
  
  kmer_rev_complementary[kmer_rev_complementary == "A"] <- "0"
  kmer_rev_complementary[kmer_rev_complementary == "T"] <- "1"
  kmer_rev_complementary[kmer_rev_complementary == "C"] <- "2"
  kmer_rev_complementary[kmer_rev_complementary == "G"] <- "3"
  kmer_rev_complementary_num <- kmer_rev_complementary
  kmer_rev_complementary_num <- as.data.frame(lapply(kmer_rev_complementary_num,as.numeric))
  
  
  kmer  <- kmer_bei
  kmer_rev_complementary <- kmer_rev_complementary_bei
  
  
  kmer_num_sum <- lapply(1:nrow(kmer_num),function(x){
    num <- kmer_num[x,1]*(4**5)+kmer_num[x,2]*(4**4)+kmer_num[x,3]*(4**3)+kmer_num[x,4]*(4**2)+kmer_num[x,5]*(4**1)+kmer_num[x,6]
  })
  
  #compute every kmer sijinzhi baishi
  kmer_num_sum <- as.numeric(unlist(kmer_num_sum))
  kmer_num_sum <- as.data.frame(kmer_num_sum)
  
  
  kmer_rev_complementary_num_sum <- lapply(1:nrow(kmer_rev_complementary_num),function(x){
    num <- kmer_rev_complementary_num[x,1]*(4**5)+kmer_rev_complementary_num[x,2]*(4**4)+kmer_rev_complementary_num[x,3]*(4**3)+kmer_rev_complementary_num[x,4]*(4**2)+kmer_rev_complementary_num[x,5]*(4**1)+kmer_rev_complementary_num[x,6]
  })
  
  #compute every fanxianghubu kmer sijinzhi baishi
  kmer_rev_complementary_num_sum <- as.numeric(unlist(kmer_rev_complementary_num_sum))
  kmer_rev_complementary_num_sum <- as.data.frame(kmer_rev_complementary_num_sum)
  
  
  return(list(
    kmer_num = kmer_num,
    kmer_rev_complementary_num = kmer_rev_complementary_num,
    kmer_rev_complementary = kmer_rev_complementary,
    kmer_num_sum = kmer_num_sum,
    kmer_rev_complementary_num_sum = kmer_rev_complementary_num_sum
  ))
  
}

