library(mgsub)   #obtain reverse complementary sequence
library(stringi) 
library(stringr)


# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  message("No or insufficient command-line parameters detected. Using default parameters.")
  
  args <- c(
    "H:/motifpair_taishi/dataset_6tao_TSS-E/GM12878/data",   # 默认输入文件夹
    "Eshijixuliequchong.txt",
    "Ebeijingxulie.txt",
    "TSSshijixuliequchong.txt",
    "TSSbeijingxulie.txt",
    "Eshijixuliequan.txt",
    "TSSshijixuliequan.txt"
  )
}

# 确保提供了输入文件夹路径和至少一个文件名参数
if (length(args) < 7) {
  stop("Error: Please provide both input folder path and at least one file name as arguments.")
}

# 输入文件夹路径
input_folder_path <- args[1]

# 创建输出文件夹路径（在输入文件夹路径基础上新建一个文件夹）
output_folder_path <- file.path(input_folder_path, "our_results")
dir.create(output_folder_path, showWarnings = FALSE)


# 构建完整的文件路径
file_path <- file.path(input_folder_path, args[2])
# 读取文件内容
Eshijixuliequchong <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  Eshijixuliequchong <- c(Eshijixuliequchong,line)
  n=n+1
}
close(con)

file_path <- file.path(input_folder_path, args[3])
Ebeijingxulie <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  Ebeijingxulie <- c(Ebeijingxulie,line)
  n=n+1
}
close(con)

file_path <- file.path(input_folder_path, args[4]) 
TSSshijixuliequchong <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  TSSshijixuliequchong <- c(TSSshijixuliequchong,line)
  n=n+1
}
close(con)

file_path <- file.path(input_folder_path, args[5]) 
TSSbeijingxulie <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  TSSbeijingxulie <- c(TSSbeijingxulie,line)
  n=n+1
}
close(con)

source("create_kmer.R")
cat("create_kmer","\n")

kmer <- create_kmer(k = 6,output_folder_path)

source("kmer_num.R")
result <- kmer_num(kmer)

kmer_num <- result[[1]]
kmer_rev_complementary_num <- result[[2]]
kmer_rev_complementary <- result[[3]]
kmer_num_sum <- result[[4]]
kmer_rev_complementary_num_sum <- result[[5]]

source("seq_del.R")

cat("seq_del","\n")

result2 <- seq_del(Eshijixuliequchong,Ebeijingxulie,TSSshijixuliequchong,TSSbeijingxulie)


Eshijixuliequchong <- result2[[1]]
Ebeijingxulie <- result2[[2]]
TSSshijixuliequchong <- result2[[3]]
TSSbeijingxuliTSS <- result2[[4]]
Eshijixulie_num_sum <- result2[[5]]
Ebeijingxulie_num_sum <- result2[[6]]
TSSshijixulie_num_sum <- result2[[7]]
TSSbeijingxulie_num_sum <- result2[[8]]

source("kmer_fre.R")
cat("kmer_fre","\n")


result3 <- kmer_fre(Eshijixulie_num_sum,kmer_num_sum,kmer_rev_complementary_num_sum,kmer_rev_freq,Ebeijingxulie_num_sum,
                   TSSshijixulie_num_sum,TSSbeijingxulie_num_sum)




Eshiji_kmer_freq <- result3[[1]]
Ebeijing_kmer_freq <- result3[[2]]
TSSshiji_kmer_freq <- result3[[3]]
TSSbeijing_kmer_freq <- result3[[4]]

source("delete_reverse_complement.R")
cat("delete_reverse_complement","\n")


result4 <- delete_reverse_complement(TSSbeijing_kmer_freq,Eshiji_kmer_freq,TSSshiji_kmer_freq,Ebeijing_kmer_freq,kmer)

Eshiji_kmer_freq_remain <- result4[[1]]
Ebeijing_kmer_freq_remain <- result4[[2]]
TSSshiji_kmer_freq_remain <- result4[[3]]
TSSbeijing_kmer_freq_remain <- result4[[4]]
p_E_whole <- result4[[5]]
p_TSS_whole <- result4[[6]]
E_core <- result4[[7]]
E_nocore <- result4[[8]]
E_core_nocore <- result4[[9]]
TSS_core <- result4[[10]]
TSS_nocore <- result4[[11]]
TSS_core_nocore <- result4[[12]]


file_path <- file.path(input_folder_path, args[6]) 
Eshijixuliequan <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  Eshijixuliequan <- c(Eshijixuliequan,line)
  n=n+1
}
close(con)


file_path <- file.path(input_folder_path, args[7]) 
TSSshijixuliequan <- list()
con <- file(file_path,open="r")
n=1
while ( TRUE ) {
  line = readLines(con, n = 1)
  
  if ( length(line) == 0 ) {
    break
  }
  TSSshijixuliequan <- c(TSSshijixuliequan,line)
  n=n+1
}
close(con)

source("seq_del2.R")
cat("seq_del2","\n")


result5 <- seq_del2(Eshijixuliequan,TSSshijixuliequan)
Eshijixuliequan <- result5[[1]]
TSSshijixuliequan <- result5[[2]]
Eshijixulie_num_sum <- result5[[3]]
TSSshijixulie_num_sum <- result5[[4]]


source("kmer_in_seq.R")
cat("kmer_in_seq","\n")

result6 <- kmer_in_seq(Eshijixulie_num_sum,E_core_nocore,TSSshijixulie_num_sum,TSS_core_nocore)
E_kmer_in_xuliehao <- result6[[1]]
TSS_kmer_in_xuliehao <- result6[[2]]

source("neighbor.R")
cat("neighbor","\n")

result7 <- neighbor(E_core,E_nocore,TSS_core,TSS_nocore,Eshijixulie_num_sum,TSSshijixulie_num_sum,E_core_nocore,kmer_num_sum,kmer,kmer_num,kmer_rev_complementary_num_sum,
                    kmer_rev_complementary, kmer_rev_complementary_num,TSS_core_nocore,E_kmer_in_xuliehao,TSS_kmer_in_xuliehao)
E_core_nocore_kmer <- result7[[1]]
E_core_nocore_kmer_num <- result7[[2]]
E_core_nocore_kmer_rev_complementary <- result7[[3]]
E_core_nocore_kmer_rev_complementary_num <- result7[[4]]
TSS_core_nocore_kmer <- result7[[5]]
TSS_core_nocore_kmer_num <- result7[[6]]
TSS_core_nocore_kmer_rev_complementary <- result7[[7]]
TSS_core_nocore_kmer_rev_complementary_num <- result7[[8]]
E_linyu <- result7[[9]]
TSS_linyu <- result7[[10]]
weight_new <- result7[[11]]
E_TSS_pairnum <- result7[[12]]
E_TSS_pairnum_core <- result7[[13]]
E_TSS_pairnum_nocore <- result7[[14]]
E_TSS_pairnum_core_nocore <- result7[[15]]


source("dense_subgraph.R")
cat("dense_subgraph","\n")

result8 <- dense_subgraph(E_TSS_pairnum_core_nocore,E_linyu,TSS_linyu)
hebing_fina_0.4 <- result8

source("integrate.R")
result9 <- integrate(hebing_fina_0.4)
E_union_0.4 <- result9[[1]]
TSS_union_0.4 <-  result9[[2]]
inters_E <- result9[[3]]
chosen_E <-  result9[[4]]
inters_TSS <- result9[[5]]
chosen_TSS <-  result9[[6]]
E_union2_0.4 <- result9[[7]]
E_motif_index <-  result9[[8]]
TSS_union2_0.4 <- result9[[9]]
TSS_motif_index <-  result9[[10]]

output_file_path_name <- paste(output_folder_path,"E_motif_index.RDS",sep="/")
saveRDS(E_motif_index,output_file_path_name)
output_file_path_name <- paste(output_folder_path,"TSS_motif_index.RDS",sep="/")
saveRDS(TSS_motif_index,output_file_path_name)


source("integrate_matrix.R")
cat("integrate_matrix","\n")

result10 <- integrate_matrix(E_union2_0.4,E_core_nocore,TSS_union2_0.4,TSS_core_nocore,E_linyu,E_core_nocore_kmer,TSS_core_nocore_kmer)

E_union3_0.4 <- result10[[1]]
TSS_union3_0.4 <- result10[[2]]
E_union2_0.4_kmer_sum <- result10[[3]]
TSS_union2_0.4_kmer_sum <- result10[[4]]
PWM_E <- result10[[5]]
PWM_TSS <- result10[[6]]

output_file_path_name <- paste(output_folder_path,"PWM_E.RDS",sep="/")
saveRDS(PWM_E,output_file_path_name)
output_file_path_name <- paste(output_folder_path,"PWM_TSS.RDS",sep="/")
saveRDS(PWM_TSS,output_file_path_name)




source("create_fimo.R")

result11 <- create_fimo(Eshijixuliequchong,TSSshijixuliequchong,PWM_E,PWM_TSS,E_motif_index,TSS_motif_index,output_folder_path)

Motif_E_FIMO <- result11[[1]]
Motif_TSS_FIMO <- result11[[2]]

file_path_E <- file.path(input_folder_path, args[8])
file_path_TSS <- file.path(input_folder_path, args[9])

fimo_path_E <- file.path(output_folder_path, "fimo_scan_E")
# 如果目录不存在，先创建
dir.create(fimo_path_E, recursive = TRUE, showWarnings = FALSE)

# 生成命令
cmd1 <- sprintf(
  "fimo --oc %s --parse-genomic-coord %s %s",
  fimo_path_E,
  Motif_E_FIMO,
  file_path_E
)
system(cmd1)

fimo_path_TSS <- file.path(output_folder_path, "fimo_scan_TSS")
dir.create(fimo_path_E, recursive = TRUE, showWarnings = FALSE)

# 生成命令
cmd2 <- sprintf(
  "fimo --oc %s --parse-genomic-coord %s %s",
  fimo_path_TSS,
  Motif_TSS_FIMO,
  file_path_TSS
)
# 执行
system(cmd2)

peaks_E_path <- file.path(input_folder_path, args[10])
peaks_E <- read.table(peaks_E_path,header=F,sep="\t")

peaks_TSS_path <- file.path(input_folder_path, args[11])
peaks_TSS <- read.table(peaks_TSS_path,header=F,sep="\t")

fimo_E_path1 <- paste(output_folder_path,"fimo_scan_E",sep="/")
fimo_E_path <- paste(fimo_E_path1,"fimo.tsv",sep="/")

fimo_TSS_path1 <- paste(output_folder_path,"fimo_scan_TSS",sep="/")
fimo_TSS_path <- paste(fimo_TSS_path1,"fimo.tsv",sep="/")

fimo_E <- read.table(fimo_E_path,header=T,sep="\t")
fimo_TSS <- read.table(fimo_TSS_path,header=T,sep="\t")

source("pair_distance.R")

result12 <- pair_distance(PWM_E,PWM_TSS,E_motif_index,TSS_motif_index,fimo_E,peaks_E,fimo_TSS,peaks_TSS)
PWM_E_zuizhong <- result12[[1]]
PWM_TSS_zuizhong <- result12[[2]]

output_file_path_name <- paste(output_folder_path,"PWM_E_zuizhong.RDS",sep="/")
saveRDS(PWM_E_zuizhong,output_file_path_name)
output_file_path_name <- paste(output_folder_path,"PWM_TSS_zuizhong.RDS",sep="/")
saveRDS(PWM_TSS_zuizhong,output_file_path_name)
