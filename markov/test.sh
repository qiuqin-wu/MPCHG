#!/usr/bin/bash
#SBATCH --job-name  markov
#SBATCH --account PCON0022
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=350GB
#SBATCH --gpus-per-node=1


module load gcc-compatibility/8.4.0
cd /fs/ess/PCON0022/wangqi/WQQ/Markov/markov2

for file in ./*.txt
do
./wqqmain  $file
done
