#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH -D /home/pbuividovich/ExDiag/
#SBATCH -J otoc4
#SBATCH --output=/home/pbuividovich/ExDiag/data/otoc4.output
#SBATCH --error=/home/pbuividovich/ExDiag/data/otoc4.errors
#SBATCH --gres=gpu:0 --cpus-per-task=8 --mem=4G --partition gtx --qos=long

cd /home/pbuividovich/ExDiag/

echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
pwd

/home/pbuividovich/ExDiag/otoc --h 6.0 --tmax 0.5 --beta 20.0 --dt 0.02 --L 14