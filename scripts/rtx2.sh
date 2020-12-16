#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH -D /home/pbuividovich/ExDiag/
#SBATCH -J trotter
#SBATCH --output=/home/pbuividovich/ExDiag/data/trotter2.output
#SBATCH --error=/home/pbuividovich/ExDiag/data/trotter2.errors
#SBATCH --gres=gpu:0 --cpus-per-task=8 --mem=4G --partition gtx

cd /home/pbuividovich/ExDiag/

echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
pwd

/home/pbuividovich/ExDiag/trotter --h 1.0 --tmax 20.0 --dt 0.025 --L 12
