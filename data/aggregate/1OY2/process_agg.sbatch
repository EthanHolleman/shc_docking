#!/bin/bash
#SBATCH -J PROCESS_AGG  # job name goes here
#SBATCH --partition=med
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etholleman@ucdavis.edu
#SBATCH --time=4:00:00
#SBATCH --mem 36000

module load R/4.0.2
Rscript /home/ethollem/jobs/random_docking_trka/aggregate/1OY2/process_agg_trka.r