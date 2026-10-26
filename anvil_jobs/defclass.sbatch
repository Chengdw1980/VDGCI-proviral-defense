#!/bin/bash
#SBATCH --job-name=defclass
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/x-dcheng3/defclass_%j.out
/anvil/projects/x-bio260593/vdgci/env/bin/python ~/defclass.py
echo "== DONE $(date -Is)"
