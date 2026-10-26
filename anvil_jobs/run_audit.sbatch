#!/bin/bash
#SBATCH --job-name=defaudit
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/x-dcheng3/defaudit_%j.out
PY=/anvil/projects/x-bio260593/vdgci/env/bin/python
$PY ~/defaudit.py && $PY ~/defrobust.py
