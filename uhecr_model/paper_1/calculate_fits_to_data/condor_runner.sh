#!/bin/bash

#eval "$(/opt/miniconda3/condabin/conda shell.bash hook)"

# conda env @ MPP
eval "$(/opt/anaconda3/condabin/conda  shell.bash hook)"
conda activate uhecr_env

python run_data_fits.py "$@"

