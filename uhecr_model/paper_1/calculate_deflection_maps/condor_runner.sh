#!/bin/bash

# source ~/.bashrc
# echo $PATH
eval "$(/opt/miniconda3/condabin/conda shell.bash hook)"
conda activate uhecr_env
# python -m site
# python -m pip list
python run_gmf_deflections.py "$@"