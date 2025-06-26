#!/bin/bash

# Request resources:
#SBATCH -c 16     # 1 entire node
#SBATCH --time=5:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=16G      # 1 GB RAM
#SBATCH -p shared


module load python
module load ffmpeg
pip install pandas
pip install vtk

rm outframes/*

python make_video.py
