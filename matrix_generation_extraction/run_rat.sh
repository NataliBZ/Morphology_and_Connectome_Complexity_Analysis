#!/bin/bash

#SBATCH -N 1                 # number of nodes
#SBATCH --mem 350G           # memory
#SBATCH --time=20:00:00      # Run for a maximum time of 0 days, 12 hours, 00 mins, 00 secs
#SBATCH --exclusive
#SBATCH --job-name=ConnMatRat
#SBATCH --account=proj159
#SBATCH --partition=prod

# Activate the virtual environment
source venv_jupyterhub/bin/activate 

# go to the path with the python script
cd /gpfs/bbp.cscs.ch/home/barros/ConnectomeUtilities/

# Run the script
python GetConnMatBIG_Rat.py


