#!/bin/bash
#SBATCH --job-name=render    # Job name
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=natkinsey@ucdavis.edu     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=12gb                     # Job memory request
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=render_%j.log   # Standard output and error log
pwd; hostname; date

# Load R and locally installed libraries
module load R/3.6.3
export R_LIBS=/share/oberbaurlab/nk/R

# Render the Rmd document
cd /share/oberbaurlab/nk/DarkMatter/scripts
Rscript -e "rmarkdown::render('screenCandidates.Rmd')"