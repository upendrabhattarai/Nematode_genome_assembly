# Step 1. Basecalling using guppy

```
#!/bin/bash -e

#SBATCH --job-name=guppy_nem                 
#SBATCH --account=uoo02752              
#SBATCH --time=10:00:00                 
#SBATCH --partition=gpu                 
#SBATCH --gres=gpu:1                    
#SBATCH --mem=6G                                
#SBATCH --ntasks=4                             
#SBATCH --cpus-per-task=1               
#SBATCH --output=%x-%j.out             
#SBATCH --error=%x-%j.err               
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz

module load ont-guppy-gpu/5.0.7
guppy_basecaller -i ../ -s . --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 4 -x auto --recursive --trim_barcodes --disable_qscore_filtering
```
