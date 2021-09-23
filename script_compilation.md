# Processing ONT reads to produce primary assembly
## Step 1. Basecalling using guppy

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
## Step 2. quality check with pycoQC (v2.5.0.3)

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name pycoqc
#SBATCH --mem=1G
#SBATCH --time=00:01:00
#SBATCH --account=uoo02752
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

pycoQC -f ../../sequencing_summary.txt -o pycoQC_output.html
```
Output -----> [PycoQC report](pycoQC_nem.html)

## Step 3. remove labda DNA with Nanolyse (v1.2.0)

Before doing this merge all the fastq files produced by guppy in to one file (i.e. Nem.ont.merged.fastq)
You also need lamda DNA fasta file in the working directory (i.e [dna_cs.fasta](dna_cs.fasta))

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name nanolyse.job
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

cat Nem.ont.merged.fastq | NanoLyse --reference ./dna_cs.fasta | gzip > m.neg_nanopore_filtered.fastq.gz
```
## Step 4. Remove barcodes with Porechop (v0.2.4)

We have removed barcodes during basecalling with guppy using `--trim_barcodes` function, however that can only remove barcodes from the read ends only. So we will use porechop to remove barcodes if any remaining within the reads. We can also use porechop for quality filtering, but we are not doing that here, because of low coverage assembly, we found that quality fitering will also reduce the assembly quality, buscowise.

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name porechop
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
porechop -i ../m.neg_nanopore_filtered.fastq.gz -o m.neg_ont_nanolyse_porechop.fastq.gz --threads 10
```

## Step 5. Quality check using NanoQC (v0.9.4)

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name nanoqc
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

nanoQC ../m.neg_ont_nanolyse_porechop.fastq.gz -o ./
```
Output ----> [NanoQC report](nanoQC.html)

## Step 6. Genome assembly with Flye (v2.7.1)

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=bigmem
#SBATCH --job-name flye.nematode
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Flye/2.7.1-gimkl-2020a-Python-2.7.18

flye --nano-hq m.neg_ont_nanolyse_porechop.fastq.gz -o ./M.neg_flye -t 10 -i 3
```
we will get an assembly file as `assembly.fasta` we will check its stats and quality with [quast](quast.sh)

## Step 7. Flye tends to produce larger assembly then its actual size. We will run purge haplotigs (v1.1.1) to remove haplotigs in our assembly.

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name purgehap.nem
#SBATCH --mem=10G
#SBATCH --time=03:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load SAMtools/1.12-GCC-9.2.0
module load minimap2/2.20-GCC-9.2.0GGq
module load BEDTools/2.29.2-GCC-9.2.0

# step 1 mapping reads back to the genome
minimap2 -t 10 -ax map-ont assembly.fasta m.neg_ont_nanolyse_porechop.fastq.gz \
--secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali

# step 2 producing histogram
export PATH="/nesi/nobackup/uoo02752/.conda/envs/purge_haplotigs_env/bin:$PATH"
purge_haplotigs hist -b aligned.bam -g assembly.fasta -t 10

# step 3 producing coverage stats
purge_haplotigs cov -i aligned.bam.gencov -l 2 -m 15 -h 190 -o coverage_stats.csv

# step 4 purging genome
purge_haplotigs purge -g assembly.fasta -c coverage_stats.csv -b aligned.bam
```

Because of low coverage for our data we didn't get two peaks on our histogram so we used `-l 2` `-m 15` and `-h 190`
We can run each steps separately, after completion we will again run [quast](quast.sh) to check assembly quality.

## Step 8. we will now scaffold genome with LRscaff (v1.1.11) (5 iterations) using the same long reads as input

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name lrscaf.nem.large
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load minimap2/2.20-GCC-9.2.0

# iteration 1
minimap2 -t 10 curated.fasta m.neg_ont_nanolyse_porechop.fastq.gz > ./aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig curated.fasta --alignedFile aln.mm -t mm -p 10 --output ./scaffolds1

# iteration 2
minimap2 -t 10 ./scaffolds1/scaffolds.fasta m.neg_ont_nanolyse_porechop.fastq.gz > ./scaffolds1/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds.fasta --alignedFile ./scaffolds1/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2

# iteration 3
minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds.fasta m.neg_ont_nanolyse_porechop.fastq.gz > ./scaffolds1/scaffolds2/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3

# iteration 4
minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds3/scaffolds.fasta m.neg_ont_nanolyse_porechop.fastq.gz > ./scaffolds1/scaffolds2/scaffolds3/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds3/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/scaffolds3/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3/scaffolds4

# iteration 5
minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds.fasta m.neg_ont_nanolyse_porechop.fastq.gz > ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5
```
After completion of `lrscaff` we will run [quast](quast.sh) for the output of each iterations in our case the fifth iterations looks better so we will use this for further processing.

## Step 9. Gap closing with lrgapcloser 

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name lr-gapNEM
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-gimkl-2017a
export PATH=/nesi/nobackup/uoo02752/bin/LR_Gapcloser/src/:$PATH

sh LR_Gapcloser.sh -i scaffolds.fasta -l m.neg_ont_nanolyse_porechop.fasta -s n -t 10 -r 10
```
This will run lrgapcloser for 10 iterations. We ran [quast](quast.sh) for the output of each iterations and in this case there was not much improvement after 8th so we took that for further processing.

## Step 10. Running Rails(v1.5.1) and Cobbler(v0.6.1) for further scaffolding and gap closing

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name rails.nem
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Perl/5.30.1-GCC-9.2.0
module load minimap2/2.20-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
export PATH="/nesi/nobackup/uoo02752/nematode/bin/RAILS/bin:$PATH"

sh runRAILSminimapSTREAM.sh ite8.assembly.fasta m.neg_ont_nanolyse_porechop.fa 250 0.80 500 2 ont \
/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/SAMtools/1.13-GCC-9.2.0/bin/samtools 10
```
We will run [quast](quast.sh) on the output scaffolded assembly.
We will keep this assembly as primary assembly. We also have 10x chromium linked reads to build upon this assembly.

# Processing 10x linked reads 

## Step 11. Assembly with Supernova assembler
We used different settings and different number of input reads. We got best results with 670milion reads input.

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --partition=hugemem
#SBATCH --job-name nem.SN.L1n2
#SBATCH --mem=350G
#SBATCH --time=80:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Supernova/2.1.1

supernova run --id=nem670Mili --fastqs=/nesi/nobackup/uoo02752/nematode/nematode_linkedreads/10x.raw/all/ --maxreads=670000000 --accept-extreme-coverage
```

We will use this assembly to scaffold the primary assembly produced with ONT reads at step 10

## Step 12. Scaffolding ONT assembly (step 10) with supernova assembly (step 11) with Ragtag (v2.0.1)

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name ragtag.nem
#SBATCH --mem=6G
#SBATCH --time=00:15:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

ragtag.py scaffold nem.670mil.SN.assembly.fasta m.neg_ont_nanolyse_porechop.fa_vs_ite8.assembly.fasta_250_0.80_rails.scaffolds.fa
```

## Step 13. Preparing 10x reads for scaffolding 
We will use 10x linked reads to scaffold assembly from step 12. To do so first we have to change the format of assembly (step 12) and align using longranger.

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name longR_NEM_align
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=/nesi/project/uoo02752/bin/longranger-2.2.2:$PATH

# formating assembly
longranger mkref ragtag.scaffold.fasta

# reads mapping
longranger align --id=Nematode \
--fastqs=/nesi/nobackup/uoo02752/nematode/nematode_linkedreads/10x.raw/all \
--reference=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/refdata-ragtag.scaffold
```
We will use the bam file produced from this step with ArbitR.

## Step 14. Scaffolding with ArbitR (v. 0.2)

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name arbitr.nem
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/ARBitR/src:$PATH"
export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

arbitr.py -i ragtag.scaffold.fasta -o output.arbitr.scaffolds ./Nematode/outs/possorted_bam.bam
```
We ran two iterations of ArbitR (step 13 and 14)

## Step 15. Scaffolding with Arks and Links

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
##SBATCH --qos=debug
#SBATCH --job-name arks.nem
#SBATCH --mem=50G
#SBATCH --time=48:00:00
##SBATCH --time=00:15:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0
module load LINKS/1.8.7-GCC-9.2.0
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Examples:$PATH
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Arks:$PATH


arks-make arks draft=output.2xarbitr.scaffolds reads=barcoded threads=10
