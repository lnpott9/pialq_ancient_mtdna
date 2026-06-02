# 04_contam

*$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/*

I copied all of the potentially contaminated files to a new directory

```
% la -1
f_19dsV2_sorted.md.bam
f_19dsV2_sorted.md.bam.bai
f_56mt_sorted.md.bam
f_56mt_sorted.md.bam.bai
f_78mt_sorted.md.bam
f_78mt_sorted.md.bam.bai
uf_19dsV2_sorted.md.bam
uf_19dsV2_sorted.md.bam.bai
uf_56mt_sorted.md.bam
uf_56mt_sorted.md.bam.bai
uf_78mt_sorted.md.bam
uf_78mt_sorted.md.bam.bai
```

## 01. Schmutzi with prediction of contaminant for 19dsV2, 56mt, and 78mt

*$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/01_wpred*

The first step is to rerun schmutzi with prediction of the contaminant to see if that changes the contamination estimates.

**dsdna (19dsV2 and 78mt)**

```
#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-3

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/dsdna
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
database=/home/mnievesc/shared/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="/users/9/pott0195/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "f_19dsV2_sorted.md.bam"
    "uf_19dsV2_sorted.md.bam"
    "uf_78mt_sorted.md.bam"
    "f_78mt_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
/home/mnievesc/shared/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library double --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi
/home/mnievesc/shared/programs/schmutzi/src/schmutzi.pl --uselength --ref "$ref" --out "${filename}_schmutzi_wpred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_wpred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"
```



**ssdna (56mt)**

```
#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-1

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/ssdna
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
database=/home/mnievesc/shared/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="/users/9/pott0195/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "f_56mt_sorted.md.bam"
    "uf_56mt_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
/home/mnievesc/shared/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library single --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi
/home/mnievesc/shared/programs/schmutzi/src/schmutzi.pl --uselength --ref "$ref" --out "${filename}_schmutzi_wpred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_wpred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"
```

I made a summary file using the following code:

```
text="wpred.txt"

echo -e "Sample\tAverage\tLower\tUpper" > "$text"
for file in *sdna/*_schmutzi/*wpred_final.cont.est; do
    sample=$(basename "$file" _schmutzi_wpred_final.cont.est)
    awk -v sample="$sample" '{print sample "\t" $0}' "$file" >> "$text"
done
```

The summary of this run is as follows:

```
Sample  Average Lower   Upper
f_19dsV2        0.01    0       0.02
f_78mt  0.99    0.98    0.99
uf_19dsV2       0.01    0       0.02
uf_78mt 0.01    0       0.02
f_56mt  0.99    0.98    0.99
uf_56mt 0.96    0.95    0.97
```

19dsV2 has a contamination of 1% when predicting with contaminants, meaning I can use the output fasta from the wpred run for this library. However, 56mt and 78mt will need to be decontaminated using PMDtools and run through schmutzi again to identify when they hit 5% contamination or below. 

## 02. Decontamination with PMDtools

*$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd*

Input files:

```
% ls -1
f_56mt_sorted.md.bam
f_56mt_sorted.md.bam.bai
f_78mt_sorted.md.bam
f_78mt_sorted.md.bam.bai
uf_56mt_sorted.md.bam
uf_56mt_sorted.md.bam.bai
uf_78mt_sorted.md.bam
uf_78mt_sorted.md.bam.bai
```

The [Github for PMDtools is here](https://github.com/pontussk/PMDtools) and I cloned the repository using `git clone` in `$SHARED/programs/PMDtools`. The path to the executable is `/home/mnievesc/shared/programs/PMDtools/pmdtools.0.60.py`. 

There are 2 approaches in PMDtools to remove contaminated sequences:

1. separate ancient DNA molecules from others based on PMD-scores of 3 or above (post mortem damage score that PMDtools calculates)
2. separate ancient DNA molecules from others based on having a C->T mismatch in the 5' and/or 3'

I will be using option 1. The specific PMDtools script I need is `pmdtools.0.60.py`, which is unfortunately written in Python 2. 

Usage for PMDtools is:

```
samtools view -h mybam.bam | python pmdtools.0.60.py --threshold 3 --header | samtools view -Sb - > mybam.pmds3filter.bam`
```



The commands I used to decontaminate 56mt and 78mt were:

```
# start interactive session
srun --mem=8gb --time=60:00 --pty bash

# load modules
module load python2/2.7.12_anaconda4.2
module load samtools/1.21

# set paths
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd
pmdtools=/home/mnievesc/shared/programs/PMDtools/pmdtools.0.60.py

# make array of intput files
files=("f_78mt_sorted.md.bam" "uf_78mt_sorted.md.bam" "f_56mt_sorted.md.bam" "uf_56mt_sorted.md.bam")

# decontaminate at different thresholds
cd "$working"

for bam in "${files[@]}"; do
	filename=$(basename "$bam" _sorted.md.bam)
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 2 --header | samtools view -Sb - > "${filename}_pmd2_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 2.5 --header | samtools view -Sb - > "${filename}_pmd2.5_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 3 --header | samtools view -Sb - > "${filename}_pmd3_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 3.5 --header | samtools view -Sb - > "${filename}_pmd3.5_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 4 --header | samtools view -Sb - > "${filename}_pmd4_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 4.5 --header | samtools view -Sb - > "${filename}_pmd4.5_sorted.md.bam"
done
```

I wanted to see how this impacted coverage, so I also ran the following:

```
for bam in *bam; do
	filename=$(basename "$bam" .bam)
	echo "$filename"
	samtools coverage "$bam" | cut -f6-7
	echo "-----------------"
done

# output
f_56mt_pmd2.5_sorted.md
coverage        meandepth
96.3909 18.4706
-----------------
f_56mt_pmd2_sorted.md
coverage        meandepth
96.409  20.8924
-----------------
f_56mt_pmd3.5_sorted.md
coverage        meandepth
94.7311 13.9645
-----------------
f_56mt_pmd3_sorted.md
coverage        meandepth
96.1494 15.623
-----------------
f_56mt_pmd4.5_sorted.md
coverage        meandepth
93.1619 9.86746
-----------------
f_56mt_pmd4_sorted.md
coverage        meandepth
93.343  12.2027
-----------------
f_56mt_sorted.md
coverage        meandepth
99.994  200.814
-----------------
f_78mt_pmd2.5_sorted.md
coverage        meandepth
74.5368 3.22868
-----------------
f_78mt_pmd2_sorted.md
coverage        meandepth
75.6775 3.5558
-----------------
f_78mt_pmd3.5_sorted.md
coverage        meandepth
69.9258 2.61199
-----------------
f_78mt_pmd3_sorted.md
coverage        meandepth
73.5772 2.92866
-----------------
f_78mt_pmd4.5_sorted.md
coverage        meandepth
55.3202 1.39791
-----------------
f_78mt_pmd4_sorted.md
coverage        meandepth
64.4457 2.14551
-----------------
f_78mt_sorted.md
coverage        meandepth
98.8895 26.6705
-----------------
uf_56mt_pmd2.5_sorted.md
coverage        meandepth
99.0947 21.7897
-----------------
uf_56mt_pmd2_sorted.md
coverage        meandepth
99.0947 24.5252
-----------------
uf_56mt_pmd3.5_sorted.md
coverage        meandepth
97.8937 16.241
-----------------
uf_56mt_pmd3_sorted.md
coverage        meandepth
98.7265 18.2539
-----------------
uf_56mt_pmd4.5_sorted.md
coverage        meandepth
97.2418 11.5101
-----------------
uf_56mt_pmd4_sorted.md
coverage        meandepth
97.4531 14.279
-----------------
uf_56mt_sorted.md
coverage        meandepth
100     223.792
-----------------
uf_78mt_pmd2.5_sorted.md
coverage        meandepth
87.6818 4.46086
-----------------
uf_78mt_pmd2_sorted.md
coverage        meandepth
88.2854 4.86384
-----------------
uf_78mt_pmd3.5_sorted.md
coverage        meandepth
84.9176 3.61199
-----------------
uf_78mt_pmd3_sorted.md
coverage        meandepth
86.7101 4.0446
-----------------
uf_78mt_pmd4.5_sorted.md
coverage        meandepth
71.5493 1.93186
-----------------
uf_78mt_pmd4_sorted.md
coverage        meandepth
81.4231 3.02444
-----------------
uf_78mt_sorted.md
coverage        meandepth
99.831  33.7385
-----------------
```

78mt does not seem like it has enough endogenous coverage to be useful after this step, but 56mt looks like a better option. 



## 03. Schmutzi (split by  library prep) with no contaminant prediction for 56mt and 78mt

*$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/03_npred*

I then repeated the earlier step of running schmutzi with no contaminant prediction on 56mt and 78mt. They are already sorted and md tagged.

**Running *schmutzi* with no contaminant prediction (--notusepredC)**

**ssdna:**

```
#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-11

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/03_npred
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
database=/home/mnievesc/shared/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="/users/9/pott0195/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "f_56mt_pmd2.5_sorted.md.bam"
    "f_56mt_pmd2_sorted.md.bam"
    "f_56mt_pmd3.5_sorted.md.bam"
    "f_56mt_pmd3_sorted.md.bam"
    "f_56mt_pmd4.5_sorted.md.bam"
    "f_56mt_pmd4_sorted.md.bam"
    "uf_56mt_pmd2.5_sorted.md.bam"
    "uf_56mt_pmd2_sorted.md.bam"
    "uf_56mt_pmd3.5_sorted.md.bam"
    "uf_56mt_pmd3_sorted.md.bam"
    "uf_56mt_pmd4.5_sorted.md.bam"
    "uf_56mt_pmd4_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
/home/mnievesc/shared/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library single --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi
/home/mnievesc/shared/programs/schmutzi/src/schmutzi.pl --notusepredC --uselength --ref "$ref" --out "${filename}_schmutzi_npred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_npred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"
```

**dsdna:**

```
#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-11

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/03_npred
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
database=/home/mnievesc/shared/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="/users/9/pott0195/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "f_78mt_pmd2.5_sorted.md.bam"
    "f_78mt_pmd2_sorted.md.bam"
    "f_78mt_pmd3.5_sorted.md.bam"
    "f_78mt_pmd3_sorted.md.bam"
    "f_78mt_pmd4.5_sorted.md.bam"
    "f_78mt_pmd4_sorted.md.bam"
    "uf_78mt_pmd2.5_sorted.md.bam"
    "uf_78mt_pmd2_sorted.md.bam"
    "uf_78mt_pmd3.5_sorted.md.bam"
    "uf_78mt_pmd3_sorted.md.bam"
    "uf_78mt_pmd4.5_sorted.md.bam"
    "uf_78mt_pmd4_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
/home/mnievesc/shared/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library double --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi
/home/mnievesc/shared/programs/schmutzi/src/schmutzi.pl --notusepredC --uselength --ref "$ref" --out "${filename}_schmutzi_npred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_npred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"
```

I created a summary text file with this:

```
text="pmd_npred.txt"

echo -e "Sample\tAverage\tLower\tUpper" > "$text"
for file in *_schmutzi/*npred_final.cont.est; do
    sample=$(basename "$file" _schmutzi_npred_final.cont.est)
    awk -v sample="$sample" '{print sample "\t" $0}' "$file" >> "$text"
done
```

Based on these results, I will use f_56mt_pmd3 and uf_56mt_pmd3 in further analyses, and not use 78mt at all (except for variant calling). A PMD score of 3 is a high-confidence threshold that a fragment is of ancient origin, and I do not want to lose more fragments than I have too by using a higher PMD score. 



78mt still has >5% contamination, and at at higher PMD score cutoffs the coverage will be way too low to create a quality consensus sequence. 

```
Sample  Average Lower   Upper
f_56mt_pmd2.5   0.04    0.03    0.05
f_56mt_pmd2     0.04    0.03    0.05
f_56mt_pmd3.5   0.04    0.03    0.05
f_56mt_pmd3     0.04    0.03    0.05
f_56mt_pmd4.5   0.03    0.02    0.04
f_56mt_pmd4     0.03    0.02    0.04
f_78mt_pmd2.5   0.11    0.09    0.13
f_78mt_pmd2     0.11    0.09    0.13
f_78mt_pmd3.5   0.09    0.07    0.11
f_78mt_pmd3     0.1     0.08    0.12
f_78mt_pmd4.5   0.09    0.07    0.11
f_78mt_pmd4     0.1     0.08    0.12
uf_56mt_pmd2.5  0.04    0.03    0.05
uf_56mt_pmd2    0.06    0.05    0.07
uf_56mt_pmd3.5  0.03    0.02    0.04
uf_56mt_pmd3    0.04    0.03    0.05
uf_56mt_pmd4.5  0.03    0.02    0.04
uf_56mt_pmd4    0.03    0.02    0.04
uf_78mt_pmd2.5  0.1     0.09    0.11
uf_78mt_pmd2    0.1     0.09    0.11
uf_78mt_pmd3.5  0.09    0.08    0.1
uf_78mt_pmd3    0.09    0.07    0.11
uf_78mt_pmd4.5  0.08    0.06    0.1
uf_78mt_pmd4    0.08    0.07    0.09
```



## 04. Schmutzi with contaminant prediction for 78mt

I will run 78mt through schmutzi with contaminant prediction to see if that helps the contamination levels like it did for 19dsV2. 

```
#!/bin/bash -l
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-11

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/04_wpred
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
database=/home/mnievesc/shared/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="/users/9/pott0195/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "f_78mt_pmd2.5_sorted.md.bam"
    "f_78mt_pmd2_sorted.md.bam"
    "f_78mt_pmd3.5_sorted.md.bam"
    "f_78mt_pmd3_sorted.md.bam"
    "f_78mt_pmd4.5_sorted.md.bam"
    "f_78mt_pmd4_sorted.md.bam"
    "uf_78mt_pmd2.5_sorted.md.bam"
    "uf_78mt_pmd2_sorted.md.bam"
    "uf_78mt_pmd3.5_sorted.md.bam"
    "uf_78mt_pmd3_sorted.md.bam"
    "uf_78mt_pmd4.5_sorted.md.bam"
    "uf_78mt_pmd4_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
/home/mnievesc/shared/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library double --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi
/home/mnievesc/shared/programs/schmutzi/src/schmutzi.pl --uselength --ref "$ref" --out "${filename}_schmutzi_wpred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_wpred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"
```



Summary file creation:

```
text="pmd_wpred.txt"

echo -e "Sample\tAverage\tLower\tUpper" > "$text"
for file in *_schmutzi/*wpred_final.cont.est; do
    sample=$(basename "$file" _schmutzi_wpred_final.cont.est)
    awk -v sample="$sample" '{print sample "\t" $0}' "$file" >> "$text"
done
```

The results are here:

```
Sample  Average Lower   Upper
f_78mt_pmd2.5   0.01    0       0.02
f_78mt_pmd2     0.01    0       0.02
f_78mt_pmd3     0.01    0       0.02
f_78mt_pmd4.5   0.13    0.11    0.15
f_78mt_pmd4     0.07    0.05    0.09
uf_78mt_pmd2.5  0.01    0       0.02
uf_78mt_pmd2    0.01    0       0.02
uf_78mt_pmd3.5  0.01    0       0.02
uf_78mt_pmd3    0.01    0       0.02
uf_78mt_pmd4.5  0.07    0.05    0.09
uf_78mt_pmd4    0.01    0       0.02
```

I will use the wpred version of 78mt_pmd3 in future analyses. 



## vcf production

```
#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# dependencies
export PATH=$SHARED/programs/nextflow:$PATH
module load singularity/current
module load java/openjdk-17.0.2

# set working directory
working=$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd
cd "$working"

# pull most recent version of eager
NXF_VER='22.10.6' nextflow pull nf-core/eager

# run EAGER script
NXF_VER='22.10.6' nextflow run nf-core/eager \
-r 2.5.0 \
-profile singularity \
--input 'f*pmd3_sorted.md.bam' \
--single_end \
--bam \
--fasta $SHARED/ref_seqs/rCRS.fasta \
--fasta_index $SHARED/ref_seqs/rCRS.fasta.fai \
--outdir './05_vcf' \
--skip_fastqc \
--skip_adapterremoval \
--skip_preseq \
--skip_deduplication \
--skip_qualimap \
--damage_calculation_tool 'mapdamage' \
--run_mapdamage_rescaling \
--run_genotyping \
--genotyping_tool 'ug' \
--genotyping_source 'rescaled' \
--gatk_ploidy 1 \
--gatk_ug_out_mode 'EMIT_ALL_CONFIDENT_SITES' \
--gatk_ug_genotype_model 'BOTH' \
--gatk_ug_keep_realign_bam

# clean 
# nextflow clean -f -k
```

cleaned vcfs and bam files are in $SHARED/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd/05_vcf/genotyping
