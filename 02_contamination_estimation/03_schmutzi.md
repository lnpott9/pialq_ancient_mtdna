# 03_schmutzi



## 05. Schmutzi (split by filtering & library prep)

"Unfiltered" will refer to the bam files that were created by subsetting Q20 reads from the sequences mapped to the full GRCh37 genome. "Filtered" will refer to bam files created by subsetting only primary, non-numt Q20 reads from the sequences mapped to the full GRCh38 genome. 

**Preparing BAM files**

I first copied over the unfiltered deduplicated bam (aka the Q20-filtered ones only) files to a new working folder:

```
# source of unfiltered bam files
/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/02_circmap/results/deduplication/*/*rmdup.bam

# new working directory for unfiltered files
/home/mnievesc/shared/projects/PIALQ/2025_sequencing5_analyses_lp/03_schmutzi/unfilt/
```

And the same for the filtered deduplicated bam (aka Q20-filtered, only primary mappings) files:

```
# source of filtered bam files
/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/02_circmap_filt/results/deduplication/*/*rmdup.bam

# new working directory for filtered files
/home/mnievesc/shared/projects/PIALQ/2025_sequencing5_analyses_lp/03_schmutzi/filt/
```

To prepare for the contamination estimation, I split the files into dsdna and ssdna folders, sorted the bam files, MD tagged them, and indexed them in the ssdna and dsdna folders:

```
ref=/home/mnievesc/shared/ref_seqs/rCRS.fasta
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/unfilt/dsdna
# working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/unfilt/ssdna
# working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/dsdna
# working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/ssdna
module load samtools/

# sort all bam files, then add MD tags to them, the index the sorted tagged files
for bam_file in "$working"/*rmdup.bam; do
	filename=$(basename "$bam_file" .mt_rmdup.bam)
	samtools sort "$bam_file" > "${filename}_sorted.bam"
	samtools calmd -b "${filename}_sorted.bam" $ref > "${filename}_sorted.md.bam"
	samtools index "${filename}_sorted.md.bam"
done
```



**Running *schmutzi* with no contaminant prediction (--notusepredC)**

The first step in the *schmutzi* pipeline is to run the program with no contaminant prediction to get a preliminary sense of contamination levels, followed by running it with contaminant prediction if contamination is more than a few percent.

Running contamination estimation with **no contaminant prediction** (dsdna):

```
#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/unfilt/dsdna
# working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/dsdna
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
    "106mt_sorted.md.bam"
    "107mt_sorted.md.bam"
    "108mt_sorted.md.bam"
    "109mt_sorted.md.bam"
    "110mt_sorted.md.bam"
    "111mt_sorted.md.bam"
    "112mt_sorted.md.bam"
    "113mt_sorted.md.bam"
    "114mt_sorted.md.bam"
    "115mt_sorted.md.bam"
    "116mt_sorted.md.bam"
    "117mt_sorted.md.bam"
    "118mt_sorted.md.bam"
    "119mt_sorted.md.bam"
    "120mt_sorted.md.bam"
    "121mt_sorted.md.bam"
    "122mt_sorted.md.bam"
    "123mt_sorted.md.bam"
    "124mt_sorted.md.bam"
    "125mt_sorted.md.bam"
    "126mt_sorted.md.bam"
    "127mt_sorted.md.bam"
    "128mt_sorted.md.bam"
    "129mt_sorted.md.bam"
    "130mt_sorted.md.bam"
    "131mt_sorted.md.bam"
    "132mt_sorted.md.bam"
    "133mt_sorted.md.bam"
    "134mt_sorted.md.bam"
    "135mt_sorted.md.bam"
    "136mt_sorted.md.bam"
    "19dsV2_sorted.md.bam"
    "33dsV2_sorted.md.bam"
    "59mt_sorted.md.bam"
    "60mt_sorted.md.bam"
    "61mt_sorted.md.bam"
    "62mt_sorted.md.bam"
    "63mt_sorted.md.bam"
    "64mt_sorted.md.bam"
    "66mt_sorted.md.bam"
    "67mt_sorted.md.bam"
    "68mt_sorted.md.bam"
    "69mt_sorted.md.bam"
    "70mt_sorted.md.bam"
    "71mt_sorted.md.bam"
    "72mt_sorted.md.bam"
    "73mt_sorted.md.bam"
    "74mt_sorted.md.bam"
    "75mt_sorted.md.bam"
    "76mt_sorted.md.bam"
    "77mt_sorted.md.bam"
    "78mt_sorted.md.bam"
    "79mt_sorted.md.bam"
    "80mt_sorted.md.bam"
    "81mt_sorted.md.bam"
    "82mt_sorted.md.bam"
    "84mt_sorted.md.bam"
    "85mt_sorted.md.bam"
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



Running contamination estimation with **no contaminant prediction** (ssdna):

```
#!/bin/bash -l
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-8

# 1. set path to reference sequences, contam database, & go to working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/unfilt/ssdna
# working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/ssdna
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
    "53mt_sorted.md.bam"
    "54mt_sorted.md.bam"
    "55mt_sorted.md.bam"
    "56mt_sorted.md.bam"
    "57mt_sorted.md.bam"
    "58mt_sorted.md.bam"
    "65mt_sorted.md.bam"
    "83mt_sorted.md.bam"
    "LQ30a.1_sorted.md.bam"
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



I created summaries using the following code:

```
text="uf_dsdna.txt"
# text="f_dsdna.txt"
# text="uf_ssdna.txt"
# text="f_ssdna.txt"

echo -e "Sample\tAverage\tLower\tUpper" > "$text"
for file in *_schmutzi/*npred_final.cont.est; do
    sample=$(basename "$file" _schmutzi_npred_final.cont.est)
    awk -v sample="$sample" '{print sample "\t" $0}' "$file" >> "$text"
done
```

The "problem" samples were 56mt, 57mt, 58mt, 60mt, 78mt, and 106mt:

```
f_dsdna.txt: 60mt/106mt can't run due to low coverage, 78mt is contaminated
uf_dsdna.txt:

f_ssdna: 57mt can't run due to low coverage, 56mt/58mt are contaminated
uf_ssdna.txt: 57mt can't run due to low coverage, 56mt/58mt are contaminated
```

57mt, 60mt, and 106mt could not run due to low coverage, and 56mt/58mt/78mt/19sV2 came back as contaminated. Gabriel Renaud highly suggests disregarding the contamination estimates for samples with less than 10x and doubting the estimate of those between 10-20X, so I checked coverage before deciding which samples to decontaminate and rerun with `samtools coverage`:

**19dsV2**

```
# filt
% samtools coverage 19dsV2_sorted.md.bam | cut -f6-7
coverage        meandepth
92.9688 12.2502

# unfilt
% samtools coverage 19dsV2_sorted.md.bam | cut -f6-7
coverage        meandepth
98.298  16.3852
```



**56mt**

```
# filt
% samtools coverage 56mt_sorted.md.bam | cut -f6-7
coverage        meandepth
99.994  200.814

# unfilt
% samtools coverage 56mt_sorted.md.bam | cut -f6-7
coverage        meandepth
100     223.792
```

**58mt**

```
# filt
% samtools coverage 58mt_sorted.md.bam | cut -f6-7
coverage        meandepth
53.92   1.64156

# unfilt
% samtools coverage 58mt_sorted.md.bam | cut -f6-7
coverage        meandepth
63.7093 1.98268
```

**78mt**

```
# filt
% samtools coverage 78mt_sorted.md.bam | cut -f6-7
coverage        meandepth
98.8895 26.6705

# unfilt
% samtools coverage 78mt_sorted.md.bam | cut -f6-7
coverage        meandepth
99.831  33.7385
```

Based on this, I will choose to only run schmutzi with contaminant prediction on 56mt, 78mt, and 19dsV2. **Libraries 57mt, 58mt, 60mt, and 106mt will be excluded from downstream analyses because a consensus genome could not be called.** This means there can only be consensus fastas for a max of 63 files. Results for the rest of the libraries are here:



**f_ssdna**

```
Sample	Average	Lower	Upper
53mt	0.01	0	0.02
54mt	0.02	0.01	0.03
55mt	0.05	0.04	0.06
56mt	0.99	0.98	0.99
58mt	0.98	0.97	0.99
65mt	0.02	0.01	0.03
83mt	0.03	0.02	0.04
LQ30a.1	0.01	0	0.02

*57mt could not run due to low coverage
```

**uf_ssdna**

```
Sample	Average	Lower	Upper
53mt	0.01	0	0.02
54mt	0.02	0.01	0.03
55mt	0.04	0.03	0.05
56mt	0.76	0.75	0.77
58mt	0.97	0.96	0.98
65mt	0.02	0.01	0.03
83mt	0.02	0.01	0.03
LQ30a.1	0.01	0	0.02

*57mt could not run due to low coverage
```

**f_dsdna**

```
Sample	Average	Lower	Upper
107mt	0.02	0.01	0.03
108mt	0.01	0	0.02
109mt	0.01	0	0.02
110mt	0.01	0	0.02
111mt	0.01	0	0.02
112mt	0.01	0	0.02
113mt	0.02	0.01	0.03
114mt	0.01	0	0.02
115mt	0.01	0	0.02
116mt	0.01	0	0.02
117mt	0.02	0.01	0.03
118mt	0.02	0.01	0.03
119mt	0.02	0.01	0.03
120mt	0.01	0	0.02
121mt	0.02	0.01	0.03
122mt	0.01	0	0.02
123mt	0.02	0.01	0.03
124mt	0.02	0.01	0.03
125mt	0.01	0	0.02
126mt	0.01	0	0.02
127mt	0.02	0.01	0.03
128mt	0.01	0	0.02
129mt	0.02	0.01	0.03
130mt	0.01	0	0.02
131mt	0.02	0.01	0.03
132mt	0.01	0	0.02
133mt	0.01	0	0.02
134mt	0.02	0.01	0.03
135mt	0.01	0	0.02
136mt	0.02	0.01	0.03
19dsV2	0.99	0.98	0.99
33dsV2	0.03	0.02	0.04
59mt	0.03	0.02	0.04
61mt	0.01	0	0.02
62mt	0.01	0	0.02
63mt	0.02	0.01	0.03
64mt	0.01	0	0.02
66mt	0.01	0	0.02
67mt	0.02	0.01	0.03
68mt	0.01	0	0.02
69mt	0.01	0	0.02
70mt	0.01	0	0.02
71mt	0.01	0	0.02
72mt	0.01	0	0.02
73mt	0.02	0.01	0.03
74mt	0.01	0	0.02
75mt	0.01	0	0.02
76mt	0.02	0.01	0.03
77mt	0.01	0	0.02
78mt	0.99	0.98	0.99
79mt	0.01	0	0.02
80mt	0.02	0.01	0.03
81mt	0.02	0.01	0.03
82mt	0.01	0	0.02
84mt	0.01	0	0.02
85mt	0.03	0.02	0.04

*60mt and 106mt could not run due to low coverage
```

**uf_dsdna**

```
Sample	Average	Lower	Upper
107mt	0.02	0.01	0.03
108mt	0.01	0	0.02
109mt	0.01	0	0.02
110mt	0.01	0	0.02
111mt	0.01	0	0.02
112mt	0.01	0	0.02
113mt	0.01	0	0.02
114mt	0.01	0	0.02
115mt	0.02	0.01	0.03
116mt	0.01	0	0.02
117mt	0.02	0.01	0.03
118mt	0.01	0	0.02
119mt	0.02	0.01	0.03
120mt	0.01	0	0.02
121mt	0.01	0	0.02
122mt	0.01	0	0.02
123mt	0.01	0	0.02
124mt	0.01	0	0.02
125mt	0.01	0	0.02
126mt	0.01	0	0.02
127mt	0.01	0	0.02
128mt	0.01	0	0.02
129mt	0.01	0	0.02
130mt	0.01	0	0.02
131mt	0.01	0	0.02
132mt	0.01	0	0.02
133mt	0.01	0	0.02
134mt	0.02	0.01	0.03
135mt	0.01	0	0.02
136mt	0.01	0	0.02
19dsV2	0.03	0.02	0.04
33dsV2	0.02	0.01	0.03
59mt	0.02	0.01	0.03
61mt	0.01	0	0.02
62mt	0.01	0	0.02
63mt	0.02	0.01	0.03
64mt	0.01	0	0.02
66mt	0.01	0	0.02
67mt	0.02	0.01	0.03
68mt	0.01	0	0.02
69mt	0.01	0	0.02
70mt	0.01	0	0.02
71mt	0.01	0	0.02
72mt	0.01	0	0.02
73mt	0.02	0.01	0.03
74mt	0.01	0	0.02
75mt	0.01	0	0.02
76mt	0.03	0.02	0.04
77mt	0.01	0	0.02
78mt	0.16	0.15	0.17
79mt	0.01	0	0.02
80mt	0.01	0	0.02
81mt	0.02	0.01	0.03
82mt	0.01	0	0.02
84mt	0.01	0	0.02
85mt	0.02	0.01	0.03

*60mt and 106mt could not run due to low coverage
```



