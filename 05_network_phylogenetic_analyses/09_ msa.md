# 09. MSA of consensus sequences

conda remake:

```
module load miniforge

conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env python numpy pandas biopython seaborn matplotlib -c conda-forge -c bioconda

source activate /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/msa_env
```



## 01. Path to consensus fasta sequences

I chose to use the more highly filtered sequences moving forward. I copied the consensus sequences produced by schmutzi to `$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_pialq`. 

I removed their "f_" prefix with the following code:

```
for file in f_*fasta; do
	basename=$(basename "$file" .fasta)
	name=${basename#f_}
	mv "$file" "$name.fasta"
done
```



## 02. Compare sequences for same individual

I also needed to decide which sequence would be used for individuals LQ02, LQ08, LQ12, LQ24, and LQ30, since each of these individuals had 2 high-quality fasta sequences representing them. To choose, I compared the sequences for each individual using the following script. If one sequence had an N at a position, I used the base at the same position in the other sequence to "fill in."

```
#!/usr/bin/env python3

import sys

# function to read in file & assign header & sequence
def read_fasta(filename):
    header = ""
    sequence = ""
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line
            else:
                sequence += line
    
    return header, sequence

# function to write new fasta file with same header & sequence
def write_fasta(filename, header, sequence):
    with open(filename, 'w') as f:
        f.write(header + '\n')
        f.write(sequence + '\n')

# function to loop through sequence 1 & compare to sequence 2
# if seq 1 has an N, use letter from other sequence
# if seq2 has an N, use letter from other sequence

def merge_sequences(seq1, seq2):
    merged_seq1 = []
    merged_seq2 = []
    
    for i in range(len(seq1)):
        base1 = seq1[i]
        base2 = seq2[i]
        
        if base1 == 'N' and base2 != 'N':
            merged_seq1.append(base2)
        else:
            merged_seq1.append(base1)
        
        if base2 == 'N' and base1 != 'N':
            merged_seq2.append(base1)
        else:
            merged_seq2.append(base2)
    
    return ''.join(merged_seq1), ''.join(merged_seq2)

file1 = sys.argv[1]
file2 = sys.argv[2]

# read in file 1
header1, seq1 = read_fasta(file1)

# read in file 2
header2, seq2 = read_fasta(file2)

merged_seq1, merged_seq2 = merge_sequences(seq1, seq2)

write_fasta(f"{file1.replace('.fasta', '')}_merged.fasta", header1, merged_seq1)
write_fasta(f"{file2.replace('.fasta', '')}_merged.fasta", header2, merged_seq2)
```

I then ran the following to generate a consensus fasta for each sequence, filling in its Ns with the base from the other sequence:

```
# LQ02
python merge.py 83mt.fasta 84mt.fasta

# LQ08
python merge.py 55mt.fasta 74mt.fasta

# LQ12
python merge.py 54mt.fasta 70mt.fasta

# LQ24
python merge.py 56mt_pmd3.fasta 78mt_pmd3.fasta

# LQ30
python merge.py 77mt.fasta LQ30a.1.fasta
```

I looked at the results of this using:

```
for file in *.fasta; do
    n_count=$(grep -v "^>" "$file" | grep -o "N" | wc -l)
    echo "$file: $n_count ns"
done
```



```
# LQ02
83mt.fasta: 6103 ns
83mt_merged.fasta: 4142 ns
84mt.fasta: 4608 ns
84mt_merged.fasta: 4142 ns

# LQ08
55mt.fasta: 1924 ns
55mt_merged.fasta: 1562 ns
74mt.fasta: 1963 ns
74mt_merged.fasta: 1562 ns

# LQ12
54mt.fasta: 178 ns
54mt_merged.fasta: 178 ns
70mt.fasta: 1225 ns
70mt_merged.fasta: 178 ns

# LQ24
56mt_pmd3.fasta: 1958 ns
56mt_pmd3_merged.fasta: 1952 ns
78mt_pmd3.fasta: 9690 ns
78mt_pmd3_merged.fasta: 1952 ns

# LQ30
77mt.fasta: 1238 ns
77mt_merged.fasta: 1214 ns
LQ30a.1.fasta: 2093 ns
LQ30a.1_merged.fasta: 1214 ns
```

The resulting files all had identical numbers of unknown bases because the missing/known bases were complementary. For consistency, I chose the file with the lower starting number of Ns as the merge file to use downstream, meaning I used

```
# LQ02
84mt_merged.fasta

# LQ08
55mt_merged.fasta

# LQ12
54mt_merged.fasta

# LQ24
56mt_pmd3_merged.fasta

# LQ30
77mt_merged.fasta: 1214 ns
```



Complete list of input PIALQ fasta files with commands to rename them:

```
% ls -1
107mt.fasta
108mt.fasta
109mt.fasta
110mt.fasta
111mt.fasta
112mt.fasta
113mt.fasta
114mt.fasta
115mt.fasta
116mt.fasta
117mt.fasta
118mt.fasta
119mt.fasta
120mt.fasta
121mt.fasta
122mt.fasta
123mt.fasta
124mt.fasta
125mt.fasta
126mt.fasta
127mt.fasta
128mt.fasta
129mt.fasta
130mt.fasta
131mt.fasta
132mt.fasta
133mt.fasta
134mt.fasta
135mt.fasta
136mt.fasta
19dsV2.fasta
33dsV2.fasta
53mt.fasta
54mt_merged.fasta
55mt_merged.fasta
56mt_pmd3_merged.fasta
59mt.fasta
61mt.fasta
62mt.fasta
63mt.fasta
64mt.fasta
65mt.fasta
66mt.fasta
67mt.fasta
68mt.fasta
69mt.fasta
71mt.fasta
72mt.fasta
73mt.fasta
75mt.fasta
76mt.fasta
77mt_merged.fasta
79mt.fasta
80mt.fasta
81mt.fasta
82mt.fasta
84mt_merged.fasta
85mt.fasta

# renaming
cp 107mt.fasta LQ48.fasta
cp 108mt.fasta LQ50.fasta
cp 109mt.fasta LQ44.fasta
cp 110mt.fasta LQ61.fasta
cp 111mt.fasta LQ64.fasta
cp 112mt.fasta LQ43.fasta
cp 113mt.fasta LQ46.fasta
cp 114mt.fasta LQ54.fasta
cp 115mt.fasta LQ40.fasta
cp 116mt.fasta LQ52.fasta
cp 117mt.fasta LQ59.fasta
cp 118mt.fasta LQ34.fasta
cp 119mt.fasta LQ39.fasta
cp 120mt.fasta LQ42.fasta
cp 121mt.fasta LQ57.fasta
cp 122mt.fasta LQ36.fasta
cp 123mt.fasta LQ63.fasta
cp 124mt.fasta LQ31.fasta
cp 125mt.fasta LQ35.fasta
cp 126mt.fasta LQ45.fasta
cp 127mt.fasta LQ33.fasta
cp 128mt.fasta LQ49.fasta
cp 129mt.fasta LQ38.fasta
cp 130mt.fasta LQ55.fasta
cp 131mt.fasta LQ62.fasta
cp 132mt.fasta LQ32.fasta
cp 133mt.fasta LQ41.fasta
cp 134mt.fasta LQ56.fasta
cp 135mt.fasta LQ58.fasta
cp 136mt.fasta LQ60.fasta
cp 19dsV2.fasta LQ19.fasta
cp 33dsV2.fasta LQ15.fasta
cp 53mt.fasta LQ11.fasta
cp 54mt_merged.fasta LQ12.fasta
cp 55mt_merged.fasta LQ08.fasta
cp 56mt_pmd3_merged.fasta LQ24.fasta
cp 59mt.fasta LQ21.fasta
cp 61mt.fasta LQ20.fasta
cp 62mt.fasta LQ10.fasta
cp 63mt.fasta LQ29.fasta
cp 64mt.fasta LQ28.fasta
cp 65mt.fasta LQ09.fasta
cp 66mt.fasta LQ16.fasta
cp 67mt.fasta LQ13.fasta
cp 68mt.fasta LQ03.fasta
cp 69mt.fasta LQ14.fasta
cp 71mt.fasta LQ18.fasta
cp 72mt.fasta LQ17.fasta
cp 73mt.fasta LQ23.fasta
cp 75mt.fasta LQ27.fasta
cp 76mt.fasta LQ25.fasta
cp 77mt_merged.fasta LQ30.fasta
cp 79mt.fasta LQ01.fasta
cp 80mt.fasta LQ22.fasta
cp 81mt.fasta LQ05.fasta
cp 82mt.fasta LQ06.fasta
cp 84mt_merged.fasta LQ02.fasta
cp 85mt.fasta LQ26.fasta
```

I also added all of these files to a folder within `$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta/pialq_ancient`. 



## 03. Metadata

I manually set up a spreadsheet with information on the following:

```
file	paper	ancient_or_modern	original_format	country	extra_info	mtdna_type
```

**file** = file name

**paper** = paper data is from, in the format <last name of first author>_<year>

**ancient_or_modern** = specifies whether sequence is ancient or modern DNA

**original_format** = original format of the data, options are fastq, mtdna_bam, wg_bam, mtdna_fasta, dloop_fasta

**country** = country individuals sampled are from

**extra_info** = for extra geographic, ethnic, or cultural identifiers authors used for sampled individuals

**mtdna_type** = full or dloop

I uploaded this spreadsheet as a tab-delimited text file called `all.meta`, only including the full mtDNA sequences. There are 18,778 sequences total.



## 04. Haplogrep input files

*04_hap.sh*

```
#!/bin/bash -l
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --tmp=2g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# 1. set paths
# parent folder is analysis folder
parent_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/"
temp="$parent_folder/tmp"
mkdir -p "$temp"

haplogrep="/projects/standard/mnievesc/shared/programs/haplogrep3"

# set paths to input fasta files
a_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta"
m_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta"

# set paths to output files & clear
summary_file="$parent_folder/all_hap_summ.txt"
short_summary_file="$parent_folder/short_all_hap_summ.txt"
> "$summary_file"
> "$short_summary_file"

# 2. add headings to full summary file
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$summary_file"

# 3. run haplogrep to generate the top hit
for folder in "$a_folder"/* "$m_folder"/*; do
    for fasta in "$folder"/*.fasta; do
        sample=$(basename "$fasta" .fasta)      
        output_tmp="$temp/${sample}.haplogrep.txt"
        "$haplogrep" classify \
            --tree phylotree-rcrs@17.2 \
            --in "$fasta" \
            --extend-report \
            --hits 1 \
            --out "$output_tmp"
        # extract second line from output and append to summary file
        if [ -s "$output_tmp" ]; then
            second_line=$(awk 'NR==2 {gsub(/"/, ""); print}' "$output_tmp")
            echo -e "$sample\t${second_line}" >> "$summary_file"
        else
            echo "⚠️⚠️ no output for $sample ⚠️⚠️"
        fi
    done
done

# 4. create short version of summary file
cut -f1-5 "$summary_file" > "$short_summary_file"
```



code for adding new sequences:

```
# 1. set paths
# parent folder is analysis folder
parent_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa"
temp="$parent_folder/tmp"

haplogrep="/projects/standard/mnievesc/shared/programs/haplogrep3"

# set paths to input fasta files
folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/bodner_2022_modern"

# set paths to output files & clear
summary_file="$parent_folder/all_hap_summ.txt"
short_summary_file="$parent_folder/short_all_hap_summ.txt"

# 3. run haplogrep to generate the top hit
for fasta in "$folder"/*.fasta; do
    sample=$(basename "$fasta" .fasta) 
    output_tmp="$temp/${sample}.haplogrep.txt"
    "$haplogrep" classify \
        --tree phylotree-rcrs@17.2 \
        --in "$fasta" \
        --extend-report \
        --hits 1 \
        --out "$output_tmp"
    # extract second line from output and append to summary file
    second_line=$(awk 'NR==2 {gsub(/"/, ""); print}' "$output_tmp")
    echo -e "$sample\t${second_line}" >> "$summary_file"
done

# 4. create short version of summary file
cut -f1-5 "$summary_file" > "$short_summary_file"
```



code for renaming pialq files in hap results:

```
full = pd.read_csv("all_hap_summ_updated.txt", sep="\t")

# separate pialq samples
pialq = full[full["Sample"].str.contains("mt", case=True)]

# delete those from text file & overwrite
updated_full = full[~full["Sample"].str.contains("mt", case=True)]

# get rid of ds samples too
pialq = updated_full[updated_full["Sample"].str.contains("dsV2", case=True)]

# write out as csv
updated_full.to_csv("all_hap_summ_updated.txt", sep="\t", index=False, header=True)
```



## 05. QC sequences

```
srun --mem=10gb --time=5:00:00 --pty bash

# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env

cd $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/

python

# 2. import pandas
import pandas as pd

# 3. read in short_all_haplogrep_summary.txt as dataframe
hg = pd.read_csv("short_all_hap_summ.txt", sep="\t")

# 4. look at head to make sure all is well
hg.head()

hg.shape
(18540, 5)

# 5. read in metadata file
md = pd.read_csv("all.meta", sep="\t")

md.shape
(18566, 7)

# 6. add new column to md df of sample name without fasta
md["Sample"] = md["file"].str.replace(".fasta", "", regex=False)

# 7. ensure haplogrep & metadata are same length
assert hg.shape[0] == md.shape[0]

# Samples in md that don't have haplogroup data
missing_in_hg = set(md["Sample"]) - set(hg["Sample"])
# print(f"Samples in md but not in hg: {len(missing_in_hg)}")
# print(missing_in_hg)
# Haplogroup data without metadata (none)
# missing_in_md = set(hg["Sample"]) - set(md["Sample"])
# print(f"Samples in hg but not in md: {len(missing_in_md)}")
# print(missing_in_md)
# remove these samples from md
md = md[~md["Sample"].isin(missing_in_hg)]
# or remove from hg
# hg = hg[~hg["Sample"].isin(missing_in_md)]
# assert hg.shape[0] == md.shape[0]

# 13. read in fasta files
from pathlib import Path
from Bio import SeqIO

# directories w fasta files
ancient_dir = Path("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta")
modern_dir = Path("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta")

# find all fasta files in both directories
fasta_files = list(ancient_dir.rglob("*.fasta")) + list(modern_dir.rglob("*.fasta"))

# Read all sequences into a dataframe
fasta_list = []
for fasta_file in fasta_files:
    df = pd.DataFrame([{
        "sequence_id": record.id,
        "description": record.description,
        "sequence": str(record.seq),
        "file": fasta_file.name
    } for record in SeqIO.parse(fasta_file, "fasta")])
    fasta_list.append(df)

# combine all dataframes
all_fasta = pd.concat(fasta_list, ignore_index=True)

# add sample name
all_fasta["Sample"] = all_fasta["file"].str.replace(".fasta", "", regex=False)

# 14. count number of characters (ATCG)
all_fasta["length"] = all_fasta["sequence"].str.upper().str.count(r"[ACGT]")
# look at distribution
all_fasta["length"].describe()
count    18541.000000
mean     16534.996171
std        518.195750
min         14.000000
25%      16564.000000
50%      16568.000000
75%      16569.000000
max      16659.000000
Name: length, dtype: float64

fasta_lengths = all_fasta[["Sample", "length"]]

# isolate first quartile (4725 seqs)
# lowest_quartile = fasta[fasta["length"] <= fasta["length"].quantile(0.25)]
# first shortest
# lowest_quartile_sorted = lowest_quartile.sort_values(by="length", ascending=True)

# 14. join cleaned haplogrep & metadata text files based on sample column
merged = pd.merge(md, hg, on="Sample", how="inner")

merged.shape
(18540, 12)

# replace Haplogroup for PIALQ C samples that were manually confirmed to be C1b
samples_to_update = ["LQ06", "LQ23", "LQ60", "LQ62"]

# update Haplogroup to C1b
merged.loc[merged["Sample"].isin(samples_to_update), "Haplogroup"] = "C1b"

# 15. add row to describe multiple character length haplogroups
def get_nchar(haplo, n):
    if pd.isna(haplo):
        return None
    return haplo[:n]
    
for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    merged[f"hap_{n}char"] = merged["Haplogroup"].apply(lambda x: get_nchar(x, n))

# add fasta lengths to df
# join fasta length with merged
missing_in_hg_id = set(merged["Sample"]) - set(all_fasta["Sample"])
missing_in_sq_id = set(all_fasta["Sample"]) - set(merged["Sample"])

combo = pd.merge(merged, fasta_lengths, on="Sample", how="inner")

# fix some of the incorrect metadata
combo['paper'] = combo['paper'].replace({
    'bodner_2023': 'bodner_2012',
    'valverde_2016': 'llamas_2016',
    'garciaolivares_2023': 'garciaolivares_2022',
    'barbieri_2012b': 'barbieri_2013a'
})

# add country labels
for country in sorted(one['country'].unique()):
    print(country)
    
# join with country_labels.txt
labels = pd.read_csv('country_labels.txt', sep='\t')
new = combo.merge(labels, on='country', how='left')

# change all that start with LQ as "Hacienda La Quebrada"
new.loc[new['file'].str.startswith('LQ'), 'country_label'] = 'Hacienda La Quebrada'
final = new

# save
final.to_csv("all_merged.meta", sep="\t", index=False, header=True)

# 16. remove anything that isnt A, B, C, D, H, L
subset = final[final["hap_2char"].isin(["A2", "B2", "C1", "D1", "H1", "L0", "L1", "L2", "L3"])].copy()

subset.shape # (10065, 21)

# 17. only keep samples longer than 14000 characters AND quality over 80%, or if it's a pialq sequence
quality = subset[((subset["length"] > 14000) & (subset["Quality"] > 0.8)) | (subset["paper"] == "pialq")].copy()

remove = subset[((subset["length"] < 14000) | (subset["Quality"] < 0.8) & (subset["paper"] != "pialq"))].copy()

remove.shape # (93, 21)

# write out low qual to a file
remove.to_csv("removed.meta", sep="\t", index=False, header=True)

# 18. only keep samples that are in the same haplogroup as mine
# A (target: A2 and A2+(64))
# B (target: B2 and B2b)
# C (target: C1 and C1b)
# D (target: D1)
# H (target: H1 and H1b)
# L0 (targets: L0a2a, L0d1a)
# L1 (targets: L1c, L1c1, L1c2'4)
# L2 (targets: L2a1, L2a1a3, L2a1d, L2a1f, L2a1q, L2b1, L2c
# L3 (targets: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b)

# define keep_list
keep_2 = ["A2", "B2", "C1", "D1", "H1", "L0", "L1", "L2", "L3"]
keep_3 = ["B2b", "C1b", "L1c", "L2c"]
keep_4 = ["H1bw", "L1c1", "L1c2", "L1c4", "L2a1", "L2b1", "L3d1", "L3d3", "L3d4", "L3e2", "L3e3"]
keep_5 = ["L0a2a", "L0d1a", "L2a1d", "L2a1f", "L2a1q", "L3b1a", "L3f1b"]
keep_6 = ["L2a1a3"]
keep_7 = ["A2+(64)", "L3e1a1a"]

# 18. subset all by samples that have a desired haplogroup
small = quality[(quality["hap_2char"].isin(keep_2)) | (quality["hap_3char"].isin(keep_3)) | (quality["hap_4char"].isin(keep_4)) | (quality["hap_5char"].isin(keep_5)) | (quality["hap_6char"].isin(keep_6)) | (quality["hap_7char"].isin(keep_7))]

small.shape
(9987, 21)

# deduplicate
small_final = small.drop_duplicates("Sample", keep="first")

# write new meta (samples of interest & outgroups)
small_final.to_csv("final_cut.meta", sep="\t", index=False, header=True)
```

country_labels.txt:

```
country	country_label
ACB	African Caribbean in Barbados (ACB)
ASW	African Ancestry in SW USA (ASW)
AbkhasianSGDP	Abkhazia or Russia
AdygeiSGDP	Russia
AlbanianSGDP	Albania
AleutSGDP	Russia
Algeria	Algeria
AltaianSGDP	Russia
AmiSGDP	Taiwan
Angola	Angola
Argentina	Argentina
ArmenianSGDP	Armenia
AtayalSGDP	Taiwan
AustralianSGDP	Australia
BEB	 Bengali in Bangladesh (BEB)
BalochiSGDP	Pakistan
BantuHereroSGDP	Botswana or Namibia
BantuKenyaHGDP	Kenya
BantuKenyaSGDP	Kenya
BantuSouthAfricaHGDP	South Africa
BantuTswanaSGDP	Botswana or Namibia
BasqueSGDP	Spain
BedouinBSGDP	Israel
Belize	Belize
BengaliSGDP	Bangladesh
BergamoItalianHGDP	Italy
BergamoSGDP	Italy
BiakaHGDP	Central African Republic
BiakaSGDP	Central African Republic
Bolivia	Bolivia
Botswana	Botswana
BougainvilleSGDP	Papua New Guinea
BrahminSGDP	India
BrahuiSGDP	Pakistan
Brazil	Brazil
BulgarianSGDP	Bulgaria
Burkina Faso	Burkina Faso
BurmeseSGDP	Myanmar
BurushoSGDP	Pakistan
CDX	Chinese Dai in Xishuangbanna China (CDX)
CEU	Utah residents with Northern and Western European ancestry (CEU)
CHB	Han Chinese in Beijing, China (CHB)
CHS	Han Chinese South (CHS)
CLM	Colombian in Medellín, Colombia (CLM)
CambodianHGDP	Cambodia
CambodianSGDP	Cambodia
Cameroon	Cameroon
Canada	Canada
Canary Islands	Canary Islands
Chad	Chad
ChaneSGDP	Argentina
ChechenSGDP	Russia
Chile	Chile
ChukchiSGDP	Russia
Colombia	Colombia
ColombianHGDP	Colombia
Comoros	Comoros
CreteSGDP	Greece
Croatia	Croatia
CzechSGDP	Czechia
Czechia	Czechia
DaiHGDP	China
DaiSGDP	China
DaurHGDP	China
DaurSGDP	China
DinkaSGDP	Sudan
Dominican Republic	Dominican Republic
DruzeSGDP	Israel
DusunSGDP	Brunei
ESN	Esan in Nigeria (ESN)
Ecuador	Ecuador
EnglishSGDP	England
Equatorial Guinea	Equatorial Guinea
EsanSGDP	Nigeria
EskimoChaplinSGDP	Russia
EskimoNaukanSGDP	Russia
EskimoSirenikiSGDP	Russia
EstonianSGDP	Estonia
Ethiopia	Ethiopia
EvenSGDP	Russian
FIN	Finnish in Finland (FIN)
FinnishSGDP	Finland
France	France
FrenchHGDP	France
FrenchSGDP	France
GBR	British from England and Scotland (GBR)
GIH	Gujarati Indians in Houston, Texas, USA (GIH)
GWD	Gambian in Western Division – Mandinka (GWD)
GambianSGDP	The Gambia
GeorgianSGDP	Georgia
Germany	Germany
Ghana	Ghana
GreekSGDP	Greece
Guinea Equatorial	Equatorial Guinea
HanHGDP	China
HanSGDP	China
HawaiianSGDP	Hawaii
HazaraSGDP	Pakistan
HezhenSGDP	China
HungarianSGDP	Hungary
IBS	Iberian populations in Spain (IBS)
ITU	Indian Telugu in the UK (ITU)
IcelandicSGDP	Iceland
IgorotSGDP	Philippines
IranianSGDP	Iran
IraqiJewSGDP	Iraq
IrulaSGDP	India
Italy	Italy
ItelmanSGDP	Russia
Ivory Coast	Côte d'Ivoire
JPT	Japanese in Tokyo, Japan (JPT)
JapaneseHGDP	Japan
JapaneseSGDP	Japan
Jordan	Jordan
JordanianSGDP	Jordan
JuhoanNorthSGDP	Namibia
KHV	Kinh in Ho Chi Minh City, Vietnam (KHV)
KalashSGDP	Pakistan
KapuSGDP	India
KaritianaHGDP	Brazil
KaritianaSGDP	Brazil
Kenya	Kenya
KhomaniSanSGDP	South Africa
KhondaDoraSGDP	India
KinhSGDP	Vietnam
KoreanSGDP	South Korea
KusundaSGDP	Nepal
KyrgyzSGDP	Kyrgyzystan
LWK	Luhya in Webuye, Kenya (LWK)
La Gomera	Canary Islands
LahuHGDP	China
LahuSGDP	China
Lanzarote	Canary Islands
LezginSGDP	Russia
Libya	Libya
LuhyaSGDP	Kenya
LuoSGDP	Kenya
MSL	Mende in Sierra Leone (MSL)
MXL	Mexican Ancestry in Los Angeles CA USA (MXL)
Madagascar	Madagascar
MadigaSGDP	India
MakraniHGDP	Pakistan
MakraniSGDP	Pakistan
MalaSGDP	India
Mali	Mali
MandenkaHGDP	Senegal
MandenkaSGDP	Senegal
MansiSGDP	Russia
MaoriSGDP	New Zealand
MasaiSGDP	Kenya
Mauritania	Mauritania
MayanSGDP	Mexico
MbutiHGDP	Congo
MbutiSGDP	Congo
MendeSGDP	Sierra Leone
Mexico	Mexico
MiaoHGDP	China
MiaoSGDP	China
MixeSGDP	Mexico
MixtecSGDP	Mexico
MongolaSGDP	China
MongolianHGDP	Mongolia
Morocco	Morocco
MozabiteHGDP	Algeria
MozabiteSGDP	Algeria
Mozambique	Mozambique
Namibia	Namibia
NaxiHGDP	China
NaxiSGDP	China
Niger	Niger
Nigeria	Nigeria
NorthOssetianSGDP	Russian
NorthernHanHGDP	China
NorwegianSGDP	Norway
OrcadianSGDP	Orkney Islands
Oromo	Ethiopia
OroqenHGDP	China
OroqenSGDP	China
PEL	Peruvian in Lima, Peru (PEL)
PJL	Punjabi in Lahore, Pakistan (PJL)
PUR	Puerto Rican in Puerto Rico (PUR)
PalestinianSGDP	Israel
Panama	Panama
PapuanSGDP	Papua New Guinea
Paraguay	Paraguay
PathanSGDP	Pakistan
Peru	Peru
PiapocoSGDP	Mexico
PimaSGDP	Mexico
PolishSGDP	Poland
Portugal	Portugal
PunjabiSGDP	India
QuechuaSGDP	Peru
RelliSGDP	India
Russia	Russia
RussianSGDP	Russia
Rwanda	Rwanda
STU	Sri Lankan Tamil in the UK (STU)
SaamiSGDP	Finland
SaharawiSGDP	Western Sahara
SamaritanSGDP	Israel
SanHGDP	South Africa
Sao Tome e Principe	São Tomé and Príncipe
SardinianHGDP	Italy
SardinianSGDP	Italy
Saudi Arabia	Saudi Arabia
Senegal	Senegal
SheHGDP	China
SheSGDP	China
SindhiSGDP	India
SomaliSGDP	Kenya
Somalia	Somalia
South Africa	South Africa
South Korea	South Korea
Spain	Spain
SpanishSGDP	Spain
St. Helena	St. Helena
St. Martin	St. Martin
Sudan	Sudan
SuruiHGDP	Brazil
SuruiSGDP	Brazil
São Tomé and Príncipe	São Tomé and Príncipe
TSI	Toscani in Italia (TSI)
TajikSGDP	Tajikistan
Tanzania	Tanzania
Tenerife	Canary Islands
ThaiSGDP	Thailand
TlingitSGDP	Russia
TuHGDP	China
TuSGDP	China
TubalarSGDP	Russia
TujiaHGDP	China
TujiaSGDP	China
Tunisia	Tunisia
TurkishSGDP	Turkey
TuscanHGDP	Italy
TuscanSGDP	Italy
USA	USA
UlchiSGDP	Russia
United States	USA
Uruguay	Uruguay
UygurHGDP	China
UygurSGDP	China
Venezuela	Venezuela
XiboHGDP	China
XiboSGDP	China
YRI	Yoruba in Ibadan, Nigeria (YRI)
YadavaSGDP	India
YakutSGDP	Russia
Yemen	Yemen
YemeniteJewSGDP	Yemen
YiHGDP	China
YiSGDP	China
YorubaHGDP	Nigeria
YorubaSGDP	Nigeria
Zambia	Zambia
ZapotecSGDP	Mexico
unknown	unknown
```



## 06. Sort files into alignment folders by 2-character haplogroup

interactively,

```
cd $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/
python 06_sort.py
```



```
# 06_sort.py

import os
import shutil
import pandas as pd

# 1. set paths
metadata_file = "final_cut.meta"

ancient_folder = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta"
modern_folder = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta"
destination_base = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft"

# 2. load metadata
df = pd.read_csv(metadata_file, sep="\t", dtype=str).fillna("")

# 3.read in each row of file, strip whitespace
for idx, row in df.iterrows():
    file = row["file"].strip()
    paper = row["paper"].strip()
    ancient_or_modern = row["ancient_or_modern"].strip().lower()
    haplo = row["hap_3char"].strip()

    # Determine source path
    if ancient_or_modern == "ancient":
        folder = f"{paper}_ancient"
        src = os.path.join(ancient_folder, folder, file)
    elif ancient_or_modern == "modern":
        folder = f"{paper}_modern"
        src = os.path.join(modern_folder, folder, file)

    # Destination folder based on Macrohaplogroup
    dest_folder = os.path.join(destination_base, haplo)
    os.makedirs(dest_folder, exist_ok=True)

    dest_file = os.path.join(dest_folder, file)

    # Copy or warn
    if os.path.isfile(src):
        shutil.copy2(src, dest_file)
    else:
        print(f"missing file: {src}")
```

output will be:

```
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669265.fasta
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669269.fasta
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669285.fasta
```

This is fine , as I deleted these three files since they were corrupted

I added outgroups by hand:

```
# neanderthal
efetch -db nucleotide -id AM948965 -format fasta > AM948965.fasta

# rsrs, supplemental file 2, named RSRS.fasta
# https://www.cell.com/AJHG/fulltext/S0002-9297(12)00146-2

# copy to alignment folders
BASE="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft"
for dir in $BASE/*/; do
cp $SHARED/ref_seqs/AM948965.fasta "$dir"
done
```



## 07. Reformat fasta headers

```
# set path to sorted fasta folder
sorted="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft"

for folder in "$sorted"/*/; do
python 07_rename.py "$folder"
done
```



```
# 07_rename.py

import os
import sys

# 1. if # arguments isn't two, fail
if len(sys.argv) != 2:
    print("Usage: 07_rename.py <fasta_directory>")
    sys.exit(1)

# 2. get directory containing fasta files from second argument
fasta_directory = sys.argv[1]

# 3. process each fasta file in directory
for fasta_file in os.listdir(fasta_directory):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(fasta_directory, fasta_file)
        new_fasta = []
        
        # extract filename without the extension
        fasta_name = os.path.splitext(fasta_file)[0]
    
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    new_fasta.append(f">{fasta_name}\n")
                else:
                    new_fasta.append(line)
        
        # overwrite original fasta
        with open(fasta_path, 'w') as f:
            f.writelines(new_fasta)
```



## 08. align all sequences

mafft how to: https://academic.oup.com/bib/article/20/4/1160/4106928

auto:

```
#!/bin/bash -l
#SBATCH --job-name=auto
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80gb
#SBATCH --tmp=80gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-8

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa"
SORTED="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# Manual folder listing
FOLDERS=(
    "A2"
    "B2"
    "C1"
    "D1"
    "H1"
    "L0"
    "L1"
    "L2"
    "L3"
)

# Folder for this task
FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID]}"
FULL_PATH="$SORTED/$FOLDER"

echo "processing folder $FOLDER at $FULL_PATH"

# Build input file
cat "$REF/rCRS.fasta" $FULL_PATH/*fasta | awk 'NF' > "$SORTED/${FOLDER}_mafft_in_auto.fasta"
echo "rcrs added"

# real mafft with auto
if mafft --auto \
    --thread 10 \
    "$SORTED/${FOLDER}_mafft_in_auto.fasta" \
    > "$SORTED/${FOLDER}_mafft_out_auto.fasta" ; then
    echo "brethren rejoice 💀"
else
	echo "once again we're fucked"
fi
```



## 09. clean MSA

on ondemand rstudio session

```
# clean.r
# load libraries
if (!require(pegas)) {
  install.packages("pegas")
  library(pegas)
}

if (!require(ape)) {
  install.packages("ape")
  library(ape)
}

if (!require(Biostrings)) {
  install.packages("Biostrings")
  library(Biostrings)
}

# set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean")

# 1. set paths to MSAs
a2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/A2_mafft_out.fasta")
b2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/B2_mafft_out.fasta")
c1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/C1_mafft_out.fasta")
d1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/D1_mafft_out.fasta")
h1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/H1_mafft_out.fasta")
l0_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L0_mafft_out.fasta")
l1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L1_mafft_out.fasta")
l2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L2_mafft_out.fasta")
l3_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L3_mafft_out.fasta")

# 2. read paths into a list
path_list <- list(a2_path, b2_path, c1_path, d1_path, h1_path, l0_path, l1_path, l2_path, l3_path)

# 3. read sequences into DNAbin objects, convert all to upper case characters
aligned_seqs <- lapply(path_list, read.dna, format = "fasta")

aligned_seqs_char <- lapply(aligned_seqs, function(x) toupper(as.character(x)))

# length/width should be pretty close to 16569 (length of rCRS)
lapply (aligned_seqs_char, function (x) dim(x))

# 4. identify rCRS sequence in each alignment
rCRS_seq <- lapply(aligned_seqs_char, function(x) {
  rcs_name <- grep("NC_012920.1", rownames(x), value = TRUE)
  x[rcs_name, ]
})

# 5. clean alignments: remove rCRS gaps, mask problem positions, replace ambiguous codes
problem_positions <- c(303:315, 515:522, 568:573, 3107, 16182:16194, 16519)
# ambiguous_codes <- c("R","Y","M","W","K","S","B","H","D","V")
haplogroups <- c("A2","B2","C1","D1","H1","L0","L1","L2","L3")

aligned_seqs_no_gaps <- lapply(seq_along(aligned_seqs_char), function(i) {
  alignment_mat <- aligned_seqs_char[[i]]
  rcs <- rCRS_seq[[i]]
  
  # mask problem positions
  rcs_coord <- 0
  mask_cols <- c()
  for (j in seq_along(rcs)) {
    if (rcs[j] != "-") {
      rcs_coord <- rcs_coord + 1
      if (rcs_coord %in% problem_positions) {
        mask_cols <- c(mask_cols, j)
      }
    }
  }
  if (length(mask_cols) > 0) {
    alignment_mat[, mask_cols] <- "N"
  }
  
  
  # replace ambiguous IUPAC codes with "N"
  #mat[mat %in% ambiguous_codes] <- "N"
  
  alignment_mat
})

names(aligned_seqs_no_gaps) <- haplogroups

# 6. Calculate missingness after cleaning
count_missing <- function(seq_matrix) {
  apply(seq_matrix, 1, function(row) sum(row == "N"))
}

missing_counts_list <- lapply(aligned_seqs_no_gaps, count_missing)

missing_df_list <- lapply(seq_along(aligned_seqs_no_gaps), function(i) {
  data.frame(
    Haplogroup = names(aligned_seqs_no_gaps)[i],
    Sample = rownames(aligned_seqs_no_gaps[[i]]),
    Missing_Count = missing_counts_list[[i]],
    Missing_Fraction = missing_counts_list[[i]] / ncol(aligned_seqs_no_gaps[[i]])
  )
})
missing_df <- do.call(rbind, missing_df_list)

# 7. Filter sequences based on missingness cutoff, set to keep in worst pialq sample
# missingness_cutoff <- 0.35
# samples_to_keep <- rownames(aligned_seqs_char)[missing_df$missing_fraction <= missingness_cutoff]
# filtered_aligned_seqs <- aligned_seqs_char[samples_to_keep, , drop = FALSE]

# 8. convert back to DNAStringSet
clean_alignments_list <- lapply(aligned_seqs_no_gaps, function(mat) {
  # Convert each row to a DNAString
  dna_seq <- apply(mat, 1, function(seq_row) {
    DNAString(paste(seq_row, collapse = ""))
  })
  
  # Assign rownames as sequence names
  names(dna_seq) <- rownames(mat)
  
  # Wrap as DNAStringSet
  DNAStringSet(dna_seq)
})

# 12. Remove the rCRS reference (NC_012920.1) from each
clean_alignments_list <- lapply(clean_alignments_list, function(dna_set) {
  names_to_drop <- grep("NC_012920.1", names(dna_set), value = TRUE)
  dna_set[!names(dna_set) %in% names_to_drop]
})

# 13. write cleaned alignments
haplogroups <- names(aligned_seqs_no_gaps)
for (i in seq_along(clean_alignments_list)) {
  writeXStringSet(
    clean_alignments_list[[i]],
    filepath = paste0("clean_", haplogroups[i], ".fasta"),
    width = 18000
  )
}

```



## 10. write cleaned alignments to tree files

splitting scripts at 

```
/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean
```



splitting for iqtree:

```
# split_iqtree.r
# 1. load libraries & install if necessary
if (!require(pegas)) {
  install.packages("pegas")
  library(pegas)
}

if (!require(ape)) {
  install.packages("ape")
  library(ape)
}

if (!require(Biostrings)) {
  install.packages("Biostrings")
  library(Biostrings)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

# 2. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree")

# 3. set paths
a2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_A2.fasta")
b2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_B2.fasta")
c1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_C1.fasta")
d1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_D1.fasta")
h1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_H1.fasta")
l0_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L0.fasta")
l1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L1.fasta")
l2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L2.fasta")
l3_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L3.fasta")

# 4. read fastas as dnabin objects
read_in <- function(path) {
  fasta <- read.dna(path, format = "fasta")
  return(fasta)
}

seqs_a2 <- read_in(a2_path)
seqs_b2 <- read_in(b2_path)
seqs_c1 <- read_in(c1_path)
seqs_d1 <- read_in(d1_path)
seqs_h1 <- read_in(h1_path)
seqs_l0 <- read_in(l0_path)
seqs_l1 <- read_in(l1_path)
seqs_l2 <- read_in(l2_path)
seqs_l3 <- read_in(l3_path)

# 5. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/final_cut.meta", sep="\t")

# 6. add 2 roots to the metadata so they're included in all subsets
ref_rows <- data.frame(
  file = "AM948965.fasta",
  paper = "REF",
  ancient_or_modern = NA,
  original_format = NA,
  country = NA,
  extra_info = "REF",
  mtdna_type = NA,
  Sample = "AM948965",
  HaplogrepSampleID = "REF",
  Haplogroup = "REF",
  Rank = NA,
  Quality = NA,
  hap_1char = "REF",
  hap_2char = "REF",
  hap_3char = "REF",
  hap_4char = "REF",
  hap_5char = "REF",
  hap_6char = "REF",
  hap_7char = "REF",
  hap_8char = "REF",
  length = NA
)

meta <- rbind(meta, ref_rows)

# 7. create objects for each tree subset
subset_dna <- function(dna_obj, meta_df, filter_expr) {
  # Evaluate the filter expression within meta_df
  rows <- meta_df %>% filter(!!enquo(filter_expr))
  # Find row indices in DNAbin
  idx <- match(rows$Sample, rownames(dna_obj))
  # Warn if any samples are missing
  if(any(is.na(idx))) {
    warning("sequences in metadata but missing in tree: ", paste(rows$Sample[is.na(idx)], collapse = ", "))
  }
  # Subset and return DNAbin (only existing sequences)
  dna_obj[idx[!is.na(idx)], ]
}

# A2: A2, A2+(64)
final_a2 <- subset_dna(
  seqs_a2,
  meta,
  (hap_2char == "A2" & !hap_3char %in% c("A20","A21","A22","A23","A24","A25","A26")) |
    hap_2char == "REF"
)

final_a264 <- subset_dna(seqs_a2, meta,(hap_7char == "A2+(64)" | hap_2char == "REF"))

# B2: B2, B2b
final_b2 <- subset_dna(seqs_b2, meta, hap_2char == "B2" | hap_2char == "REF")
final_b2b <- subset_dna(seqs_b2, meta, hap_3char == "B2b" | hap_2char == "REF" | Sample == "LQ14")

# C1: C1b
final_c1 <- subset_dna(seqs_c1, meta, hap_2char == "C1" | hap_2char == "REF")
final_c1b <- subset_dna(seqs_c1, meta, hap_3char == "C1b" | hap_2char == "REF")

# D1: D1
final_d1 <- subset_dna(seqs_d1, meta, hap_2char == "D1" | hap_2char == "REF")

# H1: H1, H1bw
final_h1 <- subset_dna(seqs_h1, meta, hap_3char %in% c("REF", "H1a", "H1b", "H1c", "H1e", "H1f", "H1g", "H1h", "H1i", "H1j", "H1k", "H1l", "H1m", "H1n", "H1o", "H1p", "H1q", "H1r", "H1s", "H1t", "H1u", "H1v", "H1w", "H1x", "H1y", "H1z"))
final_h1bw <- subset_dna(seqs_h1, meta, hap_4char == "H1bw" | hap_2char == "REF")

# L0: L0a2a, L0d1a
final_l0 <- subset_dna(seqs_l0, meta, hap_2char == "L0" | hap_2char == "REF")
final_l0a2a <- subset_dna(seqs_l0, meta, hap_5char == "L0a2a" | hap_2char == "REF")
final_l0d1a <- subset_dna(seqs_l0, meta, hap_5char == "L0d1a" | hap_2char == "REF")

# L1: L1c, L1c1, L1c2'4
final_l1 <- subset_dna(seqs_l1, meta, hap_2char == "L1" | hap_2char == "REF")
final_l1c1 <- subset_dna(seqs_l1, meta, hap_4char == "L1c1" | hap_2char == "REF")
final_l1c1a <- subset_dna(seqs_l1, meta, hap_5char == "L1c1a" | hap_2char == "REF")
final_l1c1b <- subset_dna(seqs_l1, meta, hap_5char == "L1c1b" | hap_2char == "REF")
final_l1c1d <- subset_dna(seqs_l1, meta, hap_5char == "L1c1d" | hap_2char == "REF")
final_l1c2 <- subset_dna(seqs_l1, meta, (hap_4char == "L1c2" & hap_5char != "L1c2'") | hap_2char == "REF")
final_l1c4 <- subset_dna(seqs_l1, meta, hap_4char == "L1c4" | hap_2char == "REF")
final_l1c24 <- subset_dna(seqs_l1, meta, (hap_4char == "L1c2" | hap_4char == "L1c4") | hap_2char == "REF")

# L2: L2a1, L2a1a, L2a1d, L2a1f, L2a1q, L2b1, L2c
final_l2 <- subset_dna(seqs_l2, meta, hap_2char == "L2" | hap_2char == "REF")
final_l2a1a <- subset_dna(seqs_l2, meta, hap_5char == "L2a1a" | hap_2char == "REF")
final_l2a1d <- subset_dna(seqs_l2, meta, hap_5char == "L2a1d" | hap_2char == "REF")
final_l2a1f <- subset_dna(seqs_l2, meta, hap_5char == "L2a1f" | hap_2char == "REF")
final_l2a1q <- subset_dna(seqs_l2, meta, hap_5char == "L2a1q" | hap_2char == "REF")
final_l2b1 <- subset_dna(seqs_l2, meta, hap_4char == "L2b1" | hap_2char == "REF")
final_l2c <- subset_dna(seqs_l2, meta, hap_3char == "L2c" | hap_2char == "REF")

# L3: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b
final_l3 <- subset_dna(seqs_l3, meta, hap_2char == "L3" | hap_2char == "REF")
final_l3b1a <- subset_dna(seqs_l3, meta, hap_5char == "L3b1a" | hap_2char == "REF")
final_l3d1 <- subset_dna(seqs_l3, meta, hap_4char == "L3d1" | hap_2char == "REF")
final_l3d3 <- subset_dna(seqs_l3, meta, hap_4char == "L3d3" | hap_2char == "REF")
final_l3d4 <- subset_dna(seqs_l3, meta, hap_4char == "L3d4" | hap_2char == "REF")
final_l3e1a1a <- subset_dna(seqs_l3, meta, hap_7char == "L3e1a1a" | hap_2char == "REF")
final_l3e1 <- subset_dna(seqs_l3, meta, hap_4char == "L3e1" | hap_2char == "REF")
final_l3e2 <- subset_dna(seqs_l3, meta, hap_4char == "L3e2" | hap_2char == "REF")
final_l3e3b <- subset_dna(seqs_l3, meta, hap_5char == "L3e3b" | hap_2char == "REF")
final_l3f1b <- subset_dna(seqs_l3, meta, hap_5char == "L3f1b" | hap_2char == "REF")
final_l3f1 <- subset_dna(seqs_l3, meta, hap_4char == "L3f1" | hap_2char == "REF")

# 5. write fasta files for iqtree
# A2
write.dna(final_a2, file = "iq_a2.fasta", format = "fasta")
write.dna(final_a264, file = "iq_a264.fasta", format = "fasta")

# B2
write.dna(final_b2, file = "iq_b2.fasta", format = "fasta")
write.dna(final_b2b, file = "iq_b2b.fasta", format = "fasta")

# C1
write.dna(final_c1, file = "iq_c1.fasta", format = "fasta")
write.dna(final_c1b, file = "iq_c1b.fasta", format = "fasta")

# D1
write.dna(final_d1, file = "iq_d1.fasta", format = "fasta")

# H1
write.dna(final_h1, file = "iq_h1.fasta", format = "fasta")
write.dna(final_h1bw, file = "iq_h1bw.fasta", format = "fasta")

# L0
write.dna(final_l0, file = "iq_l0.fasta", format = "fasta")
write.dna(final_l0a2a, file = "iq_l0a2a.fasta", format = "fasta")
write.dna(final_l0d1a, file = "iq_l0d1a.fasta", format = "fasta")

# L1
write.dna(final_l1, file = "iq_l1.fasta", format = "fasta")
write.dna(final_l1c1, file = "iq_l1c1.fasta", format = "fasta")
write.dna(final_l1c1a, file = "iq_l1c1a.fasta", format = "fasta")
write.dna(final_l1c1b, file = "iq_l1c1b.fasta", format = "fasta")
write.dna(final_l1c1d, file = "iq_l1c1d.fasta", format = "fasta")
write.dna(final_l1c2, file = "iq_l1c2.fasta", format = "fasta")
write.dna(final_l1c4, file = "iq_l1c4.fasta", format = "fasta")
write.dna(final_l1c24, file = "iq_l1c24.fasta", format = "fasta")

# L2
write.dna(final_l2, file = "iq_l2.fasta", format = "fasta")
write.dna(final_l2a1a, file = "iq_l2a1a.fasta", format = "fasta")
write.dna(final_l2a1d, file = "iq_l2a1d.fasta", format = "fasta")
write.dna(final_l2a1f, file = "iq_l2a1f.fasta", format = "fasta")
write.dna(final_l2a1q, file = "iq_l2a1q.fasta", format = "fasta")
write.dna(final_l2b1, file = "iq_l2b1.fasta", format = "fasta")
write.dna(final_l2c, file = "iq_l2c.fasta", format = "fasta")

# L3
write.dna(final_l3, file = "iq_l3.fasta", format = "fasta")
write.dna(final_l3b1a, file = "iq_l3b1a.fasta", format = "fasta")
write.dna(final_l3d1, file = "iq_l3d1.fasta", format = "fasta")
write.dna(final_l3d3, file = "iq_l3d3.fasta", format = "fasta")
write.dna(final_l3d4, file = "iq_l3d4.fasta", format = "fasta")
write.dna(final_l3e1a1a, file = "iq_l3e1a1a.fasta", format = "fasta")
write.dna(final_l3e1, file = "iq_l3e1.fasta", format = "fasta")
write.dna(final_l3e2, file = "iq_l3e2.fasta", format = "fasta")
write.dna(final_l3e3b, file = "iq_l3e3.fasta", format = "fasta")
write.dna(final_l3f1b, file = "iq_l3f1b.fasta", format = "fasta")
write.dna(final_l3f1, file = "iq_l3f1.fasta", format = "fasta")

# create metadata for iqtree samples
tree_samples <- unique(c(
  rownames(final_a2),
  rownames(final_b2),
  rownames(final_c1),
  rownames(final_d1),
  rownames(final_h1),
  rownames(final_l0),
  rownames(final_l1),
  rownames(final_l2),
  rownames(final_l3),
  rownames(final_a264),
  rownames(final_b2b),
  rownames(final_c1b),
  rownames(final_h1bw),
  rownames(final_l0a2a),
  rownames(final_l0d1a),
  rownames(final_l1c1a),
  rownames(final_l1c1b),
  rownames(final_l1c1d),
  rownames(final_l1c2),
  rownames(final_l1c24),
  rownames(final_l1c4),
  rownames(final_l2a1a),
  rownames(final_l2a1d),
  rownames(final_l2a1f),
  rownames(final_l2a1q),
  rownames(final_l2b1),
  rownames(final_l2c),
  rownames(final_l3b1a),
  rownames(final_l3d1),
  rownames(final_l3d3),
  rownames(final_l3d4),
  rownames(final_l3e1),
  rownames(final_l3e1a1a),
  rownames(final_l3e2),
  rownames(final_l3e3b),
  rownames(final_l3f1b),
  rownames(final_l3f1)
))

# subset metadata to only those samples
meta_iqtree <- meta %>%
  filter(Sample %in% tree_samples)

# write metadata
write.table(
  meta_iqtree,
  file = "iqtree.meta",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```



splitting for mjn:

```
# split_mjn.r
# load libraries
if (!require(pegas)) {
  install.packages("pegas")
  library(pegas)
}

if (!require(ape)) {
  install.packages("ape")
  library(ape)
}

if (!require(Biostrings)) {
  install.packages("Biostrings")
  library(Biostrings)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

# set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn")

# 1. set paths
a2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_A2.fasta")
b2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_B2.fasta")
c1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_C1.fasta")
d1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_D1.fasta")
h1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_H1.fasta")
l0_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L0.fasta")
l1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L1.fasta")
l2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L2.fasta")
l3_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L3.fasta")

# 2. read fastas as dnabin objects
read_in <- function(path) {
  fasta <- read.dna(path, format = "fasta")
  return(fasta)
}

seqs_a2 <- read_in(a2_path)
seqs_b2 <- read_in(b2_path)
seqs_c1 <- read_in(c1_path)
seqs_d1 <- read_in(d1_path)
seqs_h1 <- read_in(h1_path)
seqs_l0 <- read_in(l0_path)
seqs_l1 <- read_in(l1_path)
seqs_l2 <- read_in(l2_path)
seqs_l3 <- read_in(l3_path)

# 3. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/final_cut.meta", sep="\t")

# 4. create objects for each tree subset
subset_dna <- function(dna_obj, meta_df, filter_expr) {
  # Evaluate the filter expression within meta_df
  rows <- meta_df %>% filter(!!enquo(filter_expr))
  # Find row indices in DNAbin
  idx <- match(rows$Sample, rownames(dna_obj))
  # Warn if any samples are missing
  if(any(is.na(idx))) {
    warning("sequences in metadata but missing in tree: ", paste(rows$Sample[is.na(idx)], collapse = ", "))
  }
  # Subset and return DNAbin (only existing sequences)
  dna_obj[idx[!is.na(idx)], ]
}

# A2: A2+(64)
# crashes: final_a2 <- subset_dna(seqs_a2, meta, hap_2char == "A2")
final_a264 <- subset_dna(seqs_a2, meta, hap_7char == "A2+(64)")

# B2: B2, B2b
# crashes: final_b2 <- subset_dna(seqs_b2, meta, hap_2char == "B2")
final_b2b <- subset_dna(seqs_b2, meta, hap_3char == "B2b")

# C1: C1b
# crashes: final_c1 <- subset_dna(seqs_c1, meta, hap_2char == "C1")
final_c1b <- subset_dna(seqs_c1, meta, hap_3char == "C1b")

# D1: D1
final_d1 <- subset_dna(seqs_d1, meta, hap_2char == "D1")

# H1: H1, H1bw
# crashes: final_h1 <- subset_dna(seqs_h1, meta, hap_2char == "H1")
final_h1bw <- subset_dna(seqs_h1, meta, hap_4char == "H1bw")

# L0: L0a2a, L0d1a
# crashes: final_l0 <- subset_dna(seqs_l0, meta, hap_2char == "L0")
final_l0a2a <- subset_dna(seqs_l0, meta, hap_5char == "L0a2a")
final_l0d1a <- subset_dna(seqs_l0, meta, hap_5char == "L0d1a")

# L1: L1c, L1c1, L1c2'4
# final_l1 <- subset_dna(seqs_l1, meta, hap_2char == "L1")
final_l1c <- subset_dna(seqs_l1, meta, hap_3char == "L1c")
final_l1c1 <- subset_dna(seqs_l1, meta, hap_4char == "L1c1")
final_l1c24 <- subset_dna(seqs_l1, meta, hap_4char %in% c("L1c2", "L1c4"))

# L2: L2a1, L2a1a, L2a1d, L2a1f, L2a1q, L2b1, L2c
# final_l2      <- subset_dna(seqs_l2, meta, hap_2char == "L2")
# final_l2a1    <- subset_dna(seqs_l2, meta, hap_4char == "L2a1")
final_l2a1a   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1a")
final_l2a1d   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1d")
final_l2a1f   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1f")
final_l2a1q   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1q")
final_l2b1    <- subset_dna(seqs_l2, meta, hap_4char == "L2b1")
final_l2c     <- subset_dna(seqs_l2, meta, hap_3char == "L2c")

# L3: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b
# final_l3 <- subset_dna(seqs_l3, meta, hap_2char == "L3")
final_l3b1a <- subset_dna(seqs_l3, meta, hap_5char == "L3b1a")
final_l3d1 <- subset_dna(seqs_l3, meta, hap_4char == "L3d1")
final_l3d3 <- subset_dna(seqs_l3, meta, hap_4char == "L3d3")
final_l3d4 <- subset_dna(seqs_l3, meta, hap_4char == "L3d4")
final_l3e1a1a <- subset_dna(seqs_l3, meta, hap_7char == "L3e1a1a")
final_l3e2 <- subset_dna(seqs_l3, meta, hap_4char == "L3e2")
final_l3e3 <- subset_dna(seqs_l3, meta, hap_4char == "L3e3")
final_l3f1b <- subset_dna(seqs_l3, meta, hap_5char == "L3f1b")

# 5. write fasta files for iqtree
# A2
write.dna(final_a264, file = "mjn_a2_64.fasta", format = "fasta")

# B2
write.dna(final_b2b, file = "mjn_b2b.fasta", format = "fasta")

# C1
write.dna(final_c1b, file = "mjn_c1b.fasta", format = "fasta")

# D1
write.dna(final_d1, file = "mjn_d1.fasta", format = "fasta")

# H1
write.dna(final_h1bw, file = "mjn_h1bw.fasta", format = "fasta")

# L0
write.dna(final_l0a2a, file = "mjn_l0a2a.fasta", format = "fasta")
write.dna(final_l0d1a, file = "mjn_l0d1a.fasta", format = "fasta")

# L1
write.dna(final_l1c1, file = "mjn_l1c1.fasta", format = "fasta")
write.dna(final_l1c24, file = "mjn_l1c24.fasta", format = "fasta")

# L2
write.dna(final_l2a1a, file = "mjn_l2a1a.fasta", format = "fasta")
write.dna(final_l2a1d, file = "mjn_l2a1d.fasta", format = "fasta")
write.dna(final_l2a1f, file = "mjn_l2a1f.fasta", format = "fasta")
write.dna(final_l2a1q, file = "mjn_l2a1q.fasta", format = "fasta")
write.dna(final_l2b1, file = "mjn_l2b1.fasta", format = "fasta")
write.dna(final_l2c, file = "mjn_l2c.fasta", format = "fasta")

# L3
write.dna(final_l3b1a, file = "mjn_l3b1a.fasta", format = "fasta")
write.dna(final_l3d1, file = "mjn_l3d1.fasta", format = "fasta")
write.dna(final_l3d3, file = "mjn_l3d3.fasta", format = "fasta")
write.dna(final_l3d4, file = "mjn_l3d4.fasta", format = "fasta")
write.dna(final_l3e1a1a, file = "mjn_l3e1a1a.fasta", format = "fasta")
write.dna(final_l3e2, file = "mjn_l3e2.fasta", format = "fasta")
write.dna(final_l3e3, file = "mjn_l3e3.fasta", format = "fasta")
write.dna(final_l3f1b, file = "mjn_l3f1b.fasta", format = "fasta")

# create metadata for iqtree samples
tree_samples <- unique(c(
  rownames(final_a264),
  rownames(final_b2b),
  rownames(final_c1b),
  rownames(final_d1),
  rownames(final_h1bw),
  rownames(final_l0a2a),
  rownames(final_l0d1a),
  rownames(final_l1c1),
  rownames(final_l1c24),
  rownames(final_l2a1a),
  rownames(final_l2a1d),
  rownames(final_l2a1f),
  rownames(final_l2a1q),
  rownames(final_l2b1),
  rownames(final_l2c),
  rownames(final_l3b1a),
  rownames(final_l3d1),
  rownames(final_l3d3),
  rownames(final_l3d4),
  rownames(final_l3e1a1a),
  rownames(final_l3e2),
  rownames(final_l3e3),
  rownames(final_l3f1b)
))

# subset metadata to only those samples
meta_iqtree <- meta %>%
  filter(Sample %in% tree_samples)

# write metadata
write.table(
  meta_iqtree,
  file = "mjn.meta",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```



# sample size maps

iqtree:

```
# iqtree sample sizes

# load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)
library(stringr)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree/")

# read into data
meta <- read.csv("iqtree.meta", sep = "\t")

# clean country names & standardize
sort(unique(meta$country))

# remove empty ones
meta <- meta[!(meta$country %in% c("unknown", "")), ]

# remove parentheses
meta <- meta %>%
  mutate(country_clean = str_trim(gsub("\\s*\\([^\\)]*\\)", "", country)))

# fix 3 code 1kg & weird ones
meta <- meta %>%
  mutate(country_new = recode(country_clean,
                              "ACB" = "Barbados",
                              "ASW" = "USA",
                              "BEB" = "Bangladesh",
                              "CDX" = "China",
                              "CEU" = "USA",
                              "CHB" = "China",
                              "CHS" = "China",
                              "CLM" = "Colombia",
                              "ESN" = "Nigeria",
                              "FIN" = "Finland",
                              "GBR" = "UK",
                              "GIH" = "USA",
                              "GWD" = "Gambia",
                              "IBS" = "Spain",
                              "ITU" = "UK",
                              "JPT" = "Japan",
                              "KHV" = "Vietnam",
                              "LWK" = "Kenya",
                              "MSL" = "Sierra Leone",
                              "MXL" = "USA",
                              "PEL" = "Peru",
                              "PJL" = "Pakistan",
                              "PUR" = "Puerto Rico",
                              "STU" = "UK",
                              "TSI" = "Italy",
                              "YRI" = "Nigeria",
                              "Czechia" = "Czech Republic",
                              "United States" = "USA",
                              "Colombian" = "Colombia",
                              "Cambodian" = "Cambodia",
                              "PapuaNewGuinea" = "Papua New Guinea",     
                              "BotswanaOrNamibia" = "Namibia",
                              "SierraLeone" = "Sierra Leone",
                              "Abkhazia" = "Georgia",
                              "Congo" = "Democratic Republic of the Congo",              
                              "SouthAfrica" = "South Africa",         
                              "Kyrgyzystan" = "Kyrgyzstan",
                              "Korea" = "South Korea",
                              "Czechoslovia" = "Czech Republic",    
                              "Guinea Equatorial" = "Equatorial Guinea",
                              "Tenerife" = "Canary Islands",
                              "La Gomera" = "Canary Islands",
                              "Lanzarote" = "Canary Islands",
                              "Oromo" = "Ethiopia",
                              "São Tomé and Príncipe" = "Sao Tome and Principe",
                              "Sao Tome e Principe" = "Sao Tome and Principe",
                              "St. Helena" = "Saint Helena",
                              "St. Martin" = "Saint Martin",
                              "Scotland" = "UK",
                              "Great Britain" = "UK",
                              "United Kingdom" = "UK",
                              "England" = "UK",
                              "OrkneyIslands" = "UK",
                              "MakraniHGDP" = "Pakistan",        
                              "MbutiHGDP" = "Democratic Republic of the Congo",
                              "BiakaHGDP" = "Central African Republic",
                              "FrenchHGDP" = "France",       
                              "ColombianHGDP" = "Colombia",     
                              "SuruiHGDP" = "Brazil",
                              "MandenkaHGDP" = "Senegal",
                              "YorubaHGDP" = "Nigeria",    
                              "SanHGDP" = "Namibia",
                              "BantuSouthAfricaHGDP" = "South Africa",
                              "KaritianaHGDP" = "Brazil",
                              "BergamoItalianHGDP" = "Italy",
                              "DaurHGDP" = "China",
                              "BantuKenyaHGDP" = "Kenya",
                              "MozabiteHGDP" = "Algeria",
                              "BantuTswanaSGDP" = "Botswana",
                              "BantuHereroSGDP" = "Namibia",
                              "BantuKenyaSGDP" = "Kenya",
                              "BasqueSGDP" = "France",
                              "BiakaSGDP" = "Central African Republic",        
                              "BrahuiSGDP" = "Pakistan",  
                              "ChaneSGDP" = "Argentina",     
                              "ChukchiSGDP" = "Russia",    
                              "DinkaSGDP" = "Sudan",        
                              "EsanSGDP" = "Nigeria",        
                              "EskimoChaplinSGDP" = "Russia",   
                              "EskimoNaukanSGDP" = "Russia",  
                              "EskimoSirenikiSGDP" = "Russia", 
                              "EstonianSGDP" = "Estonia",      
                              "FrenchSGDP" = "France",         
                              "GambianSGDP" = "Gambia",       
                              "JordanianSGDP" = "Jordan",   
                              "JuhoanNorthSGDP" = "Namibia",     
                              "KaritianaSGDP" = "Brazil",      
                              "KhomaniSanSGDP" = "South Africa",    
                              "LuhyaSGDP" = "Kenya",         
                              "LuoSGDP" = "Kenya",             
                              "MakraniSGDP" = "Pakistan",      
                              "MansiSGDP" = "Russia",     
                              "MandenkaSGDP" = "Senegal",   
                              "MasaiSGDP" = "Kenya",        
                              "MayanSGDP" = "Mexico",       
                              "MbutiSGDP" = "Democratic Republic of the Congo",       
                              "MendeSGDP" = "Sierra Leone",        
                              "MixeSGDP" = "Mexico",       
                              "MixtecSGDP" = "Mexico",          
                              "MozabiteSGDP" = "Algeria",       
                              "PiapocoSGDP" = "Colombia",      
                              "PimaSGDP" = "Mexico",         
                              "QuechuaSGDP" = "Peru",       
                              "RussianSGDP" = "Russia",      
                              "SaharawiSGDP" = "Western Sahara",    
                              "SardinianSGDP" = "Italy",  
                              "SpanishSGDP" = "Spain",         
                              "SuruiSGDP" = "Brazil",           
                              "YorubaSGDP" = "Nigeria",      
                              "ZapotecSGDP" = "Mexico", 
                              .default = country_clean))  # keep any names not in codes

# compare between map data and mine
# base world map
map_data <- map_data("world")
setdiff(meta$country_new, map_data$region)

# grep("Congo", map_data$region, value = TRUE)

# make new dataframe for sample size by country
country_counts <- meta %>%
  group_by(country_new) %>%
  summarise(sample_size = n())

# get world map
ggplot(map_data, aes(x = long, y = lat, group = group, fill = region)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "none")

# merge my dataset w map data
merged_map_data <- map_data %>%
  left_join(country_counts, by = c("region" = "country_new"))

full <- ggplot(merged_map_data, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1)  +
  ggtitle("sample sizes of all used in trees")

# function to plot
plot_haplotype_map <- function(meta,
                               map_data,
                               hap_col,
                               hap_value,
                               title_prefix = "iqtree sample sizes by country") {
  
  hap_data <- meta %>%
    filter(grepl(hap_value, .data[[hap_col]]))
  
  hap_counts <- hap_data %>%
    group_by(country_new) %>%
    summarise(sample_size = n(), .groups = "drop")
  
  merged_data <- map_data %>%
    left_join(hap_counts, by = c("region" = "country_new"))
  
  plot <- ggplot(merged_data,
                 aes(x = long, y = lat, group = group, fill = sample_size)) +
    geom_polygon(color = "black", linewidth = 0.1) +
    coord_fixed(1.1) +
    theme_void() +
    theme(legend.position = "right") +
    scale_fill_viridis_c(option = "B", direction = -1) +
    ggtitle(paste(paste(hap_value, collapse = ", "), title_prefix))
  
  list(
    plot = plot,
    map_data = merged_data,
    meta_data = hap_data
  )
}


# a2+(64)
a2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value = "^A2([A-Za-z]?)$"
)

a264 <- plot_haplotype_map(
  meta = meta,
  map_data = map_data,
  hap_col = "hap_7char",
  hap_value = "A2\\+\\(64\\)"
)

# plot B2b countries
b2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "B2"
)

b2b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "B2b"
)

# plot C1b countries
c1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "C1"
)

c1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "C1b"
)

# plot D1 countries
d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "D1"
)

# plot H1bw countries
h1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "H1"
)

h1bw <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "H1bw"
)

# plot L0 countries
l0 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "L0"
)

l0a2a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0a2a"
)

l0d1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0d1a"
)

# plot L1 countries
l1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "L1"
)

l1c1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L1c1"
)

l1c24 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = c("L1c2", "L1c4")
)

# plot L2 countries
l2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "L2"
)

l2a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1a"
)

l2a1d <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1d"
)

l2a1f <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1f"
)

l2a1q <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1q"
)

l2b1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L2b1"
)

l2c <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "L2c"
)

# plot L3 countries
l3 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "L3"
)

l3b1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3b1a"
)

l3d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d1"
)

l3d3 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d3"
)

l3d4 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d4"
)

l3e1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e1"
)

l3e1a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_7char",
  hap_value  = "L3e1a1a"
)

l3e2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e2"
)

l3e3b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3e3b"
)

l3f1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3f1b"
)

# save
library(gridExtra)

pdf("maps_iqtree.pdf", width = 10, height = 8)
grid.arrange(full, ncol = 1)
grid.arrange(a2$plot, a264$plot, ncol = 1)
grid.arrange(b2$plot, b2b$plot, ncol = 1)
grid.arrange(c1$plot, c1b$plot, ncol = 1)
grid.arrange(d1$plot, h1$plot, ncol = 1)
grid.arrange(h1bw$plot, l0$plot, ncol = 1)
grid.arrange(l0a2a$plot, l0d1a$plot, ncol = 1)
grid.arrange(l1$plot, l1c1$plot, ncol = 1)
grid.arrange(l1c24$plot, l2$plot, ncol = 1)
grid.arrange(l2a1a$plot, l2a1d$plot, ncol = 1)
grid.arrange(l2a1f$plot, l2a1q$plot, ncol = 1)
grid.arrange(l2b1$plot,l2c$plot, ncol = 1)
grid.arrange(l3$plot, l3b1a$plot, ncol = 1)
grid.arrange(l3d1$plot, l3d3$plot, ncol = 1)
grid.arrange(l3d4$plot, l3e1$plot, ncol = 1)
grid.arrange(l3e1a1a$plot, l3e2$plot, ncol = 1)
grid.arrange(l3e3b$plot, l3f1b$plot, ncol = 1)
dev.off()

# also save metadata summary as text file
meta_list <- list(
  a2 = a2$meta_data,
  a264 = a264$meta_data,
  b2 = b2$meta_data,
  b2b = b2b$meta_data,
  c1 = c1$meta_data,
  c1b = c1b$meta_data,
  d1 = d1$meta_data,
  h1 = h1$meta_data,
  h1bw = h1bw$meta_data,
  l0 = l0$meta_data,
  l0a2a = l0a2a$meta_data,
  l0d1a = l0d1a$meta_data,
  l1 = l1$meta_data,
  l1c1 = l1c1$meta_data,
  l1c24 = l1c24$meta_data,
  l2 = l2$meta_data,
  l2a1a = l2a1a$meta_data,
  l2a1d = l2a1d$meta_data,
  l2a1f = l2a1f$meta_data,
  l2a1q = l2a1q$meta_data,
  l2b1  = l2b1$meta_data,
  l2c = l2c$meta_data,
  l3 = l3$meta_data,
  l3b1a = l3b1a$meta_data,
  l3d1  = l3d1$meta_data,
  l3d3 = l3d3$meta_data,
  l3d4 = l3d4$meta_data,
  l3e1 = l3e1$meta_data,
  l3e1a1a = l3e1a1a$meta_data,
  l3e2 = l3e2$meta_data,
  l3e3b = l3e3b$meta_data,
  l3f1b = l3f1b$meta_data
)

for (name in names(meta_list)) {
  paper <- meta_list[[name]] %>%
    count(paper, sort = TRUE)
  country <- meta_list[[name]] %>%
    count(country, sort = TRUE)
  extra_info <- meta_list[[name]] %>%
    count(extra_info, sort = TRUE)
  hap_summ <- meta_list[[name]] %>%
    count(hap_3char, sort = TRUE)
  hap <- meta_list[[name]] %>%
    count(Haplogroup, sort = TRUE)
  filename <- paste0(name, "_summ.txt")
  
  capture.output({
    cat("paper\n")
    print(paper)
    cat("country\n")
    print(country)
    cat("extra info from metadata\n")
    print(extra_info)
    cat("hap_3char\n")
    print(hap_summ)
    cat("haplogroup\n")
    print(hap)
  }, file = filename)
}




```



mjn:

```
# mjn sample sizes

# load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)
library(stringr)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn")

# read into data
meta <- read.csv("mjn.meta", sep = "\t")

# clean country names & standardize
sort(unique(meta$country))

# remove empty ones
meta <- meta[!(meta$country %in% c("unknown", "")), ]

# remove parentheses
meta <- meta %>%
  mutate(country_clean = str_trim(gsub("\\s*\\([^\\)]*\\)", "", country)))

# fix 3 code 1kg & weird ones
meta <- meta %>%
  mutate(country_new = recode(country_clean,
                              "ACB" = "Barbados",
                              "ASW" = "USA",
                              "BEB" = "Bangladesh",
                              "CDX" = "China",
                              "CEU" = "USA",
                              "CHB" = "China",
                              "CHS" = "China",
                              "CLM" = "Colombia",
                              "ESN" = "Nigeria",
                              "FIN" = "Finland",
                              "GBR" = "UK",
                              "GIH" = "USA",
                              "GWD" = "Gambia",
                              "IBS" = "Spain",
                              "ITU" = "UK",
                              "JPT" = "Japan",
                              "KHV" = "Vietnam",
                              "LWK" = "Kenya",
                              "MSL" = "Sierra Leone",
                              "MXL" = "USA",
                              "PEL" = "Peru",
                              "PJL" = "Pakistan",
                              "PUR" = "Puerto Rico",
                              "STU" = "UK",
                              "TSI" = "Italy",
                              "YRI" = "Nigeria",
                              "Czechia" = "Czech Republic",
                              "United States" = "USA",
                              "Colombian" = "Colombia",
                              "Cambodian" = "Cambodia",
                              "PapuaNewGuinea" = "Papua New Guinea",     
                              "BotswanaOrNamibia" = "Namibia",
                              "SierraLeone" = "Sierra Leone",
                              "Abkhazia" = "Georgia",
                              "Congo" = "Democratic Republic of the Congo",              
                              "SouthAfrica" = "South Africa",         
                              "Kyrgyzystan" = "Kyrgyzstan",
                              "Korea" = "South Korea",
                              "Czechoslovia" = "Czech Republic",    
                              "Guinea Equatorial" = "Equatorial Guinea",
                              "Tenerife" = "Canary Islands",
                              "La Gomera" = "Canary Islands",
                              "Lanzarote" = "Canary Islands",
                              "Oromo" = "Ethiopia",
                              "São Tomé and Príncipe" = "Sao Tome and Principe",
                              "Sao Tome e Principe" = "Sao Tome and Principe",
                              "St. Helena" = "Saint Helena",
                              "St. Martin" = "Saint Martin",
                              "Scotland" = "UK",
                              "Great Britain" = "UK",
                              "United Kingdom" = "UK",
                              "England" = "UK",
                              "OrkneyIslands" = "UK",
                              "MakraniHGDP" = "Pakistan",        
                              "MbutiHGDP" = "Democratic Republic of the Congo",
                              "BiakaHGDP" = "Central African Republic",
                              "FrenchHGDP" = "France",       
                              "ColombianHGDP" = "Colombia",     
                              "SuruiHGDP" = "Brazil",
                              "MandenkaHGDP" = "Senegal",
                              "YorubaHGDP" = "Nigeria",    
                              "SanHGDP" = "Namibia",
                              "BantuSouthAfricaHGDP" = "South Africa",
                              "KaritianaHGDP" = "Brazil",
                              "BergamoItalianHGDP" = "Italy",
                              "DaurHGDP" = "China",
                              "BantuKenyaHGDP" = "Kenya",
                              "MozabiteHGDP" = "Algeria",
                              "BantuTswanaSGDP" = "Botswana",
                              "BantuHereroSGDP" = "Namibia",
                              "BantuKenyaSGDP" = "Kenya",
                              "BasqueSGDP" = "France",
                              "BiakaSGDP" = "Central African Republic",        
                              "BrahuiSGDP" = "Pakistan",  
                              "ChaneSGDP" = "Argentina",     
                              "ChukchiSGDP" = "Russia",    
                              "DinkaSGDP" = "Sudan",        
                              "EsanSGDP" = "Nigeria",        
                              "EskimoChaplinSGDP" = "Russia",   
                              "EskimoNaukanSGDP" = "Russia",  
                              "EskimoSirenikiSGDP" = "Russia", 
                              "EstonianSGDP" = "Estonia",      
                              "FrenchSGDP" = "France",         
                              "GambianSGDP" = "Gambia",       
                              "JordanianSGDP" = "Jordan",   
                              "JuhoanNorthSGDP" = "Namibia",     
                              "KaritianaSGDP" = "Brazil",      
                              "KhomaniSanSGDP" = "South Africa",    
                              "LuhyaSGDP" = "Kenya",         
                              "LuoSGDP" = "Kenya",             
                              "MakraniSGDP" = "Pakistan",      
                              "MansiSGDP" = "Russia",     
                              "MandenkaSGDP" = "Senegal",   
                              "MasaiSGDP" = "Kenya",        
                              "MayanSGDP" = "Mexico",       
                              "MbutiSGDP" = "Democratic Republic of the Congo",       
                              "MendeSGDP" = "Sierra Leone",        
                              "MixeSGDP" = "Mexico",       
                              "MixtecSGDP" = "Mexico",          
                              "MozabiteSGDP" = "Algeria",       
                              "PiapocoSGDP" = "Colombia",      
                              "PimaSGDP" = "Mexico",         
                              "QuechuaSGDP" = "Peru",       
                              "RussianSGDP" = "Russia",      
                              "SaharawiSGDP" = "Western Sahara",    
                              "SardinianSGDP" = "Italy",  
                              "SpanishSGDP" = "Spain",         
                              "SuruiSGDP" = "Brazil",           
                              "YorubaSGDP" = "Nigeria",      
                              "ZapotecSGDP" = "Mexico", 
                              .default = country_clean))  # keep any names not in codes

# compare between map data and mine
# base world map
map_data <- map_data("world")
setdiff(meta$country_new, map_data$region)

# grep("Congo", map_data$region, value = TRUE)

# make new dataframe for sample size by country
country_counts <- meta %>%
  group_by(country_new) %>%
  summarise(sample_size = n())

# get world map
ggplot(map_data, aes(x = long, y = lat, group = group, fill = region)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "none")

# merge my dataset w map data
merged_map_data <- map_data %>%
  left_join(country_counts, by = c("region" = "country_new"))

full <- ggplot(merged_map_data, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1)  +
  ggtitle("sample sizes of all used in networks")

# function to plot
plot_haplotype_map <- function(meta,
                               map_data,
                               hap_col,
                               hap_value,
                               title_prefix = "mjn sample sizes by country") {
  
  hap_data <- meta %>%
    filter(grepl(hap_value, .data[[hap_col]]))
  
  hap_counts <- hap_data %>%
    group_by(country_new) %>%
    summarise(sample_size = n(), .groups = "drop")
  
  merged_data <- map_data %>%
    left_join(hap_counts, by = c("region" = "country_new"))
  
  plot <- ggplot(merged_data,
                 aes(x = long, y = lat, group = group, fill = sample_size)) +
    geom_polygon(color = "black", linewidth = 0.1) +
    coord_fixed(1.1) +
    theme_void() +
    theme(legend.position = "right") +
    scale_fill_viridis_c(option = "B", direction = -1) +
    ggtitle(paste(paste(hap_value, collapse = ", "), title_prefix))
  
  list(
    plot = plot,
    map_data = merged_data,
    meta_data = hap_data
  )
}


# plot A2+(64) countries
a264 <- plot_haplotype_map(
  meta = meta,
  map_data = map_data,
  hap_col = "hap_7char",
  hap_value = "A2\\+\\(64\\)"
)

# plot B2b countries
b2b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "B2b"
)

# plot C1b countries
c1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "C1b"
)

# plot D1 countries
d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "D1"
)

# plot H1bw countries
h1bw <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "H1bw"
)


# plot L0 countries
#l0a2a
l0a2a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0a2a"
)

l0d1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0d1a"
)

# plot L1 countries
l1c1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L1c1"
)

l1c24 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = c("L1c2", "L1c4")
)

# plot L2 countries
l2a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1a"
)

l2a1d <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1d"
)

l2a1f <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1f"
)

l2a1q <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1q"
)

l2b1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L2b1"
)

l2c <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "L2c"
)

# plot L3 countries
l3b1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3b1a"
)

l3d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d1"
)

l3d3 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d3"
)

l3d4 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d4"
)

l3e1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e1"
)

l3e1a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_7char",
  hap_value  = "L3e1a1a"
)

l3e2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e2"
)

l3e3b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3e3b"
)

l3f1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3f1b"
)

# save
library(gridExtra)

pdf("maps_mjn.pdf", width = 10, height = 8)
grid.arrange(full, ncol = 1)
grid.arrange(a264$plot, b2b$plot, ncol = 1)
grid.arrange(c1b$plot, d1$plot, ncol = 1)
grid.arrange(h1bw$plot, l0a2a$plot, ncol = 1)
grid.arrange(l0d1a$plot, l1c1$plot, ncol = 1)
grid.arrange(l1c24$plot, l2a1a$plot, ncol = 1)
grid.arrange(l2a1d$plot, l2a1f$plot, ncol = 1)
grid.arrange(l2a1q$plot, l2b1$plot, ncol = 1)
grid.arrange(l2c$plot, l3b1a$plot, ncol = 1)
grid.arrange(l3d1$plot, l3d3$plot, ncol = 1)
grid.arrange(l3d4$plot, l3e1a1a$plot, ncol = 1)
grid.arrange(l3e2$plot, l3e3b$plot, ncol = 1)
grid.arrange(l3f1b$plot, ncol = 1)
dev.off()

# also save metadata summary as text file
meta_list <- list(
  a264 = a264$meta_data,
  b2b = b2b$meta_data,
  c1b = c1b$meta_data,
  d1 = d1$meta_data,
  h1bw = h1bw$meta_data,
  l0a2a = l0a2a$meta_data,
  l0d1a = l0d1a$meta_data,
  l1c1 = l1c1$meta_data,
  l1c24 = l1c24$meta_data,
  l2a1a = l2a1a$meta_data,
  l2a1d = l2a1d$meta_data,
  l2a1f = l2a1f$meta_data,
  l2a1q = l2a1q$meta_data,
  l2b1  = l2b1$meta_data,
  l2c = l2c$meta_data,
  l3b1a = l3b1a$meta_data,
  l3d1  = l3d1$meta_data,
  l3d3 = l3d3$meta_data,
  l3d4 = l3d4$meta_data,
  l3e1 = l3e1$meta_data,
  l3e1a1a = l3e1a1a$meta_data,
  l3e2 = l3e2$meta_data,
  l3e3b = l3e3b$meta_data,
  l3f1b = l3f1b$meta_data
)

for (name in names(meta_list)) {
  paper <- meta_list[[name]] %>%
    count(paper, sort = TRUE)
  country <- meta_list[[name]] %>%
    count(country, sort = TRUE)
  extra_info <- meta_list[[name]] %>%
    count(extra_info, sort = TRUE)
  hap_summ <- meta_list[[name]] %>%
    count(hap_3char, sort = TRUE)
  hap <- meta_list[[name]] %>%
    count(Haplogroup, sort = TRUE)
  filename <- paste0(name, "_summ.txt")
  
  capture.output({
    cat("paper\n")
    print(paper)
    cat("country\n")
    print(country)
    cat("extra info from metadata\n")
    print(extra_info)
    cat("hap_3char\n")
    print(hap_summ)
    cat("haplogroup\n")
    print(hap)
  }, file = filename)
}

```

