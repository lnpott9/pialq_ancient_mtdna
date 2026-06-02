# 05_log2fasta

*$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/*

*finished 7/8/25*

Summary of schmutzi results:

| 1. Uncontaminated after initial schmutzi run with no prediction (npred) | 2. Uncontaminated after schmutzi run with prediction (wpred) | 3. Decontaminated using PMDtools | 4. Uncontaminated after step 3 and schmutzi run (npred) | 5. Uncontaminated after step 3 and schmutzi (wpred) | Excluded due to low coverage/would not run through schmutzi |
| ------------------------------------------------------------ | ------------------------------------------------------------ | -------------------------------- | ------------------------------------------------------- | --------------------------------------------------- | ----------------------------------------------------------- |
| 107mt<br/>108mt<br/>109mt<br/>110mt<br/>111mt<br/>112mt<br/>113mt<br/>114mt<br/>115mt<br/>116mt<br/>117mt<br/>118mt<br/>119mt<br/>120mt<br/>121mt<br/>122mt<br/>123mt<br/>124mt<br/>125mt<br/>126mt<br/>127mt<br/>128mt<br/>129mt<br/>130mt<br/>131mt<br/>132mt<br/>133mt<br/>134mt<br/>135mt<br/>136mt<br/>33dsV2<br/>53mt<br/>54mt<br/>55mt<br/>59mt<br/>61mt<br/>62mt<br/>63mt<br/>64mt<br/>65mt<br/>66mt<br/>67mt<br/>68mt<br/>69mt<br/>70mt<br/>71mt<br/>72mt<br/>73mt<br/>74mt<br/>75mt<br/>76mt<br/>77mt<br/>79mt<br/>80mt<br/>81mt<br/>82mt<br/>83mt<br/>84mt<br/>85mt<br/>LQ30a.1 | 19dsV2                                                       | 56mt<br/>78mt                    | 56mt_pmd3                                               | 78mt_pmd3                                           | 57mt<br/>58mt<br/>60mt<br/>106mt<br/>                       |

## 01. Consensus fasta file creation

*$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq*

The next step is to create consensus fasta files for the PIALQ individuals using the `*_schmutzi_npred_final_endo.log` files, which contain information on the quality of the consensus fasta sequence. 

File sources:

```
# 19dsV2: wpred
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/01_wpred/dsdna/uf_19dsV2_schmutzi
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/01_wpred/dsdna/f_19dsV2_schmutzi

# 56mt: decontamination with PMDtools followed by npred
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/03_npred
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/03_npred

# 78mt: decontamination with PMDtools followed by npred followed by wpred
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/04_wpred
$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/04_wpred

# filtered files (q20, primary only)
$SHARED/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/*sdna/*_schmutzi/*npred_final.cont.est 

# unfiltered files (q20)
$SHARED/projects/PIALQ/2025_ancient_lp/03_schmutzi/unfilt/*sdna/*_schmutzi/*npred_final.cont.est 
```

I copied all files to *$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus*, with filtered files in *$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/filt* and unfiltered files in *$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/unfilt* since they are still not differentiated by file name. I appended "f_" or "uf_" to the beginning of each files name using the following code:

```
# for unfiltered files
for file in *_schmutzi_npred_final_endo.log; do
    mv "$file" "uf_$file"
done

# for filtered files
for file in *_schmutzi_npred_final_endo.log; do
    mv "$file" "f_$file"
done
```

I then moved all log files to the same directory (`$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/01_log_files`) and created the consensus fastas using the following code, using a quality score threshold of 30 for individual sites and 20 for indels:

```
# start interactive session
srun --mem=4gb --time=01:00:00 --pty bash

# set paths to executable and output directory
log2fasta=/home/mnievesc/shared/programs/schmutzi/src/log2fasta
output=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq

# loop through files and call consensus
for log in *_schmutzi_npred_final_endo.log; do
	basename=$(basename "$log" _schmutzi_npred_final_endo.log)
	"$log2fasta" -q 30 -name "${basename}" -indel 20 "$log" > "$output/${basename}.fasta"
done

for log in *_schmutzi_wpred_final_endo.log; do
	basename=$(basename "$log" _schmutzi_wpred_final_endo.log)
	"$log2fasta" -q 30 -name "${basename}" -indel 20 "$log" > "$output/${basename}.fasta"
done
```

The final consensus files are here:

```
$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq
```

