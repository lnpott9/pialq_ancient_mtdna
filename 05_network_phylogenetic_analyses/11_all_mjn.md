# 11. mjns for haplogroups:

split mjn

in 09_msa

cd /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/03_all

## 06. prepare for popart

a. **make nexus files**

```
# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env
base="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn"

python $SHARED/programs/AMAS/amas/AMAS.py convert -i "$base"/mjn*fasta -f fasta -u nexus -d dna
```

b. **replace Ns with ?s**

```
for nexus in mjn*.fasta-out.nex; do
name=$(basename "$nexus" .fasta-out.nex)
python replace.py "$nexus" "rep_${name}.nex"
done
```



*replace.py*

```
#!/usr/bin/env python3
import argparse

def replace_N_with_missing(infile, outfile):
    with open(infile, "r") as f:
        lines = f.readlines()

    new_lines = []
    in_matrix = False

    for line in lines:
        stripped = line.strip()

        # Detect start/end of the MATRIX section
        if stripped.upper().startswith("MATRIX"):
            in_matrix = True
            new_lines.append(line)
            continue
        if in_matrix and stripped.upper().startswith(";"):
            in_matrix = False
            new_lines.append(line)
            continue

        if in_matrix and stripped:  # inside MATRIX
            parts = stripped.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                # Replace Ns (upper/lower) with ?
                seq = seq.replace("N", "?").replace("n", "?")
                new_line = f"{name}\t{seq}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write out modified nexus file
    with open(outfile, "w") as f:
        f.writelines(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Replace Ns in NEXUS DNA sequences with ?"
    )
    parser.add_argument("infile", help="Input NEXUS file")
    parser.add_argument("outfile", help="Output NEXUS file with Ns replaced by ?")

    args = parser.parse_args()
    replace_N_with_missing(args.infile, args.outfile)


if __name__ == "__main__":
    main()
```



c. **add traits block for paper**

```
meta="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta"

for nexus in rep_mjn*.nex; do
name=$(basename "$nexus" .nex)
python paper.py "$nexus" "$meta" popart/"${name}_paper.nex"
done
```

*paper.py*

```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = "paper"   # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")


# -----------------------------
# 2b. Clean sample names
# -----------------------------
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)


# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
with open(nexus_input_path, "r") as f:
    nexus_text = f.read()

match = re.search(r"MATRIX\s*(.*?);\s*", nexus_text, re.IGNORECASE | re.DOTALL)
if not match:
    raise ValueError("Could not find MATRIX block in the NEXUS file.")

matrix_block = match.group(1)

dna_samples = []
for line in matrix_block.strip().splitlines():
    parts = line.strip().split()
    if parts:
        dna_samples.append(parts[0])


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
trait_block_lines = []
trait_block_lines.append("BEGIN TRAITS;")
trait_block_lines.append(f"\tDimensions NTRAITS={len(trait_labels_clean)};")
trait_block_lines.append("\tFormat labels=yes missing=? separator=Comma;")
trait_block_lines.append("\tTraitLabels " + " ".join(trait_labels_clean) + ";")
trait_block_lines.append("\tMatrix")

for sample in binary_traits.index:
    values = ",".join(str(v) for v in binary_traits.loc[sample])
    trait_block_lines.append(f"\t{sample}\t{values}")

trait_block_lines.append(";")
trait_block_lines.append("END;")

traits_block = "\n".join(trait_block_lines)


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```



d. **add traits block of country**

```
meta="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta"

for nexus in rep_mjn*.nex; do
name=$(basename "$nexus" .nex)
python country.py "$nexus" "$meta" popart/"${name}_country.nex"
done
```



```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = "country_label"   # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")


# -----------------------------
# 2b. Clean sample names
# -----------------------------
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)


# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
with open(nexus_input_path, "r") as f:
    nexus_text = f.read()

match = re.search(r"MATRIX\s*(.*?);\s*", nexus_text, re.IGNORECASE | re.DOTALL)
if not match:
    raise ValueError("Could not find MATRIX block in the NEXUS file.")

matrix_block = match.group(1)

dna_samples = []
for line in matrix_block.strip().splitlines():
    parts = line.strip().split()
    if parts:
        dna_samples.append(parts[0])


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
trait_block_lines = []
trait_block_lines.append("BEGIN TRAITS;")
trait_block_lines.append(f"\tDimensions NTRAITS={len(trait_labels_clean)};")
trait_block_lines.append("\tFormat labels=yes missing=? separator=Comma;")
trait_block_lines.append("\tTraitLabels " + " ".join(trait_labels_clean) + ";")
trait_block_lines.append("\tMatrix")

for sample in binary_traits.index:
    values = ",".join(str(v) for v in binary_traits.loc[sample])
    trait_block_lines.append(f"\t{sample}\t{values}")

trait_block_lines.append(";")
trait_block_lines.append("END;")

traits_block = "\n".join(trait_block_lines)


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```



e. **add traits block of haplogroup**

```
meta="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta"

three=("rep_mjn_d1.nex")
four=("rep_mjn_l2c.nex")
five=("rep_mjn_b2b.nex" "rep_mjn_c1b.nex" "rep_mjn_h1bw.nex" "rep_mjn_l1c1.nex" "rep_mjn_l1c24.nex" "rep_mjn_l2b1.nex" "rep_mjn_l3d1.nex" "rep_mjn_l3d3.nex" "rep_mjn_l3d4.nex" "rep_mjn_l3e1.nex" "rep_mjn_l3e2.nex" "rep_mjn_l3e3.nex")
six=("rep_mjn_l0a2a.nex" "rep_mjn_l0d1a.nex" "rep_mjn_l2a1a.nex" "rep_mjn_l2a1d.nex" "rep_mjn_l2a1f.nex" "rep_mjn_l2a1q.nex" "rep_mjn_l3b1a.nex" "rep_mjn_l3e3b.nex" "rep_mjn_l3f1b.nex")
eight=("rep_mjn_a2_64.nex" "rep_mjn_l3e1a1a.nex")

for nexus in "${three[@]}"; do
name=$(basename "$nexus" .nex)
python hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_3char
done

for nexus in "${four[@]}"; do
name=$(basename "$nexus" .nex)
python hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_4char
done

for nexus in "${five[@]}"; do
name=$(basename "$nexus" .nex)
python hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_5char
done

for nexus in "${six[@]}"; do
name=$(basename "$nexus" .nex)
python hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_6char
done

for nexus in "${eight[@]}"; do
name=$(basename "$nexus" .nex)
python hap.py "$nexus" "$meta" popart/"${name}_hap.nex" Haplogroup
done
```



*hap.py*

```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = sys.argv[4] # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")


# -----------------------------
# 2b. Clean sample names
# -----------------------------
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)


# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
with open(nexus_input_path, "r") as f:
    nexus_text = f.read()

match = re.search(r"MATRIX\s*(.*?);\s*", nexus_text, re.IGNORECASE | re.DOTALL)
if not match:
    raise ValueError("Could not find MATRIX block in the NEXUS file.")

matrix_block = match.group(1)

dna_samples = []
for line in matrix_block.strip().splitlines():
    parts = line.strip().split()
    if parts:
        dna_samples.append(parts[0])


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
trait_block_lines = []
trait_block_lines.append("BEGIN TRAITS;")
trait_block_lines.append(f"\tDimensions NTRAITS={len(trait_labels_clean)};")
trait_block_lines.append("\tFormat labels=yes missing=? separator=Comma;")
trait_block_lines.append("\tTraitLabels " + " ".join(trait_labels_clean) + ";")
trait_block_lines.append("\tMatrix")

for sample in binary_traits.index:
    values = ",".join(str(v) for v in binary_traits.loc[sample])
    trait_block_lines.append(f"\t{sample}\t{values}")

trait_block_lines.append(";")
trait_block_lines.append("END;")

traits_block = "\n".join(trait_block_lines)


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```





notes on networks:

**A2**: LQ28, LQ61

```
A2: crashes

A2+(64):
- "starlike" pattern from a central node of real samples
- also no clear haplogroup clusters bc these are the A2 lineages with extra mutations
- LQ28 & LQ61 are designated as identical sequences
- cluster is connected directly to central node, includes individuals from tito_unpub, 1kg, barbieri 2017, nakatsuka_2020, all collected in Peru
```

**B2** (note: 1 B2 individual not included in B2b): LQ05, LQ49 + LQ14

```
B2: crashes

B2b:
- no clear haplogroup clusters because many sequences are just "B2b" w/ no further granularity
-LQ05 & LQ49 are not very close to each other in the network
- both are most mutationally similar to other B2b lineages from Peru
- "starlike" pattern from a reconstructed central node
- LQ49: 1kg, barbieri_2017, tito_unpub
- LQ05: tito_unpub, valverde_2016
```

**C1**: LQ06, LQ23, LQ60, LQ62

```
C1b:
- "starlike" pattern from a central node of real samples from Chile & Peru
- looks crazy
- LQ06, LQ23, and LQ62 are designated as identical sequences, directly connected to central node as well as 1 sample JX669195 from tito_unpub, Huancayo, Junín, Peru
- LQ60: more on the outskirts of the network, 7 mutations from the central node, connected to node it shares with KU523330 (llamas_2016), ancient sequence from Chancay
```

**D1**: LQ01, LQ46, LQ58

```
D1:
- "starlike" pattern from a central node of real samples
- LQ01 & LQ58 are designated as identical sequences, LQ46 is in a different branch
- LQ01 & LQ58 are in a cluster with only other Peruvian samples from huber_2025, arias_2017, barbieri_2017
- LQ46 is in a cluster with Peru & Ecuador samples from tito_unpub, huber_2025, barbieri_2017
```

**H**: LQ18

```
H1b: 
- some starlike patterns within bigger network
- mutationally close to individuals from Canary Islands, Portugal, Spain, Italy from garciaolivares_2023, silva_2021, 1kg, bodner_2022
- edge of network
```

**L0**: LQ26, LQ24, LQ27

```
L0a2a: LQ44, LQ26
- clear split into L0a2a1 and L0a2a2
- some starlike expansions
- LQ44 by itself, comes off of a node whose expansions are samples from Zambia, Kenya, US, Madagascar, Angola (L0a2a1) from barbieri_2014a, brucato_2018, just_2015, pierron_2017, barbieri_2012b
- LQ26 designated an identical sequence in node KC622056 along with 109 other sequences: majority from Madagascar Zambia, Kenya, South Africa, also sequences from Comoros, Angola, ASW, Zambia, unknown, US, Somalia, Mozambique (L0a2a2) from pierron_2017, mccrow_2016, chan_2019, brucato_2018, barbieri_2014b, rito_2013, gonder_2007, barbieri_2014a, 1kg, taylor_2020
- not to edges of network

L0d1a: LQ27 
- LQ27 connected to a starlike expansion node of several samples from Namibia, Mozambique, and South Africa form chan_2019, barbieri_2013, points include rito_2013, mccrow_2016
```

**L1**: LQ08, LQ15, LQ39, LQ41, LQ45, LQ48, LQ54

```
L1c:
- LQ45 & LQ48 are designated as identical sequences
- LQ08 by itself
- LQ39 by itself, place between L1c2 and L1c4 isn't clarified, still between clusters
- LQ15 designated as identical with 2 other sequences, KJ185518 & KJ185762, in node KJ185518
- LQ54 & LQ41 attached to the same node

L1c1: (no LQ39)
L1c1a: LQ08, LQ45, LQ48 (LQ45 & 48 identical, LQ08 off on its own)
L1c1b: LQ15 under node KJ185518
L1c1d: LQ54, LQ41 (close but not identical)
- LQ08: in a branch with ACB, United States, Paraguay, not particularly close to anything (1kg, barbieri_2012b, taylor_2020, simao_2019)
- LQ15: in node KJ185518, identical with sequences from Angola, connected to ASW & Zambia, barbieri_2014a & 1kg, oliveira_2018, barbieri_2014b
- LQ41: not super clear, connected to samples from ASW & US & Zambia, taylor_2020, 1kg, barbieri_2014a
- LQ54: connected to node that branches off to ACB & Peru, huber_2025 & 1kg
- LQ45/48: 1 mutational difference from a sample from Brazil, otherwise not very similar to anything else, from avila_2019

L1c2'4: LQ39
- still out in the middle
```

**L2**: LQ22, LQ29, LQ63, LQ59, 

```
L2a1: crashes

L2a1a: LQ22, LQ29
- to edges of network
- LQ22& LQ29 are identified as identicial sequences
- clear cluster separation into L2a1a1, 2, & 3
- in branch with Namibia & Angola samples (Barbieri_2014b, barbieri_2014a), also close to ESN, St Helena, Sao Tome & Principe, Madagascar, Nigeria, Yoruba HGDP & SGDP (1kg, hgdp, sgdp, sandovalvelasco_2023, pierron_2017, silva_2015, barbiero_2012a)
- reasonably close to STH_254 (Harney)

L2a1d: LQ63
- not really near anything, not very informative

L2a1f: LQ56
- in a branch coming off a central cluster
- other samples in the branch include STH_253, STH_215, STH_347, STH_441, STM3, and DQ304954 from locations St. Helena, St. Martin, papers (sandovalvelasco_2023, schroeder_2015, just_2008, just_2015)

L2a1q: LQ59
- only 5 sequences total, not very informative
- only sequences are form Brazil, LWK, Zambia, unknown from papers 1kh, avila_2019, barbieri 2012b, barbieri_2014b, unknown is Kwanyama (national language of Angola and Namibia)

L2b1: LQ02, LQ19, LQ16, LQ17, LQ12, LQ55, LQ31
L2b1a (star): LQ02 (outskirts of network, lots of mutational differences), LQ19 & LQ16 are designated identical sequences with 9 other sequences in the central node:
DQ304978
DQ304979
DQ304980
DQ304982
DQ304983
KJ185443
KJ185460
MF621123
NA19984
L2b1b (branches): LQ12, LQ17, LQ31, LQ55
- LQ02: connected to node from Angola by a million mutations, same node connects to St. Helena (sandovalvelasco_2012, barbieri_2014a)
- LQ16 & LQ19 (more of a diaspora node): node mostly contains US, followed by Peru, Nigeria, ASW, Zambia from papers (just_2008, 1kg, barbieri_2014a, cabrera_2018)
- L2b1b clump: connected to samples from Zambia, ESN, YRI, USA, Burkina Faso from papers (1kg, taylor_2020, barbieri_2014a, barbieri_2012a)

L2c: LQ52, LQ30, LQ03 (no sharing)
- LQ52: outskirts of network in L2c clump, no clear cluster, close to GWD, MandenkaHGDP, Burkin Faso, PUR (bring up), USA from papers (1kg, barbieri_2012a, silva_2021, hgdp)
- LQ30: outskirts of network in L2c2 clump from Zambia, ESN, Libya, USA form papers (1kg, barbieri_2014a, barbieri_2012b, colombo_2025, taylor_2020, just_2008)
- LQ03: L2c3, outskirts as well from ACB, Angola from papers (1kg, oliveira_2018, barbieri_2014a)
```

**L3**: LQ11, LQ38

```
L3b1a: LQ11, LQ38 (ugly network)
LQ38: 
- none are identical sequences
-L3b1a8, on outskirts of branch, connected to same node as 2 US other samples from taylor_2020, branch includes samples from pierron_2017, 1kg, madagascar, YRI
LQ11: 
-L3b1a1: giant starlike expansion pattern for several nodes, but LQ11 is in between all of them, connected by 2 miutatiosn to giant node with samples form pierron_2017, 1kg, just_2015, brucato_2018 (kenya, usa, comoros, madagascar)

L3d1: LQ26, LQ09, LQ24, LQ33
- clear cluster separation
- none are identical sequences
- LQ24, LQ33 separated by 11 mutations, connected to a node with samples from USA, Zambia, Angola, ASW, Canary islands, Spain from papers (barbieri_2014a, barbieri_2012b, 1kg, taylor_2020, garciaolivares_2023, cabrera-2018, hgdp)
- LQ26: not too close to anything, cameroon & us from papers ()
- LQ09 (bring up): connected to central node of STM2 from schroeder_2015 st martin, ACB, US, GWD, Burkina Faso, ALSO connected to a node w a PUR individual, IBS, Dominican, US from (schroeder_2015, 1kg, barbieri_2012a, taylor_2020, just_2015) --> node w pur individual is 1khg, torroni_2006, taylor_2020) --> labeled DQ341072

L3d3: LQ64 (bring up, seen in modern Peru today, very southern africa)
- no identical sequences
- LQ64: connected to a very large central node composed of samples from Namibia, Angola, BantuHereroSGDP, BantuSouthAfricaHGDP
- PQ827230 is also a sample form Peru connected to same large central node from huber_2025
- all spokes of star are either from Namibia, Angola, or Zambia
- papers: barbieri_2014b, oliveira_2018, barbieri_2014a, hgdp

L3d4: LQ40
- only 15 sequences, not so informative
- LQ40 connected to node with other branches ASW, US, Burkina Faso, Sudan, Tunisia
- costa_2009, soares_2012, 1kg, taylor_2020, barbieri_2012a

L3e1a1a: LQ13, LQ57, LQ10, LQ50 (maybe a network to show)
- LQ13, LQ50, LQ57 are identical sequences
- the identical sequences are 3 mutations away fro mlarge central node
- LQ10 is 1 mutation away from central node
- central node has samples from mostly madagascar, then angola, zambia, south africa, somalia, PERU
- other samples from US in plot
- central node is form papers barbieri_2014a, pierron_2017, mccrow_2016, soares_2012, barbieri_2017, peru is form barbieri_2017 and is from UtcubambaSouth

L3e2: LQ35, LQ42, LQ32, LQ36 (very 1kg heavy)
- clear cluster separation
- LQ42: outskirts of network, ttached to a node that's attached to the central L3e2b node, attached to sample from Burkina faso, central node has samples form US, Botswana, Burkina Faso, Zambia, MSL, Dominican republic from papers ()
- LQ32 & LQ36 are designated as identical sequences (node is called KC622172) - L3e2b: contains Zambia, Angola, ASW from papers (1kg, barbieri_2014b, barbieri_2014a, cabrera_2018)
- LQ35 (L3e2a, which has 3 major clusters): connected to YRI & Nigeria (hgdp, barbieri_2012a), larger node also is connected to  samples from GWD, Burkina Faso, ESN, ASW, LWK, USA from papers (1kg, taylor_2020, babieri_2012a), CHS36 is also in there

L3e3: LQ43 (diaspora node)
- clear cluster separation between l3e3a and l3e3b
- comes directly off central node of l3e3b by itself
- no identical with any others
- HG02006 from PEL is also connected to the central node
- central node has 1kg, ASW, Zambia, YRI, Namibia, Madagascar, Angola
- chad connected to network
- node: barbieri_2014a, 1kg, taylor_2020, pierron_2017, just_2015, just_2008, barbieri_2014b

L3f1b: LQ21, LQ34, LQ20
- clearish cluster separation
- LQ20 identical to 4 other sequences in node KJ185647, other sequences are from Brazil, Zambia, US, South Africa from papers (barbieri_2014b, 1kg, taylor_2020, just_2008, HARNEY_2023)
- LQ21 and LQ34 attached to the same node (large), sequences from the US, ESN, Botswana, ASW, lots of US in nodes attached to this node, LQ21 2 mutations away from a burkina faso sample from paper barbieri_2012a
```



metadata summaries:

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
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1)  +
  ggtitle("sample sizes of all used in trees")

# plot A2+(64) countries
a264_data <- meta[meta$hap_7char == "A2+(64)", ] # subset a2 data

a264_counts <- a264_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())

merged_a264 <- map_data %>%
  left_join(a264_counts, by = c("region" = "country_new"))

a264_plot <- ggplot(merged_a264, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("A2+(64) mjn sample sizes by country")

# plot B2b countries
b2b_data <- meta[meta$hap_3char == "B2b", ]
b2b_counts <- b2b_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_b2b <- map_data %>%
  left_join(b2b_counts, by = c("region" = "country_new"))
b2b_plot <- ggplot(merged_b2b, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("B2b mjn sample sizes by country")

# plot C1b countries
c1b_data <- meta[meta$hap_3char == "C1b", ]
c1b_counts <- c1b_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_c1b <- map_data %>%
  left_join(c1b_counts, by = c("region" = "country_new"))
c1b_plot <- ggplot(merged_c1b, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("C1b mjn sample sizes by country")

# plot D1 countries
d1_data <- meta[meta$hap_2char == "D1", ]
d1_counts <- d1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_d1 <- map_data %>%
  left_join(d1_counts, by = c("region" = "country_new"))
d1_plot <- ggplot(merged_d1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("D1 mjn sample sizes by country")

# plot H1b countries
h1b_data <- meta[meta$hap_3char == "H1b", ]
h1b_counts <- h1b_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_h1b <- map_data %>%
  left_join(h1b_counts, by = c("region" = "country_new"))
h1b_plot <- ggplot(merged_h1b, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("H1b iqtree sample sizes by country")


# plot L0 countries
#l0a2a
l0a2a_data <- meta[meta$hap_5char == "L0a2a", ]
l0a2a_counts <- l0a2a_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l0a2a <- map_data %>%
  left_join(l0a2a_counts, by = c("region" = "country_new"))
l0a2a_plot <- ggplot(merged_l0a2a, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L0a2a mjn sample sizes by country")

# L0d1a
l0d1a_data <- meta[meta$hap_5char == "L0d1a", ]
l0d1a_counts <- l0d1a_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l0d1a <- map_data %>%
  left_join(l0d1a_counts, by = c("region" = "country_new"))
l0d1a_plot <- ggplot(merged_l0d1a, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L0d1a mjn sample sizes per country")

# plot L1 countries
# L1c1
l1c1_data <- meta[meta$hap_4char == "L1c1", ]
l1c1_counts <- l1c1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l1c1 <- map_data %>%
  left_join(l1c1_counts, by = c("region" = "country_new"))
l1c1_plot <- ggplot(merged_l1c1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L1c1 mjn sample sizes by country")

# L1c2 or L1c4
l1c24_data <- meta[(meta$hap_4char == "L1c2" | meta$hap_4char == "L1c4" ), ]
l1c24_counts <- l1c24_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l1c24 <- map_data %>%
  left_join(l1c24_counts, by = c("region" = "country_new"))
l1c24_plot <- ggplot(merged_l1c24, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L1c'd mjn sample sizes by country")

# plot L2 countries
# L2a1
l2a1_data <- meta[meta$hap_4char == "L2a1", ]
l2a1f_data <- meta[meta$hap_5char == "L2a1f", ] # 167
l2a1q_data <- meta[meta$hap_5char == "L2a1q", ] # 4
l2a1a_data <- meta[meta$hap_5char == "L2a1a", ] # 266
l2a1d_data <- meta[meta$hap_5char == "L2a1d", ] # 34

l2a1a_counts <- l2a1a_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2a1a <- map_data %>%
  left_join(l2a1a_counts, by = c("region" = "country_new"))
l2a1a_plot <- ggplot(merged_l2a1a, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2a1a mjn sample sizes by country")

l2a1d_counts <- l2a1d_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2a1d <- map_data %>%
  left_join(l2a1d_counts, by = c("region" = "country_new"))
l2a1d_plot <- ggplot(merged_l2a1d, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2a1d mjn sample sizes by country")

l2a1f_counts <- l2a1f_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2a1f <- map_data %>%
  left_join(l2a1f_counts, by = c("region" = "country_new"))
l2a1f_plot <- ggplot(merged_l2a1f, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2a1f mjn sample sizes by country")

l2a1q_counts <- l2a1q_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2a1q <- map_data %>%
  left_join(l2a1q_counts, by = c("region" = "country_new"))
l2a1q_plot <- ggplot(merged_l2a1q, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2a1q mjn sample sizes by country")

# L2b1
l2b1_data <- meta[meta$hap_4char == "L2b1", ]
l2b1_counts <- l2b1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2b1 <- map_data %>%
  left_join(l2b1_counts, by = c("region" = "country_new"))
l2b1_plot <- ggplot(merged_l2b1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2b1 mjn sample sizes by country")

# L2c
l2c_data <- meta[meta$hap_3char == "L2c", ]
l2c_counts <- l2c_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2c <- map_data %>%
  left_join(l2b1_counts, by = c("region" = "country_new"))
l2c_plot <- ggplot(merged_l2c, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2c mjn sample sizes by country")

# plot L3 countries
# L3b1a
l3b1a_data <- meta[meta$hap_5char == "L3b1a", ]
l3b1a_counts <- l3b1a_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3b1a <- map_data %>%
  left_join(l3b1a_counts, by = c("region" = "country_new"))
l3b1a_plot <- ggplot(merged_l3b1a, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3b1a mjn sample sizes by country")

# L3d1
l3d1_data <- meta[meta$hap_4char == "L3d1", ]
l3d1_counts <- l3d1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3d1 <- map_data %>%
  left_join(l3d1_counts, by = c("region" = "country_new"))
l3d1_plot <- ggplot(merged_l3d1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3d1 mjn sample sizes by country")

# L3d3
l3d3_data <- meta[meta$hap_4char == "L3d3", ]
l3d3_counts <- l3d3_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3d3 <- map_data %>%
  left_join(l3d3_counts, by = c("region" = "country_new"))
l3d3_plot <- ggplot(merged_l3d3, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3d3 mjn sample sizes by country")

# L3d4
l3d4_data <- meta[meta$hap_4char == "L3d4", ]
l3d4_counts <- l3d4_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3d4 <- map_data %>%
  left_join(l3d4_counts, by = c("region" = "country_new"))
l3d4_plot <- ggplot(merged_l3d4, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3d4 mjn sample sizes by country")

# L3e1a1a
l3e1a1a_data <- meta[meta$hap_7char == "L3e1a1a", ]
l3e1a1a_counts <- l3e1a1a_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3e1a1a <- map_data %>%
  left_join(l3e1a1a_counts, by = c("region" = "country_new"))
l3e1a1a_plot <- ggplot(merged_l3e1a1a, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3e1a1a mjn sample sizes by country")

# L3e2
l3e2_data <- meta[meta$hap_4char == "L3e2", ]
l3e2_counts <- l3e2_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3e2 <- map_data %>%
  left_join(l3e2_counts, by = c("region" = "country_new"))
l3e2_plot <- ggplot(merged_l3e2, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3e2 mjn sample sizes by country")

# L3e3
l3e3_data <- meta[meta$hap_4char == "L3e3", ]
l3e3_counts <- l3e3_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3e3 <- map_data %>%
  left_join(l3e3_counts, by = c("region" = "country_new"))
l3e3_plot <- ggplot(merged_l3e3, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3e3 mjn sample sizes by country")

# L3f1b
l3f1b_data <- meta[meta$hap_5char == "L3f1b", ]
l3f1b_counts <- l3f1b_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3f1b <- map_data %>%
  left_join(l3f1b_counts, by = c("region" = "country_new"))
l3f1b_plot <- ggplot(merged_l3f1b, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3f1b mjn sample sizes by country")

pdf("sample_maps_mjn.pdf", width = 10, height = 8)
grid.arrange(full, ncol = 1)
grid.arrange(a264_plot, b2b_plot, ncol = 1)
grid.arrange(c1b_plot, d1_plot, ncol = 1)
grid.arrange(h1b_plot, l0a2a_plot, ncol = 1)
grid.arrange(l0d1a_plot, l1c1_plot, ncol = 1)
grid.arrange(l1c24_plot, l2a1a_plot, ncol = 1)
grid.arrange(l2a1d_plot, l2a1f_plot, ncol = 1)
grid.arrange(l2a1q_plot, l2b1_plot, ncol = 1)
grid.arrange(l2c_plot, l3b1a_plot, ncol = 1)
grid.arrange(l3d1_plot, l3d3_plot, ncol = 1)
grid.arrange(l3d4_plot, l3e1a1a_plot, ncol = 1)
grid.arrange(l3e2_plot, l3e3_plot, ncol = 1)
grid.arrange(l3f1b_plot, ncol = 1)
dev.off()

# also save metadata summary as text file
meta_list = list(a264_data, b2b_data, c1b_data, d1_data, h1b_data, l0a2a_data, l0d1a_data, l1c1_data, l1c24_data, l2a1a_data, l2a1d_data, l2a1f_data, l2a1q_data, l2b1_data, l2c_data, l3b1a_data, l3d1_data, l3d3_data, l3d4_data, l3e1a1a_data, l3e2_data, l3e3_data, l3f1b_data)
names(meta_list) <- c("a2_64", "b2b", "c1b", "d1", "h1b", "l0a2a", "l0d1a", "l1c1", "l1c24", "l2a1a", "l2a1d", "l2a1f", "l2a1q", "l2b1", "l2c", "l3b1a", "l3d1", "l3d3", "l3d4", "l3e1a1a", "l3e2", "l3e3", "l3f1b")

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



