# All conda environments used in paper created in accordance with UMN MSI conda best practices as of 2026
module load miniforge/24.3
module load conda/python3

# 1. environment for multi-sequence alignment (msa_env) + yaml creation
conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env python numpy pandas biopython seaborn matplotlib -c conda-forge -c bioconda
source activate /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/msa_env
conda env export > msa_env.yml

# 2. environment for iqtree
conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env python numpy pandas iqtree -c conda-forge -c bioconda
source activate $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env
conda env export > tree_env.yml

# 3. create environment for running interpopulation pop gen analyses
conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/12_diversity/st_env r-base=4.3 r-ape r-pegas r-hierfstat r-codetools -c conda-forge -c bioconda
source activate /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/st_env
R --vanilla
install.packages("haplotypes")
conda env export > st_env.yml

