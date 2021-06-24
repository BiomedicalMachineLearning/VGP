# Basic tutorial to run the framework

### Prepare dataset

You need to download from these ggdrive links:
- Base data: [https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view)
- Target data: [https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view](https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view)

Only uncompress the Target data (EUR.zip)

### Install the requirements

I recommend to use conda to setup the environment

```
conda install -c bioconda plink
pip install pandas-plink
```

### Run the framework
##### 1. Read files

Currently, I only support this format for sumstats then you should make sure your input format is the same.

I will update it next few weeks to support all the variation of column names but the order of columns should be the same.

```python
# Load libraries
from pandas_plink import read_plink1_bin
import pandas as pd
import transprs as tprs

# Read the population file
population = read_plink1_bin("./data/EUR.bed","./data/EUR.bim","./data/EUR.fam")
# Read sumstats file
sumstats = pd.read_table("./data/Height.gwas.txt.gz")

# Create the DataProcessor object
processor = tprs.datasets.DataProcessor(sumstats=sumstats, population=population)

```

##### 2. Run preprocessing
```python
# Run preprocessing
processor.clean_snps()
processor.filter_imputed(info=0.9)
processor.extract_intersection()
processor.check_beta_se()
processor.flip_reverse()
processor.split_chromosomes()
processor.sort_snps_chr()
```
##### 3. Run the polygenic model to adjust BETA (OR)
```python
# Run the method PRS method: clumping
tprs.methods.clumping(processor)
```
##### 4. Generate Polygenic Risk Score
```python
# Generate PRS score
tprs.scoring.generate_prs(processor,use_col="OR",method="clumping")
```

The PRS score will store in `processor.prs_results['clumping']`
