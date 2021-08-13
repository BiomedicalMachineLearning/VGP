import os
import tempfile
from pandas_plink import write_plink1_bin


def tmp_extract(processor, method=None):
    # with tempfile.NamedTemporaryFile(dir="./", mode='r+',delete=False) as temp:
    if method == None:
        processor.sumstats.to_csv("tmp_ss", sep="\t", index=False)
    else:
        processor.adjusted_ss[method].to_csv("tmp_ss", sep="\t", index=False)
    write_plink1_bin(processor.population, "tmp.bed", verbose=False)

    processor.phenotype[processor.phenotype.columns[:3]].to_csv(
        "tmp_phenotype", index=False, sep=" "
    )
