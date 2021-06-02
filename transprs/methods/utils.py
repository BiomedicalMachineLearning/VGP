import os
import tempfile
from pandas_plink import write_plink1_bin


def tmp_extract(processor):
    # with tempfile.NamedTemporaryFile(dir="./", mode='r+',delete=False) as temp:
    processor.sumstats.to_csv("tmp_ss", sep="\t", index=False)
    write_plink1_bin(processor.population, "tmp.bed", verbose=False)
