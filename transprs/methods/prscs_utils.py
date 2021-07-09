import os
import tempfile
from pandas_plink import write_plink1_bin


def tmp_extract(processor):

    processor.sumstats[["CHR", "A1", "A2", "OR", "P"]].to_csv(
        "tmp_ss", sep="\t", index=False
    )

    write_plink1_bin(processor.population, "tmp.bim", verbose=False)
