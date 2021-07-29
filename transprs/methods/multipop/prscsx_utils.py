import os
import tempfile
from pandas_plink import write_plink1_bin


def tmp_extract(processors, use_col):

    for i, processor in enumerate(processors):
        processor.sumstats[["SNP", "A1", "A2", use_col, "P"]].to_csv(
            "tmp" + str(i) + "_ss", sep="\t", index=False
        )

    write_plink1_bin(processors[0].population, "tmp.bim", verbose=False)
