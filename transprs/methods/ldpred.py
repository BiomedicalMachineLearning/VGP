import subprocess
import time
import datetime
import os
import pandas as pd
from transprs.utils import tmp_extract


def ldpred(
    processor,
    use_col,
    N,
    ldf,
    ldr,
    h2,
    fraction_causal,
    effect="BETA",
    gibbs_niter=1000,
    burnin=100,
):

    """
    LDpred method

    Args:
        bfile_ref: plink-formatted bfile for input reference genotype
        bfile_target: plink-formatted bfile for input target genotype
        sumstat: summary statistics
        N: sample size
        ldf: ld files
        h2: heritability
        fraction_causal: proportion of causal snps
        A1 (str, optional): A1 allele. Defaults to "ALT".
        A2 (str, optional): A2 allele. Defaults to "REF".
        chr (str, optional): chromosome column name. Defaults to "#CHROM".
        pos (str, optional): position column name. Defaults to "POS".
        effect (str, optional): effect size column name. Defaults to "BETA".
        snpid (str, optional): SNP ID column name. Defaults to "ID".
        ldr (int, optional): LD radius. Defaults to 500.
        ldcoord_out (str, optional): ldcoord filename out. Defaults to "ldcoord".
        ldgibbs_out (str, optional): ld gibbs out. Defaults to "ldgibbs".
        out (str, optional): output file for PRS from LDpred. Defaults to "ldpred_score".

    Example:
        bfile_ref = 'popEUR_SAS_merged_train'
        bfile_target = 'popEUR_test'
        sumstat = 'popEUR_SAS_merged_train.PHENO1.glm.linear'
        A1="ALT"
        A2="REF"
        chr="#CHROM"
        pos="POS"
        effect="BETA"
        snpid="ID"
        N=810
        ldr=500
        ldf='LDF_EUR_SAS'
        h2=0.5
        fraction_causal=0.02
        ldcoord_out="ldcoord"
        ldgibbs_out="ldgibbs"
        out="ldpred_score"

    """

    start_time = time.time()
    print("Extracting data...")
    tmp_extract(processor)
    print("Done extract data!")
    print("LDpred is running...")

    subprocess.call(
        """
            ldpred coord --gf %s --ssf %s --A1 %s --A2 %s --chr %s --pos %s --eff %s --rs %s --N %s --out %s
            """
        % (
            "tmp",
            "tmp_ss",
            "A1",
            "A2",
            "CHR",
            "BP",
            use_col,
            "SNP",
            str(N),
            "tmp_ldcoord",
        ),
        shell=True,
    )

    subprocess.call(
        """
            ldpred gibbs --cf %s --ldr %s --ldf %s --h2 %s --n-iter %s --n-burn-in %s --f %s --out %s
            """
        % (
            "tmp_ldcoord",
            ldr,
            ldf,
            h2,
            gibbs_niter,
            burnin,
            fraction_causal,
            "tmp_ldgibbs_out",
        ),
        shell=True,
    )

    res = pd.read_table("tmp_ldgibbs_out_LDpred-inf.txt", sep="\s+")
    res = res.reset_index(drop=True)
    final_snps = list(set(res["sid"]) & set(processor.sumstats["SNP"]))
    processor.adjusted_ss["ldpred"] = processor.sumstats.copy()
    processor.adjusted_ss["ldpred"] = processor.adjusted_ss["ldpred"][
        processor.adjusted_ss["ldpred"].SNP.isin(final_snps)
    ]

    processor.adjusted_ss["ldpred"][use_col] = res[res.columns[-1]].values

    processor.performance["ldpred"] = {}

    print("The ldpred result stores in .adjusted_ss['ldpred']!")

    subprocess.call(
        """
    rm ./tmp*
        """,
        shell=True,
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
