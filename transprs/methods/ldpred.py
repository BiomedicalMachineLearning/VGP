import subprocess
import time
import datetime
import os
import pandas as pd


def ldpred(
    bfile_ref,
    bfile_target,
    sumstat,
    N,
    ldf,
    ldr,
    h2,
    fraction_causal,
    A1="ALT",
    A2="REF",
    chr="#CHROM",
    pos="POS",
    effect="BETA",
    snpid="ID",
    gibbs_niter=1000,
    burnin=100,
    ldcoord_out="ldcoord",
    ldgibbs_out="ldgibbs",
    out="ldpred_score",
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
    print("LDpred is running...")

    if os.path.isfile(ldcoord_out):
        print("Output file (%s) already exists!  Deleting it" % ldcoord_out)
        os.remove(ldcoord_out)
    subprocess.call(
        """
            ldpred coord --gf %s --ssf %s --A1 %s --A2 %s --chr \"%s\" --pos %s --eff %s --rs %s --N %s --out %s
            """
        % (bfile_ref, sumstat, A1, A2, chr, pos, effect, snpid, N, ldcoord_out),
        shell=True,
    )

    if os.path.isfile(ldgibbs_out):
        print("Output file (%s) already exists!  Deleting it" % ldgibbs_out)
        os.remove(ldgibbs_out + ".pkl.gz")
    subprocess.call(
        """
            ldpred gibbs --cf %s --ldr %s --ldf %s --h2 %s --n-iter %s --n-burn-in %s --f %s --out %s
            """
        % (
            ldcoord_out,
            ldr,
            ldf,
            h2,
            gibbs_niter,
            burnin,
            fraction_causal,
            ldgibbs_out,
        ),
        shell=True,
    )

    print("Output to: %s_LDpred-inf.txt" % out)
    subprocess.call(
        """
            ldpred score --gf %s --rf %s --only-score --out %s
            """
        % (bfile_target, ldgibbs_out, out),
        shell=True,
    )

    res = pd.read_csv(out + "_LDpred-inf.txt")
    res.to_csv(out + "_LDpred.score", index=False, header=["IID", "PRS"])

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
