import subprocess
import time
import datetime
import os
import pandas as pd
import transprs as tprs


def BayesR(
    bfile_ref, bfile_target, pheno, pi, h2, mcmc_niter=1000, burnin=50, out="bayesR_out"
):
    """
    BayesR

    Args:
        bfile_ref (str): genotype of reference sample
        pheno (str): phenotype of reference sample
        pi (float): pi
        h2 (float): heritability
        mcmc_niter (int, optional): number of MCMC iterations. Defaults to 1000.
        burnin (int, optional): number of MCMC burnins. Defaults to 50.
        out (str, optional): output. Defaults to "bayesR_out".

    Example:

        bfile_ref = 'popEUR_SAS_merged_train'
        bfile_target = 'popEUR_SAS_merged_test'
        pheno = 'pheno_train_scaled.txt'
        pi = 0.1
        h2 = 0.5
        mcmc_niter = 1000
        burnin = 50
        out = "bayesR_out"
    """
    #%%
    # bfile_ref = 'popEUR_SAS_merged_train'
    # bfile_target = 'popEUR_SAS_merged_test'
    # pheno = 'pheno_train_scaled.txt'
    # pi = 0.1
    # h2 = 0.5
    # mcmc_niter = 1000
    # burnin = 50
    # out = "bayesR_out"
    #%%
    start_time = time.time()
    print("BayesR is running...")

    path = os.path.dirname(tprs.__file__)
    gctb_path = path + "/softwares/gctb"

    if os.path.isfile(out):
        print("Output file (%s) already exists!  Deleting it" % out)
        os.remove(out)

    subprocess.call(
        """
            %s --bfile %s --pheno %s --bayes S --pi %s --hsq %s --chain-length %s --burn-in %s --out %s
            """
        % (gctb_path, bfile_ref, pheno, pi, h2, mcmc_niter, burnin, out),
        shell=True,
    )

    ss = pd.read_csv("%s.snpRes" % (out), sep="\s+")
    selcol = ["Chrom", "Name", "A1", "A2", "A1Effect", "SE", "PIP"]
    ss = ss[selcol]
    ss.columns = ["CHR", "ID", "A1", "A2", "BETA", "SE", "P"]
    ss.to_csv("%s.snpRes_sumstat" % (out), index=None, sep="\t")

    subprocess.call(
        """
            plink --bfile %s --score %s.snpRes_sumstat 2 3 5 --out %s
            """
        % (bfile_target, out, out),
        shell=True,
    )

    print("Output to: " + out + ".score_BayesR")
    res = pd.read_csv(out + ".profile", sep="\s+")
    res[["IID", "SCORE"]].to_csv(
        out + "_BayesR.score", index=False, header=["IID", "PRS"]
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
