import subprocess
import time
import datetime
import os
import pandas as pd
import transprs as tprs


def BayesR(processor, pi, h2, mcmc_niter=1000, burnin=50):
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
    print("Extracting data...")
    tmp_extract(processor)
    print("Done extract data!")
    print("BayesR is running...")

    subprocess.call(
        """
            gctb --bfile %s --pheno %s --bayes S --pi %s --hsq %s --chain-length %s --burn-in %s --out %s
            """
        % ("tmp", "tmp_phenotype", pi, h2, mcmc_niter, burnin, "tmp_bayesR_out"),
        shell=True,
    )

    res = pd.read_csv("%s.snpRes" % ("tmp_bayesR_out"), sep="\s+")
    res = res.reset_index(drop=True)
    final_snps = list(set(res["Name"]) & set(processor.sumstats["SNP"]))
    processor.adjusted_ss["bayesR"] = processor.sumstats.copy()
    processor.adjusted_ss["bayesR"] = processor.adjusted_ss["bayesR"][
        processor.adjusted_ss["bayesR"].SNP.isin(final_snps)
    ]

    res = res[res.Name.isin(final_snps)]

    processor.adjusted_ss["bayesR"][use_col] = res["Effect"].values

    processor.performance["bayesR"] = {}

    print("The bayesR result stores in .adjusted_ss['bayesR']!")

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
