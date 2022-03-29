from subprocess import call, Popen, PIPE, STDOUT, CalledProcessError
import time
import datetime
import os
import pandas as pd
import numpy as np
import transprs as tprs


def SBayesS(
    processor,
    ldm,
    pi = 0.05,
    out_freq = 100,
    chain_length=1000,
    burnin=100,
    random_state=1,
    threads=1,
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

    start_time = time.time()
    path = os.path.dirname(tprs.__file__)
    gctb_path = path + "/software/gctb"

    start_time = time.time()
    print("SBayesS is running...")

    ss = pd.read_table(processor.sumstats)
    tmp = ss[["SNP","A1","A2","FRQ","BETA","SE","P","N"]]
    tmp.to_csv("tmp.cojo", sep="\t",index=False)

    process_main = Popen(
        """
            %s --sbayes S --mldm %s \
                --pi %s \
                --gwas-summary tmp.cojo \
                --chain-length %s \
                --burn-in %s \
                --exclude-mhc \
                --out-freq %s \
                --out tmp_result
        """
        % (
            gctb_path,
            ldm,
            str(pi),
            str(chain_length),
            str(burnin),
            str(out_freq)
        ),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )

    with process_main.stdout:
        try:
            for line in iter(process_main.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")
    
    sbayesS_result = pd.read_table("tmp_result.snpRes",sep="\s+")
    sbayesS_result = sbayesS_result[["Name","A1Effect"]]
    sbayesS_result.columns = ["SNP","BETA"]
    adjusted_ss = pd.merge(ss,sbayesS_result,on="SNP")[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA_y', 'FRQ']]
    adjusted_ss.columns = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'FRQ']
    save_path = processor.workdir + "/adjusted_sumstats_SBayesS"

    # Saving result
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["SBayesS"] = save_path

    processor.tuning["SBayesS"] = {}

    processor.performance["SBayesS"] = {}

    print("The SBayesS result stores in .adjusted_ss['SBayesS']!")

    call(
        """
    rm ./tmp*
        """,
        shell=True,
    )


    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )

