from subprocess import call, Popen, PIPE, STDOUT, CalledProcessError
import time
import datetime
import os
import pandas as pd
import numpy as np
import transprs as tprs


def SBayesR(
    processor,
    ldm,
    pi=[0.95, 0.02, 0.02, 0.01],
    gamma=[0.0, 0.01, 0.1, 1],
    out_freq=100,
    chain_length=1000,
    burnin=100,
    random_state=1,
    threads=1,
):
    """
    SBayesR method

    Args:
        ldm (str): LD sparse matrix
        pi (str): pi, separated by ","
        gamma (str): pi, separated by ","
        sumstat (str): GWAS summary statistics
        h2 (float): heritability
        out (str): outputfile
        mcmc_niter (int, optional): number of MCMC iterations. Defaults to 1000.
        burnin (int, optional): number of burnin iterations. Defaults to 100.
        random_state (int, optional): random state. Defaults to 1.
        threads (int, optional): number of threads to run. Defaults to 1.
    """

    path = os.path.dirname(tprs.__file__)
    gctb_path = path + "/software/gctb"

    start_time = time.time()
    print("SBayesR is running...")

    ss = pd.read_table(processor.sumstats)
    tmp = ss[["SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]]
    tmp.to_csv("tmp.cojo", sep="\t", index=False)

    process_main = Popen(
        """
            %s --sbayes R --mldm %s \
                --pi %s \
                --gamma %s \
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
            ",".join(np.array(pi).astype(str)),
            ",".join(np.array(gamma).astype(str)),
            str(chain_length),
            str(burnin),
            str(out_freq),
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

    sbayesR_result = pd.read_table("tmp_result.snpRes", sep="\s+")
    sbayesR_result = sbayesR_result[["Name", "A1Effect"]]
    sbayesR_result.columns = ["SNP", "BETA"]
    adjusted_ss = pd.merge(ss, sbayesR_result, on="SNP")[
        ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA_y", "FRQ"]
    ]
    adjusted_ss.columns = [
        "CHR",
        "BP",
        "SNP",
        "A1",
        "A2",
        "N",
        "SE",
        "P",
        "BETA",
        "FRQ",
    ]
    save_path = processor.workdir + "/adjusted_sumstats_SBayesR"

    # Saving result
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["SBayesR"] = save_path

    processor.tuning["SBayesR"] = {}

    processor.performance["SBayesR"] = {}

    print("The SBayesR result stores in .adjusted_ss['SBayesR']!")

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
