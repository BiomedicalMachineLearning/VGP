import subprocess
import time
import datetime
import os
import pandas as pd
import transprs as tprs


def get_af(snp_info, snpdict_af):
    if snp_info["SNP"] not in snpdict_af:
        return 0
    if snp_info["A1"] == snpdict_af[snp_info["SNP"]]["A1"]:
        return snpdict_af[snp_info["SNP"]]["AF"]
    else:
        return 1 - snpdict_af[snp_info["SNP"]]["AF"]


def SBayesR(
    ldm,
    pi,
    gamma,
    sumstat,
    h2,
    bfile_target,
    out,
    A1_ss="ALT",
    A2_ss="REF",
    effect_ss="BETA",
    snpid_ss="ID",
    SE_ss="SE",
    P_ss="SE",
    N_ss="OBS_CT",
    af_file="",
    SNP_af="SNP",
    A1_af="A1",
    freq_af="MAF",
    mcmc_niter=1000,
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
    gctb_path = path + "/softwares/gctb"

    start_time = time.time()
    print("SBayesR is running...")

    # if os.path.isfile(out):
    #     print('Output file (%s) already exists!  Deleting it' % out)
    #     os.remove(out)
    if af_file != "":
        ss = pd.read_csv(sumstat, sep="\t")
        af = pd.read_csv(af_file, sep="\s+")
        af.rename({freq_af: "AF", A1_af: "A1", SNP_af: "SNP"}, axis=1, inplace=True)
        snpdict_af = af.set_index(SNP_af).to_dict("index")
        ss["freq"] = ss.apply(
            lambda row: get_af({"SNP": row[snpid_ss], "A1": row[A1_ss]}, snpdict_af),
            axis=1,
        )
        ss = ss[[snpid_ss, A1_ss, A2_ss, "freq", effect_ss, SE_ss, P_ss, N_ss]]
        ss.to_csv("%s_sbayesr" % (sumstat), index=None, sep=" ", header=None)

    subprocess.call(
        """
            %s --sbayes R --ldm %s --pi %s --gamma %s --gwas-summary %s --exclude-mhc --hsq %s --chain-length %s --burn-in %s --seed %s --thread %s --out %s
            """
        % (
            gctb_path,
            ldm,
            pi,
            gamma,
            sumstat + "_sbayesr",
            h2,
            mcmc_niter,
            burnin,
            random_state,
            threads,
            out,
        ),
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

    print("Output to: " + out + ".score_SBayesR")
    res = pd.read_csv(out + ".profile", sep="\s+")
    res[["IID", "SCORE"]].to_csv(
        out + "_SBayesR.score", index=False, header=["IID", "PRS"]
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
