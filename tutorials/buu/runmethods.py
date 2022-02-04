# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import os
import subprocess
import pandas as pd
import numpy as np
import transprs as tprs

os.chdir("/home/ubuntu/temp")


# %%
tprs.methods.ldpred(
    bfile_ref="popEUR_SAS_merged_train",
    bfile_target="popEUR_test",
    sumstat="popEUR_SAS_merged_train.PHENO1.glm.linear",
    N=810,
    ldr=500,
    h2=0.5,
    freq=0.02,
    out="ldpred_score",
)


# %%
tprs.methods.BayesR(
    bfile_ref="popEUR_SAS_merged_train",
    bfile_test="popEUR_SAS_merged_test",
    pheno="pheno_train_scaled.txt",
    pi=0.1,
    h2=0.5,
    mcmc_length=1000,
    burnin=50,
    out="bayesR_out",
)
# %%
tprs.methods.SBayesR(
    ldm="popSAS_sparseldm.ldm.sparse",
    pi="0.1,0.2,0.3,0.4",
    gamma="0.04,0.05,0.2,0.1",
    sumstat="popEUR_SAS_merged_train.PHENO1.glm.linear",
    bfile_target="popEUR_SAS_merged_test",
    h2=0.5,
    mcmc_length=1000,
    burnin=50,
    out="SbayesR_out",
)

# %%
ldm = "popSAS_sparseldm.ldm.sparse"
pi = "0.1,0.2,0.3,0.4"
gamma = "0.04,0.05,0.2,0.1"
sumstat = "popEUR_SAS_merged_train.PHENO1.glm.linear"
h2 = 0.5
bfile_target = "popEUR_SAS_merged_test"
out = "SBayesR"
A1_ss = "ALT"
A2_ss = "REF"
chr_ss = "#CHROM"
pos_ss = "POS"
effect_ss = "BETA"
SE_ss = "SE"
P_ss = "SE"
N_ss = "OBS_CT"
snpid_ss = "ID"
af_file = "popSAS.frq"
SNP_af = "SNP"
A1_af = "A1"
freq_af = "MAF"
mcmc_length = 1000
burnin = 100
random_state = 1
threads = 1

#%%
ss = pd.read_csv(sumstat, sep="\t")
af = pd.read_csv(af_file, sep="\s+")
af.rename({freq_af: "AF", A1_af: "A1", SNP_af: "SNP"}, axis=1, inplace=True)
snpdict_af = af.set_index(SNP_af).to_dict("index")
# %%
def get_af(snp_info, snpdict_af):
    if snp_info["SNP"] not in snpdict_af:
        return 0
    if snp_info["A1"] == snpdict_af[snp_info["SNP"]]["A1"]:
        return snpdict_af[snp_info["SNP"]]["AF"]
    else:
        return 1 - snpdict_af[snp_info["SNP"]]["AF"]


# %%
ss_c = ss.copy()
# %%
ss = ss_c.copy()
ss["freq"] = ss.apply(
    lambda row: get_af({"SNP": row[snpid_ss], "A1": row[A1_ss]}, snpdict_af), axis=1
)
ss = ss[[snpid_ss, A1_ss, A2_ss, "freq", effect_ss, SE_ss, P_ss, N_ss]]
ss.to_csv("%s_sbayesr" % (sumstat), index=None, sep=" ", header=None)
# %%
subprocess.call(
    """
        gctb --sbayes R --ldm %s --pi %s --gamma %s --gwas-summary %s --exclude-mhc --hsq %s --chain-length %s --burn-in %s --seed %s --thread %s --out %s
        """
    % (
        ldm,
        pi,
        gamma,
        sumstat + "_sbayesr",
        h2,
        mcmc_length,
        burnin,
        random_state,
        threads,
        out,
    ),
    shell=True,
)
# %%

ss = pd.read_csv("%s.snpRes" % (out), sep="\s+")
selcol = ["Chrom", "Name", "A1", "A2", "A1Effect", "SE", "PIP"]
ss = ss[selcol]
ss.columns = ["CHR", "ID", "A1", "A2", "BETA", "SE", "P"]
ss.to_csv("%s.snpRes_sumstat" % (out), index=None, sep="\t")

#%%
subprocess.call(
    """
        plink --bfile %s --score %s.snpRes_sumstat 2 3 5 --out %s
        """
    % (bfile_target, out, out),
    shell=True,
)

# %%
tprs.methods.BayesR(
    bfile_ref="popEUR_SAS_merged_train",
    bfile_target="popEUR_SAS_merged_test",
    pheno="pheno_train_scaled.txt",
    pi=0.1,
    h2=0.5,
    mcmc_length=1000,
    burnin=50,
    out="asdfadbayesR_out",
)


# %%
out = "ldpred_score"
res = pd.read_csv(out + "_LDpred-inf.txt")
res.to_csv(out + "_LDpred-inf.txt", index=False, header=["IID", "PRS"])
# %%
out = "SbayesR_out"
res = pd.read_csv(out + ".profile", sep="\s+")
res[["IID", "SCORE"]].to_csv(out + ".score_SBayesR", index=False, header=["IID", "PRS"])

# %%
tprs.methods.SBayesR(
    ldm="popSAS_sparseldm.ldm.sparse",
    pi="0.1,0.2,0.3,0.4",
    gamma="0.04,0.05,0.2,0.1",
    sumstat="popEUR_SAS_merged_train.PHENO1.glm.linear",
    bfile_target="popEUR_SAS_merged_test",
    h2=0.5,
    mcmc_length=1000,
    burnin=50,
    out="SbayesR_out",
)

#%%
path = os.path.dirname(tprs.__file__)
print(path)
prscs_path = path + "/methods/prscs/PRScs.py"
# %%
