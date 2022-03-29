from subprocess import call, Popen, PIPE, STDOUT, CalledProcessError
import time
import datetime
import os
import pandas as pd
import time
import datetime
import pandas as pd
from rpy2 import robjects


def ldpred2(
    processor,
    ldref_path,
    map_file,
):

    """
    LDpred-2 method

    """

    start_time = time.time()

    print("LDpred-2 is running...")

    # Setup input
    robjects.globalenv["ldref_path"] = ldref_path
    robjects.globalenv["map_file"] = map_file
    robjects.globalenv["population_path"] = processor.validation + ".bed"
    robjects.globalenv["workdir"] = processor.workdir

    ss = pd.read_table(processor.sumstats)

    renamed_ss = ss.rename(columns={'CHR': 'chr', 'SNP': 'rsid',
                  'BP': 'pos', 'A1': 'a0',
                  'A2': 'a1', 'BETA': 'beta',
                  'P': 'p', 'N': 'n_eff',
                  'SE': 'beta_se', 'FRQ':'frq'})

    renamed_ss.to_csv("tmp_ss.ss",sep="\t",index=False)

    robjects.globalenv["CHRs"] = renamed_ss.chr.unique()

    # Main processing
    robjects.r(
    """
    options(stringsAsFactors = FALSE)
    NCORES=1
    library(bigsnpr)
    library(bigreadr)
    map_ldref <- readRDS(paste0(ldref_path,map_file))
    path2lnogwas<-'tmp_ss.ss'
    sumstats<-bigreadr::fread2(path2lnogwas)

    # Filter out hapmap SNPs

    info_snp <- snp_match(sumstats, map_ldref, match.min.prop = 0.01)
    sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))

    #is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0s.1 | sd_ldref < 0.05
    df_beta <- info_snp

    tmp <- tempfile(tmpdir = "tmp-data")
    for (chr in CHRs) {

        cat(chr, ".. ", sep = "")

        ## indices in 'df_beta'
        ind.chr <- which(df_beta$chr == chr)
        ## indices in 'map_ldref'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        ## indices in 'corr_chr'
        ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

        corr_chr <- readRDS(paste0(ldref_path,"/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]
 

        if (chr == CHRs[1]) {
          corr <- as_SFBM(corr_chr, tmp)
        } else {
          corr$add_columns(corr_chr, nrow(corr))
        }
      }
      
    (ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                  chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff,
                                  ncores = NCORES)))
    h2_est <- ldsc[["h2"]]
    
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
  
    """
    )

    df_beta = robjects.r["df_beta"]

    ldpred_ss = ss[ss.SNP.isin(df_beta["rsid.ss"].values)]

    ldpred_ss.BETA = robjects.r["beta_inf"]

    # Saving result

    save_path = processor.workdir + "/adjusted_sumstats_ldpred2"

    ldpred_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["ldpred2"] = save_path

    processor.tuning["ldpred2"] = {}
    
    processor.performance["ldpred2"] = {}

    print("The Ldpred-2 result stores in .adjusted_ss['ldpred2']!")

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

