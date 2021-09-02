from transprs.combine import estimate_weighting
import numpy as np
import pandas as pd
import time
import datetime
import subprocess


def combine_methods(
    processor, methods, trait_col, key_ss, use_col, model="ols", prs_col="SCORE"
):

    start_time = time.time()
    print("Estimating mixing weights...")
    mixing_weight = estimate_weighting(processor, methods, trait_col, model, prs_col)
    print("Mixing weights are: " + str(mixing_weight))
    print("Estimating mixing is done!")

    print("Adjusted BETA...")
    mixing_results = []
    for method in zip(methods, mixing_weight):
        adjusted_ss = pd.read_table(processor.adjusted_ss[method[0]])
        mixing_results.append(adjusted_ss[use_col].values * method[1])

    snps_list = []
    for method in methods:
        adjusted_ss = pd.read_table(processor.adjusted_ss[method])
        snps_list.append(set(adjusted_ss["SNP"].values))

    final_snps = set.intersection(*snps_list)

    mixing_results = []
    for method in zip(methods, mixing_weight):

        adjusted_ss = pd.read_table(processor.adjusted_ss[method[0]])

        mixing_results.append(
            adjusted_ss[adjusted_ss.SNP.isin(final_snps)][use_col].values * method[1]
        )

    adjusted_beta = np.sum(mixing_results, axis=0)

    print("Adjusted BETA is done")

    if key_ss != None:
        name_combine = key_ss
    else:
        name_combine = "+".join(methods)

    sumstats = pd.read_table(processor.sumstats)

    adjusted_ss_combined = sumstats.copy()

    adjusted_ss_combined = adjusted_ss_combined[
        adjusted_ss_combined.SNP.isin(final_snps)
    ]

    adjusted_ss_combined[use_col] = adjusted_beta

    save_path = processor.workdir + "/adjusted_sumstats_" + name_combine

    adjusted_ss_combined.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss[name_combine] = save_path

    print("The clumping result stores in .adjusted_ss['" + name_combine + "']!")

    processor.performance[name_combine] = {}

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
