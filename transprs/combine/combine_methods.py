from transprs.combine import estimate_weighting
from transprs.combine import estimate_weighting_multipop, weighting_prs_multipop
import numpy as np
import pandas as pd
import time
import datetime
import subprocess



def combine_multipop(
    processors,
    methods,
    trait_col,
    key_ss,
    use_col,
    model="nnls",
    prs_col="SCORESUM",
    causal_snps=None,
):

    start_time = time.time()
    print("Estimating mixing weights...")

    # Calculate mixing weight

    mixing_weight = estimate_weighting_multipop(
        processors, methods, trait_col, model, prs_col
    )
    print("Mixing weights are: " + str(mixing_weight))
    print("Estimating mixing is done!")

    print("Adjusted " + use_col + "...")

    # Create snps list

    if causal_snps == None:
        causal_snps = []

    for processor, method in zip(processors, methods):

        adjusted_ss = pd.read_table(processor.adjusted_ss[method])

        causal_snps.append(set(adjusted_ss["SNP"].values))

    final_snps = set.intersection(*causal_snps)

    mixing_results = []

    betas_df = []

    for processor, method, weight in zip(processors, methods, mixing_weight):

        adjusted_ss = pd.read_table(processor.adjusted_ss[method])

        tmp = adjusted_ss[adjusted_ss.SNP.isin(final_snps)][use_col].values * weight

        mixing_results.append(tmp)

        betas_df.append(adjusted_ss[~adjusted_ss.SNP.isin(final_snps)])

    adjusted_beta = np.sum(mixing_results, axis=0)

    print("Adjusted BETA is done")

    if key_ss != None:
        name_combine = key_ss
    else:
        name_combine = "+".join(methods)

    sumstats = pd.read_table(processors[0].sumstats)

    adjusted_ss = {}

    adjusted_ss[name_combine] = sumstats.copy()

    adjusted_ss[name_combine] = adjusted_ss[name_combine][
        adjusted_ss[name_combine].SNP.isin(final_snps)
    ]

    adjusted_ss[name_combine][use_col] = adjusted_beta

    betas_df.append(adjusted_ss[name_combine])

    adjusted_ss[name_combine] = pd.concat(betas_df).sort_values(["CHR", "BP", "A1"])

    save_path = processor.workdir + "/adjusted_sumstats_" + name_combine
    adjusted_ss[name_combine].to_csv(save_path, sep="\t", index=False)

    print("The clumping result stores in .adjusted_ss['" + name_combine + "']!")

    processors[0].adjusted_ss[name_combine] = save_path

    processors[0].performance[name_combine] = {}

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )


def combine_prs_multipop(
    processors,
    methods,
    trait_col,
    key_ss,
    use_col,
    prs_col="SCORESUM",
):
    start_time = time.time()
    print("Estimating mixing weights...")
    
    # Calculate mixing weight

    adjusted_prs = weighting_prs_multipop(
        processors, methods, trait_col, prs_col
    )
    
    
    if key_ss != None:
        name_combine = key_ss
    else:
        name_combine = "+".join(methods)
        
    processors[0].prs_test[name_combine] = {}
    processors[0].prs_test[name_combine]["combine"] = processors[0].prs_test[methods[0]]["best_fit"].copy()
    processors[0].prs_test[name_combine]["combine"][prs_col] = adjusted_prs
    
    processors[0].performance[name_combine] = {}
    
    print("The " + name_combine + " result stores in .prs_test['" + name_combine + "']!")
    
    
    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
