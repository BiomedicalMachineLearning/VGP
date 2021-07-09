from transprs.combine import estimate_weighting
import numpy as np
import time
import datetime
import subprocess


def combine_methods(processor, methods, trait_col, key_ss, use_col, prs_col="SCORE"):

    print("Estimating mixing weights...")
    mixing_weight = estimate_weighting(processor, methods, trait_col, prs_col)
    print("Mixing weights are: " + str(mixing_weight))
    print("Estimating mixing is done!")

    print("Adjusted BETA...")
    mixing_results = []
    for method in zip(methods, mixing_weight):
        mixing_results.append(
            processor.adjusted_ss[method[0]][use_col].values * method[1]
        )

    snps_list = []
    for method in methods:
        snps_list.append(set(processor.adjusted_ss[method]["SNP"].values))

    final_snps = set.intersection(*snps_list)

    mixing_results = []
    for method in zip(methods, mixing_weight):
        mixing_results.append(
            processor.adjusted_ss[method[0]][
                processor.adjusted_ss[method[0]].SNP.isin(final_snps)
            ][use_col].values
            * method[1]
        )

    adjusted_beta = np.sum(mixing_results, axis=0)

    print("Adjusted BETA is done")

    if key_ss != None:
        name_combine = key_ss
    else:
        name_combine = "+".join(methods)

    processor.adjusted_ss[name_combine] = processor.sumstats.copy()

    processor.adjusted_ss[name_combine] = processor.adjusted_ss[name_combine][
        processor.adjusted_ss[name_combine].SNP.isin(final_snps)
    ]

    processor.adjusted_ss[name_combine][use_col] = adjusted_beta

    print("The clumping result stores in .adjusted_ss['" + name_combine + "']!")

    processor.performance[name_combine] = {}

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
