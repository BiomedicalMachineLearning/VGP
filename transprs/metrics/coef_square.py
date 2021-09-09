from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler
from transprs.metrics.utils import correlation_evaluation
import numpy as np
import pandas as pd


def coef_squared_evaluation(
    processor,
    method,
    trait_col,
    prs_col,
    best_fit_key="best_fit",
    id_col="FID",
    scale=True,
    optimal_pval_key="optimal_pval",
    use_pca=True,
):

    try:
        del processor.prs_results[method][best_fit_key]
        del processor.prs_results[method][optimal_pval_key]
    except:
        pass

    phenotype = pd.read_table(processor.phenotype)

    # Do repeated k-fold
    results = {}
    for pval in processor.prs_results[method].keys():

        merged_df = pd.merge(processor.prs_results[method][pval], phenotype)
        if scale:
            scaler = MinMaxScaler()
            scaler.fit(merged_df[[prs_col, trait_col]])
            merged_df[[prs_col, trait_col]] = scaler.transform(
                merged_df[[prs_col, trait_col]]
            )
        score_list = correlation_evaluation(
            processor,
            merged_df,
            trait_col=trait_col,
            prs_col=prs_col,
            id_col=id_col,
        )
        results[pval] = score_list[0]

    # Get mean each fold
    mean_results = {}
    for key in results.keys():
        mean_results[key] = results[key]

    # Get best p-value
    best_key = max(mean_results, key=mean_results.get)

    print("The best fit p-value is " + best_key)

    processor.prs_results[method][best_fit_key] = processor.prs_results[method][
        best_key
    ]

    processor.prs_results[method][optimal_pval_key] = best_key

    print(
        "The best fit result is stored in processor.prs_results['"
        + method
        + "']['"
        + best_fit_key
        + "']"
    )

    processor.performance[method]["coef_square"] = results[best_key]

    print(
        "The best fit result is stored in processor.performance['"
        + method
        + "']['coef_square']"
    )
