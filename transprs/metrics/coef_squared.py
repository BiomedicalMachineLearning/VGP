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
    validate=True,
):
    if validate:
        prs_results = processor.prs_validation
        phenotype = pd.read_table(processor.phenotype_val)
    else:
        prs_results = processor.prs_test
        phenotype = pd.read_table(processor.phenotype)

    try:
        del prs_results[method][best_fit_key]
        del prs_results[method][optimal_pval_key]
    except:
        pass

    # Do repeated k-fold
    results = {}
    for pval in prs_results[method].keys():

        merged_df = pd.merge(prs_results[method][pval], phenotype)
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
            use_pca=use_pca,
        )
        results[pval] = score_list[0]

    # Get mean each fold
    mean_results = {}
    for key in results.keys():
        mean_results[key] = results[key]

    # Get best p-value
    best_key = max(mean_results, key=mean_results.get)

    print("The best fit p-value is " + best_key)

    prs_results[method][best_fit_key] = prs_results[method][best_key]

    prs_results[method][optimal_pval_key] = best_key

    if validate:

        print(
            "The best fit result is stored in processor.prs_validation['"
            + method
            + "']['"
            + best_fit_key
            + "']"
        )
        processor.tuning[method]["coef_squared"] = results[best_key]

        print(
            "The best fit result is stored in processor.tuning['"
            + method
            + "']['coef_squared']"
        )
    else:
        print(
            "The best fit result is stored in processor.prs_test['"
            + method
            + "']['"
            + best_fit_key
            + "']"
        )
        processor.performance[method]["coef_squared"] = results[best_key]

        print(
            "The best fit result is stored in processor.performance['"
            + method
            + "']['coef_squared']"
        )
