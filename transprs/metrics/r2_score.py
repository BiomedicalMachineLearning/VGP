from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from .utils import model_based_evaluation
import numpy as np


def r2_score_evaluation(
    processor, method, trait_col, prs_col, best_fit_key="best_fit", id_col="FID"
):

    # Do repeated k-fold
    results = {}
    for pval in processor.prs_results[method].keys():
        merged_df = pd.merge(processor.prs_results[method][pval], processor.phenotype)
        score_list = model_based_evaluation(
            processor,
            merged_df,
            trait_col=trait_col,
            model=LinearRegression(),
            metric=r2_score,
            prs_col=prs_col,
            id_col=id_col,
        )
        results[pval] = score_list

    # Get mean each fold
    mean_results = {}
    for key in results.keys():
        mean_results[key] = np.mean(results[key])

    # Get best p-value
    best_key = min(mean_results, key=mean_results.get)

    processor.prs_results[method][best_fit_key] = processor.prs_results[method][
        best_key
    ]
