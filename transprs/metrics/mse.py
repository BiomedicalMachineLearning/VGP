from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error
import pandas as pd


def mse(processor, method, use_phenotype, scale=False, verbose=True):
    merged_df = pd.merge(
        processor.prs_results[method], processor.phenotype[use_phenotype]
    )

    if scale:
        scaler = MinMaxScaler()
        scaler.fit(merged[["SCORE", use_phenotype]])

    processor.performance[method]["mse"] = mean_squared_error(
        merged["SCORE"], merged[use_phenotype]
    )

    if verbose:
        print(
            "The MSE performance stored in processor.performance['"
            + method
            + "']['mse']"
        )
