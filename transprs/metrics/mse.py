from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np


def mse(processor, method, use_phenotype, scale=False, verbose=True):
    merged_df = pd.merge(processor.prs_results[method], processor.phenotype)

    if scale:
        scaler = MinMaxScaler()
        scaler.fit(merged_df[["SCORE", use_phenotype]])
        merged_df[["SCORE", use_phenotype]] = scaler.transform(
            merged_df[["SCORE", use_phenotype]]
        )

    mse = mean_squared_error(merged_df["SCORE"], merged_df[use_phenotype])

    processor.performance[method]["mse"] = mse

    processor.performance[method]["rmse"] = np.sqrt(mse)

    if verbose:
        print(
            "The MSE performance stored in processor.performance['"
            + method
            + "']['mse'] and processor.performance['"
            + method
            + "']['rmse']"
        )
