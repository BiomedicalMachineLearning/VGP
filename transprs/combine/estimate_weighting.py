from sklearn.preprocessing import MinMaxScaler
from .utils import nonneg_lstsq
import pandas as pd


def estimate_weighting(processor, methods, trait_col, prs_col="SCORE"):

    prs_results = []
    for method in methods:
        prs_results.append(processor.prs_results[method]["best_fit"][prs_col])

    df_prs_all = pd.concat(prs_results, axis=1).values

    scaler = MinMaxScaler()
    scaler.fit(df_prs_all)
    df_prs_all = scaler.transform(df_prs_all)

    df_pheno = processor.phenotype[trait_col].values.reshape(-1, 1)
    scaler = MinMaxScaler()
    scaler.fit(df_pheno)
    df_pheno = scaler.transform(df_pheno)

    mixing_weight, intercepts = nonneg_lstsq(df_prs_all, df_pheno.reshape(1, -1)[0])

    return mixing_weight
