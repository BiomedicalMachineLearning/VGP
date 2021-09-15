from sklearn.preprocessing import MinMaxScaler
from transprs.combine.utils import nonneg_lstsq, ols
import pandas as pd


def estimate_weighting(processor, methods, trait_col, model="ols", prs_col="SCORESUM"):

    prs_results = []
    for method in methods:
        prs_results.append(processor.prs_results[method]["best_fit"][prs_col])

    df_prs_all = pd.concat(prs_results, axis=1).values

    scaler = MinMaxScaler()
    scaler.fit(df_prs_all)
    df_prs_all = scaler.transform(df_prs_all)

    phenotype = pd.read_table(processor.phenotype)

    df_pheno = phenotype[trait_col].values.reshape(-1, 1)
    # scaler = MinMaxScaler()
    # scaler.fit(df_pheno)
    # df_pheno = scaler.transform(df_pheno)

    if model == "nnls":
        mixing_weight, intercepts = nonneg_lstsq(df_prs_all, df_pheno.reshape(1, -1)[0])
    else:
        mixing_weight, intercepts = ols(df_prs_all, df_pheno.reshape(1, -1)[0])

    return mixing_weight


def estimate_weighting_multipop(
    processors, methods, trait_col, model="ols", prs_col="SCORESUM"
):

    prs_results = []
    for processor, method in zip(processors, methods):
        prs_results.append(processor.prs_results[method]["best_fit"][prs_col])

    df_prs_all = pd.concat(prs_results, axis=1).values

    scaler = MinMaxScaler()
    scaler.fit(df_prs_all)
    df_prs_all = scaler.transform(df_prs_all)

    phenotype = pd.read_table(processor.phenotype)

    df_pheno = phenotype[trait_col].values.reshape(-1, 1)
    # scaler = MinMaxScaler()
    # scaler.fit(df_pheno)
    # df_pheno = scaler.transform(df_pheno)

    if model == "nnls":
        mixing_weight, intercepts = nonneg_lstsq(df_prs_all, df_pheno.reshape(1, -1)[0])
    else:
        mixing_weight, intercepts = ols(df_prs_all, df_pheno.reshape(1, -1)[0])

    return mixing_weight
