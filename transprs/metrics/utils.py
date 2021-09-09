import numpy as np


def model_based_evaluation(
    processor, merged_df, trait_col, prs_col, model, metric, use_pca, id_col="FID"
):
    score_list = []
    for i in range(len(processor.dataset_repeated_split)):
        for j in range(len(processor.dataset_repeated_split[0])):
            train_inds = [
                x for k, x in enumerate(processor.dataset_repeated_split[i]) if k != j
            ]
            train_inds = [item for sublist in train_inds for item in sublist]
            test_inds = processor.dataset_repeated_split[i][j]

            if use_pca:
                X_train = np.array(
                    [
                        merged_df[merged_df[id_col].isin(train_inds)][
                            [prs_col] + [i for i in merged_df.columns if "PC" in i]
                        ]
                    ]
                )[0]
                X_test = np.array(
                    [
                        merged_df[merged_df[id_col].isin(test_inds)][
                            [prs_col] + [i for i in merged_df.columns if "PC" in i]
                        ]
                    ]
                )[0]

            else:
                X_train = np.array(
                    [merged_df[merged_df[id_col].isin(train_inds)][prs_col]]
                ).T

                X_test = np.array(
                    [merged_df[merged_df[id_col].isin(test_inds)][prs_col]]
                ).T

            Y_train = np.array(
                [merged_df[merged_df[id_col].isin(train_inds)][trait_col]]
            ).T

            Y_test = np.array(
                [merged_df[merged_df[id_col].isin(test_inds)][trait_col]]
            ).T

            model.fit(X_train, Y_train)
            Y_pred = model.predict(X_test)

            score_list.append(metric(Y_test, Y_pred))

    return score_list


def correlation_evaluation(processor, merged_df, trait_col, prs_col, id_col="FID"):
    score_list = []

    from scipy.stats import pearsonr

    R2 = pearsonr(list(merged_df[prs_col]), list(merged_df[trait_col]))[0] ** 2

    score_list.append(R2)
    score_list.append(R2)

    return score_list
