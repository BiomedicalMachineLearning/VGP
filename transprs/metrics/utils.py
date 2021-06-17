import numpy as np


def model_based_evaluation(
    processor, merged_df, trait_col, prs_col, model, metric, id_col="FID"
):
    score_list = []
    for i in range(len(processor.dataset_repeated_split)):
        for j in range(len(processor.dataset_repeated_split[0])):
            train_inds = [
                x for k, x in enumerate(processor.dataset_repeated_split[i]) if k != j
            ]
            train_inds = [item for sublist in train_inds for item in sublist]
            test_inds = processor.dataset_repeated_split[i][j]

            X_train = np.array(
                [merged_df[merged_df[id_col].isin(train_inds)][prs_col]]
            ).T

            Y_train = np.array(
                [merged_df[merged_df[id_col].isin(train_inds)][trait_col]]
            ).T

            X_test = np.array([merged_df[merged_df[id_col].isin(test_inds)][prs_col]]).T
            Y_test = np.array(
                [merged_df[merged_df[id_col].isin(test_inds)][trait_col]]
            ).T

            model.fit(X_train, Y_train)
            Y_pred = model.predict(X_test)

            score_list.append(metric(Y_test, Y_pred))

    return score_list
