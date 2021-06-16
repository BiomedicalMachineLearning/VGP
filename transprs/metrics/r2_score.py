from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np


def r2_score(processor):

    score_list = []
    for i in range(len(processor.dataset_repeated_split)):
        for j in range(len(processor.dataset_repeated_split[0])):
            train_inds = [
                x for k, x in enumerate(processor.dataset_repeated_split[i]) if k != j
            ]
            train_inds = [item for sublist in train_inds for item in sublist]
            test_inds = processor.dataset_repeated_split[i][j]

            X_train = np.array(
                [
                    processor.phenotype[processor.phenotype["FID"].isin(train_inds)][
                        "Height"
                    ]
                ]
            ).T
            Y_train = np.array(
                [
                    processor.phenotype[processor.phenotype["FID"].isin(train_inds)][
                        "Height"
                    ]
                ]
            ).T

            X_test = np.array(
                [
                    processor.phenotype[processor.phenotype["FID"].isin(test_inds)][
                        "Height"
                    ]
                ]
            ).T
            Y_test = np.array(
                [
                    processor.phenotype[processor.phenotype["FID"].isin(test_inds)][
                        "Height"
                    ]
                ]
            ).T

            model = LinearRegression()
            model.fit(X_train, Y_train)
            Y_pred = model.predict(X_test)
            score_list.append(r2_score(Y_test, Y_pred))
