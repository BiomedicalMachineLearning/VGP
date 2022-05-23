import transprs as tprs


def Inner_evaluate(
    processor, methods, metric, trait_col, prs_col, use_pca=True, validate=False
):
    for method in methods:
        if validate:
            assert method in processor.prs_test.keys(), (
                method + "is not in the .prs_test"
            )
        else:
            assert method in processor.prs_validation.keys(), (
                method + "is not in the .prs_validation"
            )

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processor,
                method=method,
                trait_col=trait_col,
                prs_col=prs_col,
                use_pca=use_pca,
                validate=validate,
            )
        elif metric == "coef_squared":
            tprs.metrics.coef_squared_evaluation(
                processor,
                method=method,
                trait_col=trait_col,
                prs_col=prs_col,
                use_pca=use_pca,
                validate=validate,
            )
