import transprs as tprs


def Inner_evaluate(processor, methods, metric, trait_col, prs_col):
    for method in methods:
        assert method in processor.adjusted_ss.keys(), (
            method + "is not in the .adjusted_ss"
        )

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processor, method=method, trait_col=trait_col, prs_col=prs_col
            )
