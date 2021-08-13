import transprs as tprs


def Combine_inpop_methods(
    processor, methods, model, prs_col, trait_col, use_col, metric="r2_score"
):

    import itertools

    subsets = []
    for L in range(0, len(methods) + 1):
        for subset in itertools.combinations(methods, L):
            if len(subset) > 1:
                subsets.append(list(subset))

    for subset in subsets:
        tprs.combine.combine_methods(
            processor,
            methods=subset,
            trait_col=trait_col,
            key_ss="+".join(subset),
            use_col=use_col,
            prs_col=prs_col,
            model=model,
        )
        tprs.scoring.generate_prs(processor, method="+".join(subset))

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processor, method="+".join(subset), trait_col=trait_col, prs_col=prs_col
            )
