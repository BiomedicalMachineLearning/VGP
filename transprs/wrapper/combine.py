import transprs as tprs
import itertools


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
        )
        tprs.scoring.generate_prs(processor, method="+".join(subset))

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processor, method="+".join(subset), trait_col=trait_col, prs_col=prs_col
            )
        elif metric == "coef_squared":
            tprs.metrics.coef_squared_evaluation(
                processor, method="+".join(subset), trait_col=trait_col, prs_col=prs_col
            )


def Combine_multipop_methods(
    processors, methods, prs_col, trait_col, use_col, metric="r2_score"
):

    subsets = list(itertools.product(*[methods] * len(processors)))

    new_subsets = []
    for subset in subsets:
        track = 0
        for i, method in enumerate(list(subset)):
            if method in processors[i].adjusted_ss.keys():
                track += 1
        if track == len(processors):
            new_subsets.append(subset)
    for subset in new_subsets:
        new_subset = []
        for i, method in enumerate(subset):
            new_subset.append(method + str(i))

        tprs.combine.combine_prs_multipop(
            processors,
            methods=subset,
            trait_col=trait_col,
            prs_col=prs_col,
            key_ss="+".join(new_subset),
            use_col=use_col,
        )
        # tprs.scoring.generate_prs(processors[0], method="+".join(subset))

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processors[0],
                method="+".join(new_subset),
                trait_col=trait_col,
                prs_col=prs_col,
                validate=False,
            )
        elif metric == "coef_squared":
            tprs.metrics.coef_squared_evaluation(
                processors[0],
                method="+".join(new_subset),
                trait_col=trait_col,
                prs_col=prs_col,
                validate=False,
            )
