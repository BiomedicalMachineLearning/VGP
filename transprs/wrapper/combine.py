import transprs as tprs
from transprs.combine import multipop


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


def Combine_multipop_methods(
    processors, methods, model, prs_col, trait_col, use_col, metric="r2_score"
):

    subsets = []
    for L in range(0, len(methods)):
        for H in range(0, len(methods)):
            subsets.append([methods[L], methods[H]])

    for subset in subsets:
        multipop.combine_multipop(
            [processor, processor2],
            methods=subset,
            trait_col=trait_col,
            prs_col=prs_col,
            key_ss="+".join(subset),
            use_col=use_col,
        )
        tprs.scoring.generate_prs(processor, method="+".join(subset))

        if metric == "r2_score":
            tprs.metrics.r2_score_evaluation(
                processor, method="+".join(subset), trait_col=trait_col, prs_col=prs_col
            )
