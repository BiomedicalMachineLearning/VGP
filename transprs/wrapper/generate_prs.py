import transprs as tprs


def Generate_PRS(processor, methods):
    for method in methods:
        assert method in processor.adjusted_ss.keys(), (
            method + "is not in the .adjusted_ss"
        )
        tprs.scoring.generate_prs(processor, method=method)
