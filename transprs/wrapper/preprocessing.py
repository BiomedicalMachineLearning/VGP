def Preprocessing(processor, n_components, k_folds, n_repeats, info=0.9, id_col="FID"):
    # Run preprocessing
    processor.clean_snps()
    processor.extract_intersection()
    if "INFO" in processor.sumstats.columns:
        processor.filter_imputed(info=info)
    processor.check_beta_se()
    processor.flip_reverse()
    processor.compute_pca(n_components=n_components)
    processor.cross_validation_split(
        id_col=id_col, k_folds=k_folds, n_repeats=n_repeats
    )
