def Preprocessing(processor, n_components, info=0.9, id_col="FID"):
    # Run preprocessing
    processor.clean_snps(processor.test)

    if hasattr(processor, "validation"):
        processor.clean_snps(processor.validation)

    processor.extract_intersection()
    if "INFO" in processor.sumstats.columns:
        processor.filter_imputed(info=info)
    processor.check_beta_se()

    processor.flip_reverse(processor.test)
    processor.compute_pca(processor.test, n_components=n_components)

    if hasattr(processor, "validation"):
        processor.flip_reverse(processor.validation)
        processor.compute_pca(processor.validation, n_components=n_components)

    processor.store_path()
