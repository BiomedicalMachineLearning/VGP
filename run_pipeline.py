import os
import yaml
import transprs as tprs
import pandas as pd

def main():

    # Load settings
    assert os.path.exists('./settings.yml'), "Cannot find the settings file!"
    with open("settings.yml", "r") as stream:
        try:
            settings = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    # Check general format
    assert set(settings.keys()).issubset(['inputs', 
                                        'methods', 
                                        'metrics', 
                                        'visualization',
                                        'workdir',
                                        'preprocessing']), "Wrong format of .yml setting files!"
    
    # Check inputs
    list_headers = []
    for x in settings["inputs"]["list_sumstats"]:
        list_headers.append(pd.read_table(x, index_col=0, nrows=0).columns.tolist())
        
    for header in list_headers:
        assert set(header).issubset(["CHR",
                                    "BP",
                                    "SNP",
                                    "A1",
                                    "A2",
                                    "N",
                                    "SE",
                                    "P",
                                    "OR",
                                    "BETA",
                                    "INFO",
                                    "FRQ"]), "Wrong format of sumstats!"

    processors = []
    pre_params = settings["preprocessing"]

    for i, ss in enumerate(settings["inputs"]["list_sumstats"]):
        processor = tprs.read_input(prefix_bed=settings["inputs"]["genotype_prefix"], 
                            sumstats_path=ss,
                            workdir=settings["workdir"][i])
        phenotype = pd.read_table(settings["inputs"]["phenotype"],sep=" ")
        # phenotype.Height = np.random.randint(2, size=len(phenotype)).astype(str)
        processor.add_phenotype(phenotype)
        
        tprs.Preprocessing(processor, 
                        n_components=int(pre_params["n_components"]), 
                        k_folds=int(pre_params["k_folds"]), 
                        n_repeats=int(pre_params["n_repeats"]))
        
        processors.append(processor)

    if "clumping" in settings["methods"].keys():
        # Run the method PRS method: clumping
        for processor in processors:
            params = settings["methods"]["clumping"]
            tprs.methods.clumping(processor,
                                clump_p1=float(params["clump_p1"]),
                                clump_r2=float(params["clump_r2"]),
                                clump_kb=float(params["clump_kb"]))
            
            tprs.scoring.generate_prs(processor,method="clumping")
            
            if "coef_squared" in settings["metrics"]:
                tprs.metrics.coef_squared_evaluation(processor,
                                                    method="clumping",
                                                    trait_col=settings["inputs"]["trait_col"], 
                                                    prs_col=settings["inputs"]["prs_method"])
    
    if "double_weight" in settings["methods"].keys():
        # Run the method PRS method: clumping
        for processor in processors:
            params = settings["methods"]["double_weight"]
            tprs.methods.double_weight(processor,
                                top_choice=int(params["top_choice"]))
            
            tprs.scoring.generate_prs(processor,method="double_weight")
            
            if "coef_squared" in settings["metrics"]:
                tprs.metrics.coef_squared_evaluation(processor,
                                                    method="double_weight",
                                                    trait_col=settings["inputs"]["trait_col"], 
                                                    prs_col=settings["inputs"]["prs_method"])
    
    if len(processors) >= 2:
    
        if "coef_squared" in settings["metrics"]:
            tprs.Combine_multipop_methods(processors, methods=list(settings["methods"].keys()),
                                        trait_col=settings["inputs"]["trait_col"], 
                                        prs_col=settings["inputs"]["prs_method"],
                                        use_col="BETA",
                                        model="nnls",
                                        metric="coef_squared")
    
    # Transfer data to main processor

    for i, processor in enumerate(processors[1:]):
        for method in list(processor.adjusted_ss.keys()):
            new_name = method + "_" + str(i + 1)
            processor.adjusted_ss[new_name] = processor.adjusted_ss[method]
            processor.prs_results[new_name] = processor.prs_results[method]
            processor.performance[new_name] = processor.performance[method]
            del processor.adjusted_ss[method]
            del processor.prs_results[method]
            del processor.performance[method]

    for processor in processors[1:]:
        processors[0].adjusted_ss = {**processors[0].adjusted_ss, **processor.adjusted_ss}
        processors[0].prs_results = {**processors[0].prs_results, **processor.prs_results}
        processors[0].performance = {**processors[0].performance, **processor.performance}

    tprs.Generate_PRS(processors[0], methods=list(processors[0].adjusted_ss.keys()))
    
    if "coef_squared" in settings["metrics"]:
        tprs.Inner_evaluate(processors[0], methods=list(processors[0].adjusted_ss.keys()),
                trait_col=settings["inputs"]["trait_col"], 
                prs_col=settings["inputs"]["prs_method"],
                metric="coef_squared")

        if "box_plot" in settings["visualization"].keys():
            tprs.visualization.visualize_performance(processors[0], 
                                                    metric="coef_squared",
                                                    plot_type="box_plot",
                                                    figsize=eval(settings["visualization"]["box_plot"]["figsize"]),
                                                    fname=processors[0].workdir + "/coef_squared_evaluation_" + "box_plot" + ".png"
                                                    )

    print("======================================================")
    print("The pipeline is done!")


if __name__ == "__main__":
    main()
    
