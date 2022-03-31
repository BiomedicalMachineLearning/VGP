import os
import yaml
import transprs as tprs
import pandas as pd

def run_clumping(processor, population, settings):
    params = settings[population]["methods"]["clumping"]
    tprs.methods.clumping(processor,
                             clump_p1=float(params["clump_p1"]),
                             clump_r2=float(params["clump_r2"]),
                             clump_kb=float(params["clump_kb"]))
    
    tprs.scoring.generate_prs(processor,method="clumping")
    
    if settings["GENERAL_SETTING"]["metrics"] == "coef_squared":
        tprs.metrics.coef_squared_evaluation(processor,
                                                 method="clumping",
                                                 trait_col=settings["GENERAL_SETTING"]["inputs"]["trait_col"], 
                                                 prs_col=settings["GENERAL_SETTING"]["inputs"]["prs_method"])
                                   
    tprs.scoring.generate_prs(processor,method="clumping",validate=False)
                                   
    if settings["GENERAL_SETTING"]["metrics"] == "coef_squared":
        tprs.metrics.coef_squared_evaluation(processor,
                                                 method="clumping",
                                                 trait_col=settings["GENERAL_SETTING"]["inputs"]["trait_col"], 
                                                 prs_col=settings["GENERAL_SETTING"]["inputs"]["prs_method"],
                                                 validate=False)
        
def run_double_weight(processor, population, settings):
    params = settings[population]["methods"]["double_weight"]
    tprs.methods.double_weight(processor,
                             top_choice=int(params["top_choice"]))
    
    tprs.scoring.generate_prs(processor,method="double_weight")
    
    if settings["GENERAL_SETTING"]["metrics"] == "coef_squared":
        tprs.metrics.coef_squared_evaluation(processor,
                                                 method="double_weight",
                                                 trait_col=settings["GENERAL_SETTING"]["inputs"]["trait_col"], 
                                                 prs_col=settings["GENERAL_SETTING"]["inputs"]["prs_method"])
                                   
    tprs.scoring.generate_prs(processor,method="double_weight",validate=False)
                                   
    if settings["GENERAL_SETTING"]["metrics"] == "coef_squared":
        tprs.metrics.coef_squared_evaluation(processor,
                                                 method="double_weight",
                                                 trait_col=settings["GENERAL_SETTING"]["inputs"]["trait_col"], 
                                                 prs_col=settings["GENERAL_SETTING"]["inputs"]["prs_method"],
                                                 validate=False)
        
def run_pipeline(population,index, settings):
    processor = tprs.read_input(prefix_test=settings["GENERAL_SETTING"]["inputs"]["test_prefix"][index],
                    prefix_validation=settings["GENERAL_SETTING"]["inputs"]["validation_prefix"][index],
                    sumstats_path=settings[population]["inputs"]["sumstats"],
                    workdir=settings[population]["workdir"])
    phenotype = pd.read_table(settings["GENERAL_SETTING"]["inputs"]["phenotype_test"][index],sep="\s+")
    # phenotype.Height = np.random.randint(2, size=len(phenotype)).astype(str)
    processor.add_phenotype(processor.test, phenotype=phenotype)
    
    phenotype_val = pd.read_table(settings["GENERAL_SETTING"]["inputs"]["phenotype_validation"][index],sep="\s+")
    processor.add_phenotype(processor.validation, phenotype=phenotype_val) 
    
    tprs.Preprocessing(processor, 
                       n_components=int(settings[population]["preprocessing"]["n_components"])
                      )
    for method in settings["GENERAL_SETTING"]["inputs"]["methods"]:
        if method == "clumping":
            run_clumping(processor,population,settings)
        if method == "double_weight":
            run_double_weight(processor,population,settings)
    return processor
        
    
                    
        
def main():

    # Load settings
    assert os.path.exists('./settings.yml'), "Cannot find the settings file!"
    with open("settings.yml", "r") as stream:
        try:
            settings = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    populations = settings["GENERAL_SETTING"]["inputs"]["run_populations"]
    
    assert set(populations).issubset(settings.keys()), "Please check `run_populations` in the settings file again!"
    
    # Check inputs
    list_headers = []
    for x in populations:
        list_headers.append(pd.read_table(settings[x]["inputs"]["sumstats"], index_col=0, nrows=0,sep="\s+").columns.tolist())

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
                                    "FRQ","MAF"]), "Wrong format of sumstats!"

   
    all_processors = []
    for population in populations:
        processor = run_pipeline(population,index=0,settings=settings)
        all_processors.append(processor)
        
    import _pickle as cPickle
    with open(r"%s/all_processors.pickle" % settings["GENERAL_SETTING"]["inputs"]["output_folder"], "wb") as output_file:
        cPickle.dump(all_processors, output_file)


    print("======================================================")
    print("The pipeline is done!")


if __name__ == "__main__":
    main()
    
