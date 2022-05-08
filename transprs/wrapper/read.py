from pandas_plink import read_plink1_bin
import pandas as pd
import transprs as tprs
def read_input(prefix_test, test_phenotype, sumstats_path, workdir,prefix_validation=None, validation_phenotype=None):
    # Read the population file

    print("Reading the Test genotype...")

    test = read_plink1_bin(prefix_test + ".bed",
                                 prefix_test + ".bim",
                                 prefix_test + ".fam")
    test.name = 'test'

    

    if prefix_validation is not None:
        print("Reading the Validation genotype...")
        validation = read_plink1_bin(prefix_validation + ".bed",
                                    prefix_validation + ".bim",
                                    prefix_validation + ".fam")
        validation.name = 'validation'

        
    
    sumstats = pd.read_table(sumstats_path,sep="\s+",index_col=0)

    if sumstats.index.name != None:
        sumstats = pd.read_table(sumstats_path,sep="\s+")

    processor = tprs.datasets.DataProcessor(sumstats=sumstats, test=test, validation=validation,workdir=workdir)

    phenotype = pd.read_table(test_phenotype,sep=" ")
    processor.add_phenotype(processor.test, phenotype=phenotype)
    
    if validation_phenotype is not None:
        phenotype_val = pd.read_table(validation_phenotype,sep=" ")
        processor.add_phenotype(processor.validation, phenotype=phenotype_val) 
        
    return processor