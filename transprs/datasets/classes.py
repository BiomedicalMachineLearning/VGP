import pandas as pd
import numpy as np
from xarray import DataArray
from collections import OrderedDict
from random import randrange
from transprs.datasets.utils import pre_reader, complement
import subprocess


class DataProcessor(object):
    """
    Main processor class which contains all preprocessing method
    """

    def __init__(
        self,
        sumstats: pd.DataFrame,
        population: DataArray,
    ):

        assert type(sumstats) == pd.DataFrame, "Sumstats is not DataFrame!"
        assert type(population) == DataArray, "Imputed snps is not an DataArray!"

        sumstats = pre_reader(sumstats)

        self.sumstats = sumstats
        self.population = population
        self.adjusted_ss = OrderedDict()
        self.prs_results = OrderedDict()
        self.performance = OrderedDict()

    def clean_snps(self):
        """
        Clean both sumstats and imputed data
        """
        # Remove duplicate SNPs
        self.sumstats = self.sumstats.drop_duplicates("SNP")

        self.imputed_snps = np.unique(self.population["variant"]["snp"])

        # Remove ambigous SNPs
        self.sumstats = self.sumstats[
            self.sumstats["A1"].apply(lambda x: len(x) == 1)
            & self.sumstats["A2"].apply(lambda x: len(x) == 1)
        ]

        self.sumstats = self.sumstats[
            self.sumstats["A1"].isin(["A", "C", "G", "T"])
            & self.sumstats["A2"].isin(["A", "C", "G", "T"])
        ]

        self.sumstats[
            ~(
                ((self.sumstats["A1"] == "A") & (self.sumstats["A2"] == "T"))
                | ((self.sumstats["A1"] == "T") & (self.sumstats["A2"] == "A"))
                | ((self.sumstats["A1"] == "G") & (self.sumstats["A2"] == "C"))
                | ((self.sumstats["A1"] == "C") & (self.sumstats["A2"] == "G"))
            )
        ]

        self.sumstats["CHR"][
            ~self.sumstats["CHR"].isin(np.array(range(1, 23)).astype(int))
        ] = "23"
        self.sumstats["CHR"].astype(int)

    def filter_imputed(self, info=0.9):
        if "INFO" in self.sumstats:
            self.sumstats = self.sumstats[self.sumstats["INFO"] >= info]
        else:
            print("There is no 'INFO' column in sumstats!")

    def extract_intersection(self):
        """
        Extract intersection based on sumstats and imputed data
        """
        intersected_snps = np.intersect1d(self.sumstats["SNP"], self.imputed_snps)
        self.sumstats = self.sumstats[self.sumstats.SNP.isin(intersected_snps)]

    def check_beta_se(self):
        """
        Check the form of beta/se
        """
        if "OR" in self.sumstats.columns:
            self.sumstats["OR"] = np.log(self.sumstats["OR"])

    def split_chromosomes(self):
        """
        Split chromosomes
        """
        self.splited_sumstats = []
        for chromosome in self.sumstats["CHR"].unique():
            self.splited_sumstats.append(
                self.sumstats[self.sumstats["CHR"].isin([chromosome])]
            )

    def sort_snps_chr(self):
        """
        Sort SNPs by chr
        """
        for i, chr_df in enumerate(self.splited_sumstats):
            self.splited_sumstats[i] = chr_df.sort_values("SNP")

    def flip_reverse(self):
        """
        Flip/inverse snps
        """

        # Extract bim file
        bim = self.population.variant.to_dataframe()
        bim.columns = bim.columns.str.upper()
        bim = bim.rename(columns={"A0": "B.A1", "A1": "B.A2"})

        # Merge with sumstats
        info = pd.merge(bim, self.sumstats, on="SNP", how="inner")

        # Indentify match SNPs
        info_match = info[(info["A1"] == info["B.A1"]) & (info["A2"] == info["B.A2"])][
            "SNP"
        ]

        # Identify complementary SNPs
        ## Identify complementary SNPs
        com_snps = info[
            (info["B.A1"].apply(complement) == info["A1"])
            & (info["B.A2"].apply(complement) == info["A2"])
        ]["SNP"]
        com_ind = np.intersect1d(bim["SNP"], com_snps, return_indices=True)[1]
        com_ind = np.core.defchararray.add("variant", com_ind.astype(str))

        ## Update bim file
        bim.loc[com_ind, "B.A1"] = bim.loc[com_ind, "B.A1"].apply(complement)
        bim.loc[com_ind, "B.A2"] = bim.loc[com_ind, "B.A2"].apply(complement)

        # Identify SNPs requiring recoding
        ## Identify SNPs requiring recoding
        recode_snps = info[(info["B.A1"] == info["A2"]) & (info["B.A2"] == info["A1"])][
            "SNP"
        ]
        recode_ind = np.intersect1d(bim["SNP"], recode_snps, return_indices=True)[1]
        recode_ind = np.core.defchararray.add("variant", recode_ind.astype(str))

        ## Update bim file
        bim_copy = bim.copy()
        bim.loc[recode_ind, "B.A1"] = bim_copy.loc[recode_ind, "B.A2"]
        bim.loc[recode_ind, "B.A2"] = bim_copy.loc[recode_ind, "B.A1"]

        # Identify SNPs requiring recoding & complement
        ## Identify SNPs requiring recoding & complement
        com_recode_snps = info[
            (info["B.A1"].apply(complement) == info["A2"])
            & (info["B.A2"].apply(complement) == info["A1"])
        ]["SNP"]
        com_recode_ind = np.intersect1d(
            bim["SNP"], com_recode_snps, return_indices=True
        )[1]
        com_recode_ind = np.core.defchararray.add("variant", com_recode_ind.astype(str))

        ## Update bim file
        bim_copy = bim.copy()
        bim.loc[com_recode_ind, "B.A1"] = bim_copy.loc[com_recode_ind, "B.A2"].apply(
            complement
        )
        bim.loc[com_recode_ind, "B.A2"] = bim_copy.loc[com_recode_ind, "B.A1"].apply(
            complement
        )

        # Update result to population
        self.population.variant["a0"] = bim["B.A1"]
        self.population.variant["a1"] = bim["B.A2"]

        ## Indentify mismatch SNPs
        mismatch = np.setdiff1d(
            bim["SNP"].tolist(),
            info_match.tolist()
            + com_snps.tolist()
            + recode_snps.tolist()
            + com_recode_snps.tolist(),
        )

        self.sumstats = self.sumstats[~self.sumstats["SNP"].isin(mismatch)].reset_index(
            drop=True
        )

    def add_phenotype(self, phenotype: pd.DataFrame, id_col="FID"):

        self.phenotype = phenotype
        print("Phenotype stored in .phenotype")
        self.population = self.population.where(
            self.population.fid.isin(phenotype[id_col]), drop=True
        )

    def cross_validation_split(self, k_folds=5, n_repeats=10, id_col="FID"):

        self.dataset_repeated_split = []

        for i in range(n_repeats):
            dataset_split = []
            dataset_copy = list(self.phenotype[id_col])
            fold_size = int(len(self.phenotype[id_col]) / k_folds)
            for j in range(k_folds):
                fold = []
                while len(fold) < fold_size:
                    index = randrange(len(dataset_copy))
                    fold.append(dataset_copy.pop(index))
                dataset_split.append(fold)

            self.dataset_repeated_split.append(dataset_split)

        print("The splitted indexes are stored in .dataset_repeated_split")

    def compute_pca(self, n_components, id_col="FID"):
        from transprs.utils import tmp_extract

        tmp_extract(self)

        subprocess.call(
            """
        plink \
            --bfile tmp \
            --indep-pairwise 200 50 0.25 \
            --out tmp

        plink \
            --bfile tmp \
            --extract tmp.prune.in \
            --pca %s \
            --out tmp
        """
            % (str(n_components)),
            shell=True,
        )

        pca = pd.read_table("tmp.eigenvec", sep=" ", header=None)

        pca.columns = [id_col, "IID"] + ["PC" + str(x) for x in range(0, n_components)]

        self.phenotype = pd.merge(self.phenotype, pca)

        print("PCA result is stored in .phenotype")

        subprocess.call(
            """
            rm ./tmp*
                """,
            shell=True,
        )

    # def estimate_heritability(self):
    #     """
    #     Estimate heritability
    #     """
    #     return result

    # def estimate_genetic_corr(self):
    #     """
    #     Estimate genetic correlation
    #     """
    #     return result
