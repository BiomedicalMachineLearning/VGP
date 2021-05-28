import pandas as pd
import numpy as np


class DataProcessor(object):
    """
    Main processor class which contains all preprocessing method
    """

    def __init__(
        self,
        sumstats: pd.DataFrame,
        imputed_snps: np.ndarray,
    ):

        assert type(sumstats) == pd.DataFrame, "Sumstats is not DataFrame!"
        assert type(imputed_snps) == np.ndarray, "Imputed snps is not an numpy array!"

        self.sumstats = sumstats
        self.imputed_snps = imputed_snps

    def clean_snps(self):
        """
        Clean both sumstats and imputed data
        """
        # Remove duplicate SNPs
        self.sumstats = self.sumstats.drop_duplicates("SNP")

        self.imputed_snps = np.unique(self.imputed_snps)

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
            ~self.sumstats["CHR"].isin(np.array(range(1, 23)).astype(str))
        ] = "23"
        self.sumstats["CHR"].astype(int)

    def filter_imputed(self, info=0.9):
        self.sumstats = self.sumstats[self.sumstats["INFO"] >= info]

    def extract_intersection(self):
        """
        Extract intersection based on sumstats and imputed data
        """
        intersected_snps = np.intersect1d(self.sumstats["SNP"], self.imputed_snps)
        self.sumstats = self.sumstats[self.sumstats.SNP.isin(intersected_snps)]

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

    def check_beta_se(self):
        """
        Check the form of beta/se
        """
        return None

    def flip_reverse(self):
        """
        Flip/inverse snps
        """
        return None

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
