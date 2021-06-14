import argparse
import pandas as pd
import numpy as np
import subprocess
import os
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

def update_effect_size(base_path):
    base_df = pd.read_csv(base_path, delim_whitespace=True)
    base_df['BETA'] = np.log(base_df['OR'])
    out_file_name = "plink/" + os.path.basename(base_path) + ".transformed"
    base_df.to_csv(out_file_name, sep=' ', index=None)
    print('\nEffect size updated.')
    return out_file_name 

def clumping(target_path, updated_base, clump_p1, clump_r2, clump_kb):
    target = os.path.basename(target_path)
    subprocess.call(
        f"plink \
            --bfile {target_path} \
            --clump-p1 {clump_p1} \
            --clump-r2 {clump_r2} \
            --clump-kb {clump_kb} \
            --clump {updated_base} \
            --clump-snp-field SNP \
            --clump-field P \
            --out plink/{target}", shell=True
    )
    subprocess.call(
        "awk 'NR!=1{print $3}' " + 
        f"plink/{target}.clumped > plink/{target}.valid.snp",
        shell=True
    )
    print('\nClumping done.')

def generate_prs(target_path, updated_base):
    target = os.path.basename(target_path)
    subprocess.run(
        f"plink \
            --bfile {target_path} \
            --score {updated_base} 3 4 12 header \
            --extract plink/{target}.valid.snp \
            --out plink/{target}", shell=True
    )
    print("\nPRS generated.")

def evaluate(target_path, cov_file, pheno_file, pc_file):
    target = os.path.basename(target_path)
    phenotype = pd.read_csv(pheno_file, delim_whitespace=True)
    pcs = pd.read_csv(pc_file,header=None, delim_whitespace=True)
    pcs = pcs.rename(columns={**{0:'FID', 1:'IID'},**{i:'PC'+str(i-1) for i in range(2,8)}})
    covariate = pd.read_csv(cov_file, delim_whitespace=True)
    pheno = pd.merge(phenotype,covariate).merge(pcs)

    model0 = LinearRegression()
    model0.fit(pheno.drop(columns=['FID','IID','Height']),pheno['Height'])
    r2_0 = r2_score(pheno['Height'],model0.predict(pheno.drop(columns=['FID','IID','Height'])))

    prs = pd.read_csv(f"plink/{target}.profile",delim_whitespace=True)
    pheno_prs = prs[['FID','IID','SCORE']].merge(pheno)

    model_prs = LinearRegression()
    model_prs.fit(pheno_prs.drop(columns=['FID','IID','Height']),pheno_prs['Height'])
    r2_prs = r2_score(pheno_prs['Height'],model_prs.predict(pheno_prs.drop(columns=['FID','IID','Height'])))
    return r2_prs - r2_0

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-t","--target", 
        dest="target", type=str, required=True,
        help="Target input prefix (.bed, .bim, .fam files, should be QC-ed)")
    
    parser.add_argument("-b","--base", 
        dest="base", type=str, required=True,
        help="Base input (should be QC-ed)")
    
    parser.add_argument("-p","--pheno", 
        dest="pheno", type=str, required=True,
        help="Phenotype of samples")

    parser.add_argument("-c","--cov", 
        dest="cov", type=str, required=True,
        help="Covariates of samples")

    parser.add_argument("-pc","--pc", 
        dest="pc", type=str, required=True,
        help="Principal components")

    parser.add_argument("--clump-p1",
        dest="clump_p1", default=1, required=False,
        help="p-value threshold for a SNP to be included as an index SNP.\
            1 means all SNPs are included")
    
    parser.add_argument("--clump-r2",
        dest="clump_r2", default=0.1, required=False,
        help="SNPs having r2 higher than 0.1 with the index SNPs will be removed")
    
    parser.add_argument("--clump-kb",
        dest="clump_kb", default=250, required=False,
        help="SNPs within 250k of the index SNP are considered for clumping")

    args = parser.parse_args()

    os.makedirs('plink',exist_ok=True)

    updated_base_data = update_effect_size(args.base)

    clumping(target_path = args.target, updated_base = updated_base_data,
        clump_p1 = 1, clump_r2 = 0.1, clump_kb = 250)

    generate_prs(target_path = args.target, updated_base = updated_base_data)

    print("\nScore: " + str(evaluate(args.target, args.cov, args.pheno, args.pc)))