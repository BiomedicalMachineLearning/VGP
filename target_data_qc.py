import pandas as pd
import numpy as np
import os
import argparse
import subprocess

def QC(target_path, base_path):
    target = os.path.basename(target_path).split('.')[0]

    # Unzip the data
    subprocess.call(f'unzip {target_path} -d {target}.QC',shell=True)
    print("\nTarget data unzipped")

    # Standard GWAS QC
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target} \
            --maf 0.01 \
            --hwe 1e-6 \
            --geno 0.01 \
            --mind 0.01 \
            --write-snplist \
            --make-just-fam \
            --out {target}.QC/{target}.QC", shell=True
        )
    print("\nStandard GWAS QC-ed")

    # Remove highly correlated SNPs
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target} \
            --keep {target}.QC/{target}.QC.fam \
            --extract {target}.QC/{target}.QC.snplist \
            --indep-pairwise 200 50 0.25 \
            --out {target}.QC/{target}.QC", shell=True
    )
    print("\nHighly correlated SNPs removed")

    # Remove samples with high heterozygosity rates
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target} \
            --extract {target}.QC/{target}.QC.prune.in \
            --keep {target}.QC/{target}.QC.fam \
            --het \
            --out {target}.QC/{target}.QC", shell=True
    )
    het = pd.read_csv(f"{target}.QC/{target}.QC.het", delim_whitespace=True)
    valid = het[
        (het['F'] <= het['F'].mean() + 3*het['F'].std()) &
        (het['F'] >= het['F'].mean() - 3*het['F'].std())
        ]
    valid[['FID','IID']].to_csv(f'{target}.QC/{target}.valid.sample', sep='\t', index=False)
    print("\nSamples with high heterozygosity rates removed")

    # Flip strand
    ## Load sumstats & QC SNP list
    ### Read bim file
    bim = pd.read_csv(f"{target}.QC/{target}.bim",delim_whitespace=True, header=None)
    bim = bim.rename(columns={0:'CHR', 1:'SNP', 2:'CM', 3:'BP', 4:'B.A1', 5:'B.A2'})
    bim['B.A1'] = [nuc.upper() for nuc in bim['B.A1']]
    bim['B.A2'] = [nuc.upper() for nuc in bim['B.A2']]
    ### Read sumstats
    height = pd.read_csv(base_path, delim_whitespace=True)
    height['A1'] = [nuc.upper() for nuc in height['A1']]
    height['A2'] = [nuc.upper() for nuc in height['A2']]
    ### Read QC-ed list
    qc = pd.read_csv(f"{target}.QC/{target}.QC.snplist",delim_whitespace=True, header=None)
    qc = qc[0]

    ## Identify SNPs requiring strand flipping
    ### Merge summary statistic with target
    info = pd.merge(bim, height)
    qced_ind = np.intersect1d(info['SNP'],qc,return_indices=True)[1]
    qced_ind.sort()
    info = info.loc[qced_ind]
    ### Function for calculating the complementary allele
    def complement(x):
        switch = {"A":"T","T":"A","G":"C","C":"G"}
        if x in switch:
            return switch[x]
        else:
            return x
    ### Indentify match SNPs
    info_match = info[(
        info['A1'] == info['B.A1']) & 
        (info['A2'] == info['B.A2'])]['SNP']
    ### Identify complementary SNPs
    com_snps = info[
        (info['B.A1'].apply(complement) == info['A1']) &
        (info['B.A2'].apply(complement) == info['A2'])]['SNP']
    com_ind = np.intersect1d(bim['SNP'],com_snps,return_indices=True)[1]
    ### Update bim file
    bim.loc[com_ind,'B.A1'] = bim.loc[com_ind,'B.A1'].apply(complement)
    bim.loc[com_ind,'B.A2'] = bim.loc[com_ind,'B.A2'].apply(complement)

    ## Identify SNPs requiring recoding
    ### Identify SNPs requiring recoding
    recode_snps = info[
        (info['B.A1'] == info['A2']) & 
        (info['B.A2'] == info['A1'])]['SNP']
    recode_ind = np.intersect1d(bim['SNP'],recode_snps,return_indices=True)[1]
    ### Update bim file
    bim_copy = bim.copy()
    bim.loc[recode_ind,'B.A1'] = bim_copy.loc[recode_ind,'B.A2']
    bim.loc[recode_ind,'B.A2'] = bim_copy.loc[recode_ind,'B.A1']
    ### Identify SNPs requiring recoding & complement
    com_recode_snps = info[
        (info['B.A1'].apply(complement) == info['A2']) &
        (info['B.A2'].apply(complement) == info['A1'])]['SNP']
    com_recode_ind = np.intersect1d(bim['SNP'],com_recode_snps,return_indices=True)[1]
    ### Update bim file
    bim_copy = bim.copy()
    bim.loc[com_recode_ind,'B.A1'] = bim_copy.loc[com_recode_ind,'B.A2'].apply(complement)
    bim.loc[com_recode_ind,'B.A2'] = bim_copy.loc[com_recode_ind,'B.A1'].apply(complement)
    ### Write the updated bim file to .a1 file
    bim.to_csv(f'{target}.QC/{target}.a1', sep='\t', columns=['SNP','B.A1'], header=None, index=False)
    ## Indentify mismatch SNPs
    mismatch = np.setdiff1d(
        bim['SNP'].tolist(),
        info_match.tolist() + com_snps.tolist() + recode_snps.tolist() + com_recode_snps.tolist()
        )
    pd.Series(mismatch).to_csv(f'{target}.QC/{target}.mismatch',header=None,index=False)
    print(f'\n Strand flipped')

    # Remove samples with mismatched sex information
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target} \
            --extract {target}.QC/{target}.QC.prune.in \
            --keep {target}.QC/{target}.valid.sample \
            --check-sex \
            --out {target}.QC/{target}.QC", shell=True
    )
    dat = pd.read_csv(f'{target}.QC/{target}.QC.sexcheck', delim_whitespace=True)
    dat[dat['STATUS']=='OK'][['FID','IID']].to_csv(f'{target}.QC/{target}.QC.valid', sep='\t', header=None, index=False)
    print('\nSamples with mismatched sex information')

    # Remove related samples
    subprocess.call(
        f"plink \
        --bfile {target}.QC/{target} \
        --extract {target}.QC/{target}.QC.prune.in \
        --keep {target}.QC/{target}.QC.valid \
        --rel-cutoff 0.125 \
        --out {target}.QC/{target}.QC", shell=True
    )
    print('\nRelated samples removed')

    # Generate final QC-ed target data file
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target} \
            --make-bed \
            --keep {target}.QC/{target}.QC.rel.id \
            --out {target}.QC/{target}.QC \
            --extract {target}.QC/{target}.QC.snplist \
            --exclude {target}.QC/{target}.mismatch \
            --a1-allele {target}.QC/{target}.a1", shell=True
    )
    print(f'\nQC-ed target data written in {target}.QC folder.')

    # Compute principal components
    ## Prune
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target}.QC \
            --indep-pairwise 200 50 0.25 \
            --out {target}.QC/{target}", shell=True
    )
    ## Calculate 6 principal components
    subprocess.call(
        f"plink \
            --bfile {target}.QC/{target}.QC \
            --extract {target}.QC/{target}.prune.in \
            --pca 6 \
            --out {target}.QC/{target}", shell=True
    )
    print('\nPrincipal components computed.')

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input',
        dest='input',required=True,
        help='Input target data path')

    parser.add_argument("-b","--base", 
        dest="base", type=str, required=True,
        help="Base input (should be QC-ed)")
    
    args = parser.parse_args()

    QC(args.input, args.base)