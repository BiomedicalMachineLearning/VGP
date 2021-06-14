import pandas as pd
import os
import argparse

def qc(in_path):

    # Standard GWAS QC
    df = pd.read_csv(in_path, delim_whitespace=True)
    df = df[(df["INFO"] > 0.8) & (df["MAF"] > 0.01)]
    print("\nStandard GWAS QC-ed")

    # Remove duplicate SNPs
    df = df.drop_duplicates(subset='SNP',keep='first')
    print("\nDuplicate SNPs removed")
    
    # Remove ambigous SNPs
    df = df[
        ~((df['A1'] == 'A') & (df['A2'] == 'T') |
        (df['A1'] == 'T') & (df['A2'] == 'A') |
        (df['A1'] == 'G') & (df['A2'] == 'C') |
        (df['A1'] == 'C') & (df['A2'] == 'G'))
        ]
    print("\nAmbigous SNPs removed")
    
    # Write QC-ed file
    outfile = os.path.basename("../QC/orig/Height.gwas.txt.gz").split('.gwas.txt')[0] + '.QC.gz'
    df.to_csv(outfile, sep='\t', index=None)

    print(f"\nQC-ed base data writen to {outfile}.")

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input',
        dest='input',required=True,
        help='Input base data file or path')

    args = parser.parse_args()
    
    qc(args.input)
