from transprs.methods.prscs_utils import tmp_extract
import subprocess
import pandas as pd
import time
import datetime
import os
import transprs


def prscsx(
    processors,
    populations,
    use_col,
    ldref_dir,
    N,
    phi=1e-2,
):
    start_time = time.time()
    print("Extracting data...")
    tmp_extract(processors, use_col=use_col)
    print("Done extract data!")
    print("PRScsx is running...")

    n_processors = len(processors)

    subprocess.call(
        """
        for chr in {1..22}; do \
            plink --bim tmp.bim --chr $chr --make-just-bim --out tmp${chr}; \
        done
    """,
        shell=True,
    )

    path = os.path.dirname(transprs.__file__)
    prscs_path = path + "/methods/multipop/prscsx/PRScsx.py"

    ss = ""
    for i in range(0, len(processors)):
        ss = ss + ",./tmp" + str(i) + "_ss"
    ss = ss[1:]

    populations = ",".join(populations)

    N = ",".join(map(str, N))

    outdir = "./tmp"

    outname = "tmp"

    for i in range(1, 23):
        print("Applying PRScs for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR = i
        subprocess.call(
            """
            mkdir tmp
            python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --pop=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s --outname=%s"
            """
            % (
                prscs_path,
                ldref_dir,
                bim,
                ss,
                populations,
                N,
                str(CHR),
                str(phi),
                outdir,
                outname,
            ),
            shell=True,
        )
        subprocess.call(
            "mv tmp_pst*chr%s.txt tmp_pst_chr%s.txt" % (str(i), str(i)), shell=True
        )
        print("PRScs for CHR " + str(i) + " is done!")

    print("Get adjusted_beta...")

    df_adj_ss = pd.read_table("./tmp_pst_chr1.txt", header=None)

    for i in range(2, 23):
        try:
            df_next = pd.read_table("./tmp_pst_chr" + str(i) + ".txt", header=None)
            df_adj_ss = pd.concat([df_adj_ss.reset_index(drop=True), df_next], axis=0)
        except:
            pass

    df_adj_ss = df_adj_ss.reset_index(drop=True)
    final_snps = list(set(df_adj_ss[1]) & set(processor.sumstats["SNP"]))
    processor.adjusted_ss["PRScs"] = processor.sumstats.copy()
    processor.adjusted_ss["PRScs"] = processor.adjusted_ss["PRScs"][
        processor.adjusted_ss["PRScs"].SNP.isin(final_snps)
    ]

    processor.adjusted_ss["PRScs"][use_col] = df_adj_ss[5].values

    processor.performance["PRScs"] = {}

    print("The clumping result stores in .adjusted_ss['PRScs']!")

    subprocess.call(
        """
    rm ./tmp*
        """,
        shell=True,
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
