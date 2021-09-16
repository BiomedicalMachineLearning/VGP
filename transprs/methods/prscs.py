import subprocess
import pandas as pd
import time
import datetime
import os
import transprs


def prscs(
    processor,
    use_col,
    ldref_dir,
    N,
    phi=1e-2,
):
    start_time = time.time()
    # print("Extracting data...")
    # tmp_extract(processor, use_col=use_col)
    # print("Done extract data!")
    print("PRScs is running...")

    subprocess.call(
        """
        for chr in {21..22}; do \
            plink --bim %s.bim --chr $chr --make-just-bim --out tmp${chr}; \
        done
    """
        % (processor.population,),
        shell=True,
    )

    path = os.path.dirname(transprs.__file__)
    prscs_path = path + "/methods/prscs/PRScs.py"

    outdir = "tmp"

    try:
        os.mkdir(outdir)
    except:
        pass

    ss = pd.read_table(processor.sumstats)
    ss[["SNP", "A1", "A2", use_col, "P"]].to_csv("tmp_ss", sep="\t", index=False)

    for i in range(21, 23):
        print("Applying PRScs for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR = i
        subprocess.call(
            "python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s"
            % (
                prscs_path,
                ldref_dir,
                bim,
                "tmp_ss",
                str(N),
                str(CHR),
                str(phi),
                outdir,
            ),
            shell=True,
        )

        subprocess.call(
            "mv tmp_pst*chr%s.txt tmp_pst_chr%s.txt" % (str(i), str(i)), shell=True
        )
        print("PRScs for CHR " + str(i) + " is done!")

    print("Get adjusted_beta...")

    df_adj_ss = pd.read_table("./tmp_pst_chr21.txt", header=None)

    for i in range(22, 23):
        try:
            df_next = pd.read_table("./tmp_pst_chr" + str(i) + ".txt", header=None)
            df_adj_ss = pd.concat([df_adj_ss.reset_index(drop=True), df_next], axis=0)
        except:
            pass

    df_adj_ss = df_adj_ss.reset_index(drop=True)
    final_snps = list(set(df_adj_ss[1]) & set(ss["SNP"]))
    processor.adjusted_ss["PRScs"] = ss.copy()
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
