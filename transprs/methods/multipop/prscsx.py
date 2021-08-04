from transprs.methods.multipop.prscsx_utils import tmp_extract
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
    prscsx_path = path + "/methods/multipop/prscsx/PRScsx.py"

    ss = ""
    for i in range(0, len(processors)):
        ss = ss + ",./tmp" + str(i) + "_ss"
    ss = ss[1:]

    populations = ",".join(populations)

    N = ",".join(map(str, N))

    outdir = "./tmp"

    outname = "tmp"

    subprocess.call("mkdir tmp", shell=True)

    for i in range(1, 23):
        print("Applying PRScsx for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR = i
        subprocess.call(
            """
            python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --pop=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s --out_name=%s
            """
            % (
                prscsx_path,
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

        for pop in populations.split(","):
            subprocess.call(
                "mv tmp/tmp_%s_pst*chr%s.txt tmp/tmp_%s_pst_chr%s.txt"
                % (pop, str(i), pop, str(i)),
                shell=True,
            )
        print("PRScsx for CHR " + str(i) + " is done!")

    print("Get adjusted_beta...")

    for processor, pop in zip(processors, populations.split(",")):

        df_adj_ss = pd.read_table("tmp/tmp_" + pop + "_pst_chr1.txt", header=None)

        for i in range(2, 23):
            try:
                df_next = pd.read_table(
                    "tmp/tmp_" + pop + "_pst_chr" + str(i) + ".txt", header=None
                )
                df_adj_ss = pd.concat(
                    [df_adj_ss.reset_index(drop=True), df_next], axis=0
                )
            except:
                pass

        df_adj_ss = df_adj_ss.reset_index(drop=True)
        final_snps = list(set(df_adj_ss[1]) & set(processor.sumstats["SNP"]))
        processor.adjusted_ss["PRScsx"] = processor.sumstats.copy()
        processor.adjusted_ss["PRScsx"] = processor.adjusted_ss["PRScsx"][
            processor.adjusted_ss["PRScsx"].SNP.isin(final_snps)
        ]

        processor.adjusted_ss["PRScsx"][use_col] = df_adj_ss[5].values

        processor.performance["PRScsx"] = {}

    print("The clumping result stores in .adjusted_ss['PRScsx']!")

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
