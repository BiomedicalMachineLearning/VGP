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
    print("PRScsx is running...")

    n_processors = len(processors)

    subprocess.call(
        """
        for chr in {1..22}; do \
            plink --bim %s.bim --chr $chr --make-just-bim --out tmp${chr}; \
        done
        """
        % processors[0].population,
        shell=True,
    )

    path = os.path.dirname(transprs.__file__)
    prscsx_path = path + "/methods/multipop/prscsx/PRScsx.py"

    ss = ""
    for i in range(0, len(processors)):
        tmp = pd.read_table(processors[i].sumstats)
        tmp[["SNP", "A1", "A2", "BETA", "P"]].to_csv(
            "tmp_sumstat" + str(i), sep=" ", index=False
        )
        ss = ss + ",tmp_sumstat" + str(i)
    ss = ss[1:]

    tmp_populations = ",".join(populations)

    tmp_N = ",".join(map(str, N))

    outdir = "./tmp"

    outname = "tmp"

    subprocess.call("mkdir tmp", shell=True)

    for i in range(1, 23):
        print("Applying PRScsx for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR = i

        print(
            """
            python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --pop=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s --out_name=%s
            """
            % (
                prscsx_path,
                ldref_dir,
                bim,
                ss,
                tmp_populations,
                tmp_N,
                str(CHR),
                str(phi),
                outdir,
                outname,
            ),
        )
        subprocess.call(
            """
            python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --pop=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s --out_name=%s
             """
            % (
                prscsx_path,
                ldref_dir,
                bim,
                ss,
                tmp_populations,
                tmp_N,
                str(CHR),
                str(phi),
                outdir,
                outname,
            ),
            shell=True,
        )

        for pop in populations:
            subprocess.call(
                "mv tmp/tmp_%s_pst*chr%s.txt tmp/tmp_%s_pst_chr%s.txt"
                % (pop, str(i), pop, str(i)),
                shell=True,
            )
        print("PRScsx for CHR " + str(i) + " is done!")

    print("Get adjusted_beta...")

    for processor, pop in zip(processors, populations):

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
        sumstats = pd.read_table(processor.sumstats)
        final_snps = list(set(df_adj_ss[1]) & set(sumstats["SNP"]))

        adjusted_ss = sumstats.copy()
        adjusted_ss = adjusted_ss[adjusted_ss.SNP.isin(final_snps)]

        adjusted_ss[use_col] = df_adj_ss[5].values

        save_path = processor.workdir + "/adjusted_sumstats_PRSCSx"
        adjusted_ss.to_csv(save_path, sep="\t", index=False)

        processor.adjusted_ss["PRSCSx"] = save_path

        processor.performance["PRSCSx"] = {}

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
