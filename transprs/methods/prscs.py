from subprocess import Popen, PIPE, STDOUT, CalledProcessError
import subprocess
import pandas as pd
import time
import datetime
import os
import transprs


def prscs(
    processor,
    ldref_dir,
    N,
    use_col="BETA",
    phi=1e-2,
):
    start_time = time.time()
    
    print("PRScs is running...")

    # Get number of chr
    CHR = pd.read_table(processor.sumstats).CHR.unique().astype(str)
    CHR_input = " ".join(CHR)

    # Split genotype to each chr
    process_split = Popen(
        """
        for chr in %s; do \
            plink --bim %s.bim --chr $chr --make-just-bim --out tmp${chr}; \
        done
        """
        % (CHR_input, processor.validation),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )

    with process_split.stdout:
        try:
            for line in iter(process_split.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")


    path = os.path.dirname(transprs.__file__)
    prscs_path = path + "/methods/prscs/PRScs.py"

    # Store temporary sumstats for each data object to prepare input
    # for PRSCS

    ss = pd.read_table(processor.sumstats)
    ss[["SNP", "A1", "A2", use_col, "P"]].to_csv("tmp_ss", sep="\t", index=False)

    # Prepare output folder
    outdir = "./tmp"

    outname = "tmp"

    Popen("mkdir tmp", shell=True)

    # Apply PRSCS script for each chr
    for i in CHR:
        print("Applying PRScs for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR_no = i

        # Main run
        process_run = Popen(
            """
            python %s --ref_dir=%s --bim_prefix=%s --sst_file=%s --n_gwas=%s --chrom=%s --phi=%s --out_dir=%s
            """
            % (
                prscs_path,
                ldref_dir,
                bim,
                "tmp_ss",
                str(N),
                str(CHR_no),
                str(phi),
                outdir+outname,
            ),
            shell=True,
            stdout=PIPE,
            stderr=STDOUT,
        )

        with process_run.stdout:
            try:
                for line in iter(process_run.stdout.readline, b""):
                    print(line.decode("utf-8").strip())

            except CalledProcessError as e:
                print(f"{str(e)}")

        # Moving the result to tmp output dir
        process_move = Popen(
            "mv tmp/tmp_pst*chr%s.txt tmp/tmp_pst_chr%s.txt"
            % (str(i),str(i)),
            shell=True,
            stdout=PIPE,
            stderr=STDOUT,
        )
        with process_move.stdout:
            try:
                for line in iter(process_move.stdout.readline, b""):
                    print(line.decode("utf-8").strip())

            except CalledProcessError as e:
                print(f"{str(e)}")

        print("PRScs for CHR " + str(i) + " is done!")


    print("Get adjusted_beta...")

    df_adj_ss = pd.DataFrame(columns=[0, 1, 2, 3, 4, 5])
    for i in CHR:
        try:
            df_next = pd.read_table(
                "tmp/tmp_pst_chr" + str(i) + ".txt", header=None
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

    df_adj_ss.columns = ["0", "SNP", "2", "3", "4", "5"]
    adjusted_ss = pd.merge(adjusted_ss, df_adj_ss)[
        ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "5"]
    ]
    adjusted_ss.columns = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA"]
    save_path = processor.workdir + "/adjusted_sumstats_PRSCS"
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["PRScs"] = save_path

    processor.tuning["PRScs"] = {}

    processor.performance["PRScs"] = {}

    print("The PRScs result stores in .adjusted_ss['PRScs']!")

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
