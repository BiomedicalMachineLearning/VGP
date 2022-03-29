from subprocess import Popen, PIPE, STDOUT, CalledProcessError
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

    # Set time
    start_time = time.time()
    print("PRScsx is running...")

    # Check number of data object
    n_processors = len(processors)

    # Get number of chr
    CHR = pd.read_table(processors[0].sumstats).CHR.unique().astype(str)
    CHR_input = ",".join(CHR)

    # Split genotype to each chr
    process_split = Popen(
        """
        for chr in %s; do \
            plink --bim %s.bim --chr $chr --make-just-bim --out tmp${chr}; \
        done
        """
        % (CHR_input, processors[0].validation),
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

    # Track method source code
    path = os.path.dirname(transprs.__file__)
    prscsx_path = path + "/methods/multipop/prscsx/PRScsx.py"

    # Store temporary sumstats for each data object to prepare input
    # for PRSCSx
    ss = ""
    for i in range(0, len(processors)):
        tmp = pd.read_table(processors[i].sumstats)
        tmp[["SNP", "A1", "A2", "BETA", "P"]].to_csv(
            "tmp_sumstat" + str(i), sep=" ", index=False
        )
        ss = ss + ",tmp_sumstat" + str(i)
    ss = ss[1:]

    # Prepare population names
    tmp_populations = ",".join(populations)

    # Prepare N samples
    tmp_N = ",".join(map(str, N))

    # Prepare output folder
    outdir = "./tmp"

    outname = "tmp"

    Popen("mkdir tmp", shell=True)

    # Apply PRSCSx script for each chr
    for i in CHR:
        print("Applying PRScsx for CHR " + str(i) + "...")
        bim = "./tmp" + str(i)
        CHR_no = i

        # Main run
        process_run = Popen(
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
                str(CHR_no),
                str(phi),
                outdir,
                outname,
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
        for pop in populations:
            process_move = Popen(
                "mv tmp/tmp_%s_pst*chr%s.txt tmp/tmp_%s_pst_chr%s.txt"
                % (pop, str(i), pop, str(i)),
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

        print("PRScsx for CHR " + str(i) + " is done!")

    print("Get adjusted_beta...")

    for processor, pop in zip(processors, populations):

        df_adj_ss = pd.DataFrame(columns=[0, 1, 2, 3, 4, 5])

        for i in CHR:
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

        df_adj_ss.columns = ["0", "SNP", "2", "3", "4", "5"]
        adjusted_ss = pd.merge(adjusted_ss, df_adj_ss)[
            ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "5"]
        ]
        adjusted_ss.columns = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA"]

        save_path = processor.workdir + "/adjusted_sumstats_PRScsx"
        adjusted_ss.to_csv(save_path, sep="\t", index=False)

        processor.adjusted_ss["PRScsx"] = save_path

        processor.tuning["PRScsx"] = {}

        processor.performance["PRScsx"] = {}

    print("The clumping result stores in .adjusted_ss['PRScsx']!")

    subprocess.call(
        """
    rm ./tmp*
    rm ./tmp/*
        """,
        shell=True,
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
