from subprocess import Popen, PIPE, STDOUT, CalledProcessError
import pandas as pd
import time
import datetime
import os
import transprs
import pandas as pd


def polyfun(
    processor,
    N=None,
    use_col="BETA",
    ref_dir=None,
    max_num_causal=2,
):
    start_time = time.time()
    print("PolyFun is running...")

    # Track method source code
    path = os.path.dirname(transprs.__file__)
    polypred_path = path + "/methods/polypred"

    poly_dir = processor.workdir + "/polyfun/"

    sumstat_file = poly_dir + "sumstats.munged.parquet"
    N = pd.read_table(processor.sumstats).N[0]

    try:
        os.mkdir(poly_dir)
    except:
        print(poly_dir + "is already created")
        pass

    process_munge = Popen(
        """
        python %s/munge_polyfun_sumstats.py \
            --sumstats %s \
            --out %s
    """
        % (polypred_path, processor.sumstats, sumstat_file),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )
    with process_munge.stdout:
        try:
            for line in iter(process_munge.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")

    process_finemap = Popen(
        """
            python %s/create_finemapper_jobs.py \
                --sumstats %s  \
                --n %s \
                --method susie \
                --non-funct \
                --max-num-causal %s \
                --allow-missing \
                --jobs-file %s/jobs.txt \
                --out-prefix %s/polyfun_output
        """
        % (
            polypred_path,
            sumstat_file,
            str(N),
            str(max_num_causal),
            poly_dir,
            poly_dir,
        ),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )

    with process_finemap.stdout:
        try:
            for line in iter(process_finemap.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")

    if ref_dir != None:
        tmp = pd.read_table(poly_dir + "/jobs.txt", header=None).values
        jobs = list(
            map(
                lambda x: x[0].replace(
                    "https://data.broadinstitute.org/alkesgroup/UKBB_LD/", ref_dir
                ),
                tmp,
            )
        )

        with open(poly_dir + "/jobs.txt", "w") as f:
            for item in jobs:
                f.write("%s\n" % item)

    print("Finemapping is running. It takes a while...")
    # process_run_finemap = Popen(
    #     """
    #     bash %s/jobs.txt
    # """
    #     % (poly_dir),
    #     shell=True,
    #     stdout=PIPE,
    #     stderr=STDOUT,
    # )

    # with process_run_finemap.stdout:
    #     try:
    #         for line in iter(process_run_finemap.stdout.readline, b""):
    #             print(line.decode("utf-8").strip())

    #     except CalledProcessError as e:
    #         print(f"{str(e)}")

    process_aggregate = Popen(
        """
            python %s/aggregate_finemapper_results.py \
                --out-prefix %s/polyfun_output \
                --sumstats %s \
                --out %s/polyfun_output.agg.txt.gz \
                --allow-missing-jobs
        """
        % (polypred_path, poly_dir, sumstat_file, poly_dir),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )

    with process_aggregate.stdout:
        try:
            for line in iter(process_aggregate.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")

    print("Get adjusted_beta...")
    tmp = pd.read_table(poly_dir + "/polyfun_output.agg.txt.gz")
    ss = pd.read_table(processor.sumstats)
    adjusted_ss = pd.merge(tmp, ss, on="SNP")[
        ["CHR_x", "BP_x", "SNP", "A1_x", "A2_x", "N_x", "SE", "P_y", "BETA_MEAN"]
    ]
    adjusted_ss.columns = ss.columns

    save_path = processor.workdir + "/adjusted_sumstats_polyfun"
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["polyfun"] = save_path

    processor.performance["polyfun"] = {}

    print("The polyfun result stores in .adjusted_ss['polyfun']!")

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
