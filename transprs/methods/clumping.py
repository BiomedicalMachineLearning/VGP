from transprs.utils import tmp_extract
import subprocess
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
import time
import datetime
import pandas as pd


def clumping(
    processor,
    clump_p1=1,
    clump_r2=0.5,
    clump_kb=250,
):

    start_time = time.time()
    print("Clumping is running...")
    process_clumping = Popen(
        """
        plink \
            --bfile %s \
            --clump-p1 %s \
            --clump-r2 %s \
            --clump-kb %s \
            --clump %s \
            --clump-snp-field SNP \
            --clump-field P \
            --out tmp_out


        awk 'NR!=1{print $3}' tmp_out.clumped >  tmp_out.valid.snp
        """
        % (
            processor.population,
            str(clump_p1),
            str(clump_r2),
            str(clump_kb),
            processor.sumstats,
        ),
        shell=True,
        stdout=PIPE,
        stderr=STDOUT,
    )

    with process_clumping.stdout:
        try:
            for line in iter(process_clumping.stdout.readline, b""):
                print(line.decode("utf-8").strip())

        except CalledProcessError as e:
            print(f"{str(e)}")

    print("Done clumping!")
    valid_snps = list(pd.read_table("tmp_out.valid.snp", header=None)[0])

    sumstats = pd.read_table(processor.sumstats)

    adjusted_ss = sumstats[sumstats.SNP.isin(valid_snps)]

    save_path = processor.workdir + "/adjusted_sumstats_clumping"
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["clumping"] = save_path

    processor.performance["clumping"] = {}

    print("The clumping result stores in .adjusted_ss['clumping']!")

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
