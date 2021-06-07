from transprs.utils import tmp_extract
import subprocess
import pandas as pd
import time
import datetime


def generate_prs(processor, use_col, method):

    assert method in processor.adjusted_ss.keys(), (
        "Please run " + method + " before calculate PRS!"
    )
    start_time = time.time()
    print("Extracting adjusted sumstats from " + method + " method...")
    tmp_extract(processor, method=method)
    print("Done extract data!")

    if use_col == "OR":
        effect_ind = "9"
    else:
        effect_ind = "12"

    print("Generating PRS...")

    subprocess.call(
        """
    plink \
        --bfile tmp \
        --score tmp_ss 3 4 %s header \
        --out tmp_results
    """
        % (str(effect_ind)),
        shell=True,
    )

    print("PRS is generated!")

    prs_result = pd.read_csv("tmp_results.profile", delim_whitespace=True)

    processor.prs_results[method] = prs_result

    print("The PRS result stores in .prs_results['" + method + "']!")

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
