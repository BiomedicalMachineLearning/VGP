from transprs.utils import tmp_extract
import subprocess
import pandas as pd
import time
import datetime
import glob


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
    awk '{print $3,$8}' tmp_ss > tmp_SNP.pvalue

    echo "0.001 0 0.001" > tmp_range_list
    echo "0.05 0 0.05" >> tmp_range_list
    echo "0.1 0 0.1" >> tmp_range_list
    echo "0.2 0 0.2" >> tmp_range_list
    echo "0.3 0 0.3" >> tmp_range_list
    echo "0.4 0 0.4" >> tmp_range_list
    echo "0.5 0 0.5" >> tmp_range_list

    plink \
        --bfile tmp \
        --score tmp_ss 3 4 %s header sum \
        --q-score-range tmp_range_list tmp_SNP.pvalue \
        --out tmp_results
    """
        % (str(effect_ind)),
        shell=True,
    )

    print("PRS is generated!")

    processor.prs_results[method] = {}
    snp_results_files = glob.glob("./*.profile")
    import re

    for filename in snp_results_files:
        prs_result = pd.read_csv(filename, delim_whitespace=True)
        new_key = re.sub(".profile", "", filename[14:])
        processor.prs_results[method][new_key] = prs_result

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
