from transprs.utils import tmp_extract
import subprocess
import pandas as pd
import time
import datetime
import glob


def generate_prs(processor, method, use_sum=True):

    assert method in processor.adjusted_ss.keys(), (
        "Please run " + method + " before calculate PRS!"
    )
    start_time = time.time()
    print("Extracting adjusted sumstats from " + method + " method...")
    # tmp_extract(processor, method=method)
    # print("Done extract data!")

    effect_ind = "9"

    print("Generating PRS...")
    if use_sum:
        subprocess.call(
            """
        awk '{print $3,$8}' %s > tmp_SNP.pvalue

        echo "0.000000001 0 0.00000001" > tmp_range_list
        echo "0.00000001 0 0.0000001" > tmp_range_list
        echo "0.0000001 0 0.0000001" > tmp_range_list
        echo "0.000003 0 0.000003" > tmp_range_list
        echo "0.000001 0 0.000001" > tmp_range_list
        echo "0.00003 0 0.00003" > tmp_range_list
        echo "0.00001 0 0.00001" > tmp_range_list
        echo "0.0003 0 0.0003" > tmp_range_list
        echo "0.0001 0 0.0001" > tmp_range_list
        echo "0.003 0 0.003" > tmp_range_list
        echo "0.001 0 0.001" > tmp_range_list
        echo "0.05 0 0.05" >> tmp_range_list
        echo "0.1 0 0.1" >> tmp_range_list
        echo "0.2 0 0.2" >> tmp_range_list
        echo "0.3 0 0.3" >> tmp_range_list
        echo "0.4 0 0.4" >> tmp_range_list
        echo "0.5 0 0.5" >> tmp_range_list

        plink \
            --bfile %s \
            --score %s 3 4 %s header sum \
            --q-score-range tmp_range_list tmp_SNP.pvalue \
            --out tmp_results
        """
            % (
                processor.adjusted_ss[method],
                processor.population,
                processor.adjusted_ss[method],
                str(effect_ind),
            ),
            shell=True,
        )
    else:
        subprocess.call(
            """
        awk '{print $3,$8}' %s > tmp_SNP.pvalue

        echo "0.000000001 0 0.00000001" > tmp_range_list
        echo "0.00000001 0 0.0000001" > tmp_range_list
        echo "0.0000001 0 0.0000001" > tmp_range_list
        echo "0.000003 0 0.000003" > tmp_range_list
        echo "0.000001 0 0.000001" > tmp_range_list
        echo "0.00003 0 0.00003" > tmp_range_list
        echo "0.00001 0 0.00001" > tmp_range_list
        echo "0.0003 0 0.0003" > tmp_range_list
        echo "0.0001 0 0.0001" > tmp_range_list
        echo "0.003 0 0.003" > tmp_range_list
        echo "0.001 0 0.001" > tmp_range_list
        echo "0.05 0 0.05" >> tmp_range_list
        echo "0.1 0 0.1" >> tmp_range_list
        echo "0.2 0 0.2" >> tmp_range_list
        echo "0.3 0 0.3" >> tmp_range_list
        echo "0.4 0 0.4" >> tmp_range_list
        echo "0.5 0 0.5" >> tmp_range_list

        plink \
            --bfile %s \
            --score %s 3 4 %s header \
            --q-score-range tmp_range_list tmp_SNP.pvalue \
            --out tmp_results
        """
            % (
                processor.adjusted_ss[method],
                processor.population,
                processor.adjusted_ss[method],
                str(effect_ind),
            ),
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
