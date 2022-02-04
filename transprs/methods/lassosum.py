import time
import datetime
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()


def lassosum(
    processor,
    reference,
    region_file,
    use_col,
    effect_ss,
    s=[0.1, 0.4, 0.8],
    lambda_value=[0.002, 0.004, 0.007, 0.1],
):

    start_time = time.time()

    print("Lassosum method is running...")

    sumstats = pd.read_table(processor.sumstats)

    sumstats = pandas2ri.py2rpy(sumstats)

    robjects.globalenv["ss"] = sumstats
    robjects.globalenv["reference"] = reference
    robjects.globalenv["region_file"] = region_file
    robjects.globalenv["s"] = s
    robjects.globalenv["lambda"] = lambda_value

    robjects.r(
        """
        library(lassosum)

        cor <- p2cor(p = ss$P, n = ss$N[1], sign=ss$BETA)

        out <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$BP,
                                 A1=ss$A1, A2=ss$A2, # A2 is not required but advised
                                 ref.bfile=reference, LDblocks = region_file, sample = 5000,
                                 s = s, lambda = lambda)

        new_ss <- ss[ss$BP %in% out$sumstats[,2],]
        print("Adjusted BETA is done!")
    """
    )

    print("Done lassosum!")

    adjusted_ss = robjects.globalenv["ss"]

    save_path = processor.workdir + "/adjusted_sumstats_double_weight"
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["lassosum"] = save_path

    processor.performance["lassosum"] = {}

    print("The double weight result stores in .adjusted_ss['lassosum']!")

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
