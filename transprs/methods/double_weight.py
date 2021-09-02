import time
import datetime
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()


def double_weight(processor, top_choice=100000):

    start_time = time.time()

    print("Double weight method is running...")

    sumstats = pd.read_table(processor.sumstats)

    sumstats = pandas2ri.py2rpy(sumstats)

    robjects.globalenv["ss"] = sumstats
    robjects.globalenv["top_choice"] = top_choice

    robjects.r(
        """
        norm_samples <- apply(ss, 1, function(x) rnorm(100, mean = as.numeric(x[9]), sd = as.numeric(x[7])))

        rank_samples <- apply(norm_samples, 1, function(x) rank(x * -1))

        #calculate the number of time for each SNP the rank is less than the desired top number of SNPs
        prob_ranks <- apply(rank_samples, 1, function(x) sum(x < top_choice)/nrow(rank_samples))

        #multiple the normal beta values by the winner's curse probability to get an updated beta value
        ss[,9] <- ss[,9] * prob_ranks
        print("Adjusted BETA is done!")
    """
    )

    print("Done Double weight!")

    adjusted_ss = robjects.globalenv["ss"]

    save_path = processor.workdir + "/adjusted_sumstats_double_weight"
    adjusted_ss.to_csv(save_path, sep="\t", index=False)

    processor.adjusted_ss["double_weight"] = save_path

    processor.performance["double_weight"] = {}

    print("The double weight result stores in .adjusted_ss['double_weight']!")

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
