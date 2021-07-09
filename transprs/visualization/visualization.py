import seaborn as sns
import matplotlib as mpl
import pandas as pd


def visualize_performance(processor, metric, cmap="Dark2", plot_type="box_plot"):

    allow_type = ["box_plot", ""]

    assert plot_type in allow_type, "Please choose the right plot type: " + allow_type

    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False

    performance_df = pd.DataFrame.from_dict(processor.performance)
    selected_metric_df = pd.DataFrame(performance_df.loc[metric]).T

    if plot_type == "box_plot":
        sns.boxplot(data=selected_metric_df.apply(pd.Series.explode), palette=cmap)

        sns.stripplot(
            data=selected_metric_df.apply(pd.Series.explode), color="#333333", alpha=0.5
        )

        mpl.pyplot.xticks(rotation=45)
