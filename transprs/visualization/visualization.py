import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


def visualize_performance(
    processor,
    metric,
    validate=False,
    cmap="Dark2",
    plot_type="box_plot",
    figsize=(7, 7),
    fname=None,
    dpi=96,
):

    if validate:
        performance = processor.tuning
    else:
        performance = processor.performance

    sns.set(rc={"figure.figsize": figsize})

    allow_type = ["box_plot", "bar_plot"]

    assert plot_type in allow_type, "Please choose the right plot type: " + allow_type

    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False

    performance_df = pd.DataFrame.from_dict(performance)
    selected_metric_df = pd.DataFrame(performance_df.loc[metric]).T

    label = selected_metric_df.columns

    refined_label = list(map(lambda x: refine_labels(x), label))

    selected_metric_df.columns = refined_label

    if plot_type == "box_plot":
        sns.boxplot(data=selected_metric_df.apply(pd.Series.explode), palette=cmap)

        sns.stripplot(
            data=selected_metric_df.apply(pd.Series.explode), color="#333333", alpha=0.5
        )

        mpl.pyplot.xticks(rotation=45)
        plt.xticks(ha="right")

    if plot_type == "bar_plot":
        sns.barplot(data=selected_metric_df.apply(pd.Series.explode), palette=cmap)

        sns.stripplot(
            data=selected_metric_df.apply(pd.Series.explode), color="#333333", alpha=0.5
        )

        mpl.pyplot.xticks(rotation=45)
        plt.xticks(ha="right")

    if fname != None:
        mpl.pyplot.savefig(fname, bbox_inches="tight", pad_inches=0, dpi=dpi)


def refine_labels(label):
    if "+" in label:
        splitted = label.split("+")
        new_label = []
        for i, pop in enumerate(splitted):
            new_label.append(pop + "_" + str(i))
        new_label = "+".join(new_label)
        return new_label
    return label
