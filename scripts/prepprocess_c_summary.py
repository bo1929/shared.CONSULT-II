import pandas as pd
from collections import defaultdict

c00 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c00.tsv", sep="\t", header=None)
c01 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c01.tsv", sep="\t", header=None)
c03 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c03.tsv", sep="\t", header=None)


def compute_metrics(df, name=""):
    scores = defaultdict(list)
    for bin_rank, sub_df in df.groupby([0, 1]):
        values_TP = sub_df.loc[sub_df[2] == "TP"][3].values
        values_TN = sub_df.loc[sub_df[2] == "TN"][3].values
        values_FP = sub_df.loc[sub_df[2] == "FP"][3].values
        values_FN = sub_df.loc[sub_df[2] == "FN"][3].values
        if values_TP.size > 0:
            count_TP = values_TP[0]
        else:
            count_TP = 0
        if values_TN.size > 0:
            count_TN = values_TN[0]
        else:
            count_TN = 0
        if values_FP.size > 0:
            count_FP = values_FP[0]
        else:
            count_FP = 0
        if values_FN.size > 0:
            count_FN = values_FN[0]
        else:
            count_FN = 0
        precision = (count_TP) / (count_TP+count_FP+10**-5)
        recall = (count_TP) / (count_TP+count_FN+10**-5)
        scores["Taxonomic_Rank"].append(bin_rank[1])
        scores["Distance_to_closest"].append(f"[{eval(bin_rank[0])[0]}, {eval(bin_rank[0])[1]})")
        scores["Precision"].append(precision)
        scores["Recall"].append(recall)
        scores["F1"].append((2*precision*recall)/(precision+recall+10**-5))
        scores["Method"].append(name)
    return pd.DataFrame(scores)
scores_00 = compute_metrics(c00, "(0.00)")
scores_01 = compute_metrics(c01, "(0.01)")
scores_03 = compute_metrics(c03, "(0.03)")
scores = pd.concat((scores_00, scores_01, scores_03))

scores.to_csv("../results/summary_scores_c.csv")
