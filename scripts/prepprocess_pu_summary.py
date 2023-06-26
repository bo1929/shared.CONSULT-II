import pandas as pd
from collections import defaultdict

w2s5 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-w2s5_th05_c00.tsv", sep="\t", header=None)
w4s5 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c00.tsv", sep="\t", header=None)
wINFs5 = pd.read_csv("../results/CONSULTII-summary_evaluation-bacteria-wINFs5_th05_c00.tsv", sep="\t", header=None)


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
scores_wINFs5 = compute_metrics(wINFs5, "non-probabilistic")
scores_w4s5 = compute_metrics(w4s5, "w=4, s=5")
scores_w2s5 = compute_metrics(w2s5, "w=2, s=5")
scores = pd.concat((scores_w4s5, scores_w2s5, scores_wINFs5))

scores.to_csv("../results/summary_scores_pu.csv")
