import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing

tax_order = list(reversed(["superkingdom", "phylum", "class","order", "family","genus","species"]))
tax_order_rank = {"superkingdom":6, "phylum":5, "class":4, "order":3, "family":2, "genus":1, "species":0}
tax_order_rank_r = {val:key for key,val in tax_order_rank.items()}

q_ranks_tid = pd.read_csv("./udance_rank_tid.tsv", sep="\t")
r_ranks_tid = pd.read_csv("./10k_ranks_tid.tsv", sep="\t")

r_ranks_tid = r_ranks_tid.rename(columns={"kingdom": "superkingdom"})
q_ranks_tid = q_ranks_tid.rename(columns={"kingdom": "superkingdom"})

q_ranks_tid.index = q_ranks_tid["genome"]
r_ranks_tid.index = r_ranks_tid["genome"]
r_ranks_tid = r_ranks_tid.drop("genome", axis=1)
q_ranks_tid = q_ranks_tid.drop("genome", axis=1)

consult_result_dir = Path("..results/consult-predictions/")

from collections import defaultdict

def calc(result_file):
    results_consult = defaultdict(list)
    genome = result_file.stem.split("_")[-1]
    true_tax = q_ranks_tid.loc[genome]
    with open(result_file) as f:
        lines = f.readlines()
    for line in (lines):
        read, pred = line.split(" ")
        read_name, form = read.split(":")
        tax_pred, w = pred.split(":")
        pred_v = {}
        if tax_pred != "NA":
            tax_pred, w = int(tax_pred), float(w)
            tax_lvl = r_ranks_tid.columns[np.max(np.nonzero((r_ranks_tid == tax_pred).any(axis=0).values)[0])]
            row = r_ranks_tid.iloc[np.nonzero((r_ranks_tid == tax_pred).any(axis=1).values)[0][0]]

            for lvl in tax_order[tax_order_rank[tax_lvl]:]:
                pred_v[lvl] = row[lvl]

        results_consult["genome"].append(genome)
        results_consult["read"].append(read_name)
        seen = False
        missed = False
        for idx, taxa in enumerate(tax_order):
            if true_tax[taxa] == 0:
                if (idx > 0) :
                    results_consult[taxa].append(results_consult[tax_order[idx-1]][-1])
                else:
                    results_consult[taxa].append("=")
            elif taxa in pred_v.keys():
                if pred_v[taxa] == true_tax[taxa]:
                    for i in range(idx, len(tax_order)):
                        results_consult[tax_order[i]].append("TP")
                    break
                else:
                    results_consult[taxa].append("FP")
                seen = True
            elif not seen:
                if (not missed) and (true_tax[taxa]!=r_ranks_tid[taxa]).all():
                    results_consult[taxa].append("TN")
                else:
                    results_consult[taxa].append("FN")
                    missed = True
            else:
                results_consult[taxa].append("?")
    return pd.DataFrame(results_consult)

pool = multiprocessing.Pool(96)
all_results = [pool.map(calc, [file for file in list(consult_result_dir.iterdir()) if not file.stem.startswith(".")])]
pd.concat(all_results[0]).to_csv("../results/consult_results.csv")
print("SAVED")
