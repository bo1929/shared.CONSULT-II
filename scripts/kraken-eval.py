import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing

tax_order = list(reversed(["superkingdom", "phylum", "class","order", "family","genus","species"]))
tax_order_rank = {"superkingdom":6, "phylum":5, "class":4, "order":3, "family":2, "genus":1, "species":0}
tax_order_rank_r = {val:key for key,val in tax_order_rank.items()}

q_ranks_tid = pd.read_csv("../data/udance_rank_tid.tsv", sep="\t")
r_ranks_tid = pd.read_csv("../data/10k_ranks_tid.tsv", sep="\t")
taxdump = pd.read_csv("../data/nodes.dmp", sep="|", header=None)

r_ranks_tid = r_ranks_tid.rename(columns={"kingdom": "superkingdom"})
q_ranks_tid = q_ranks_tid.rename(columns={"kingdom": "superkingdom"})

taxdump[2] = taxdump[2].apply(lambda x: x.strip())
taxdump.index = taxdump[0] 

q_ranks_tid.index = q_ranks_tid["genome"]
r_ranks_tid.index = r_ranks_tid["genome"]
r_ranks_tid = r_ranks_tid.drop("genome", axis=1)
q_ranks_tid = q_ranks_tid.drop("genome", axis=1)

kraken_result_dir = Path("../results/bac_output_kraken_results/")

from collections import defaultdict

def calc(result_file):
    results_kraken = defaultdict(list)
    genome = result_file.stem.split("_")[-1]
    true_tax = q_ranks_tid.loc[genome]
    with open(result_file) as f:
        lines = f.readlines()
    for line in (lines):
        c, tax_pred = line.split("\t")[0], int(line.split("\t")[2])
        pred_v = {}
        lvl = None
        if c != "U":
            parent, lvl = taxdump.at[tax_pred, 1], taxdump.at[tax_pred, 2]
            pred_v[lvl] = tax_pred
            while (parent != 1):
                tax_pred = parent
                parent, lvl = taxdump.at[parent, 1], taxdump.at[parent, 2]
                pred_v[lvl] = tax_pred

        results_kraken["genome"].append(genome)
        results_kraken["read"].append(line.split("\t")[1])
        seen = False
        missed = False
        for idx, taxa in enumerate(tax_order):
            if true_tax[taxa] == 0:
                if (idx > 0) :
                    results_kraken[taxa].append(results_kraken[tax_order[idx-1]][-1])
                else:
                    results_kraken[taxa].append("=")
            elif taxa in pred_v.keys():
                if pred_v[taxa] == true_tax[taxa]:
                    for i in range(idx, len(tax_order)):
                        results_kraken[tax_order[i]].append("TP")
                    break
                else:
                    results_kraken[taxa].append("FP")
                seen = True
            elif not seen:
                if (not missed) and (true_tax[taxa]!=r_ranks_tid[taxa]).all():
                    results_kraken[taxa].append("TN")
                else:
                    results_kraken[taxa].append("FN")
                    missed = True
            else:
                results_kraken[taxa].append("?")
    return pd.DataFrame(results_kraken)

pool = multiprocessing.Pool(96)
all_results = [pool.map(calc, [file for file in list(kraken_result_dir.iterdir()) if not file.stem.startswith(".")])]
pd.concat(all_results[0]).to_csv("../results/kraken_results.csv")
print("SAVED")
