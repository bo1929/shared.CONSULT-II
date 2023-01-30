import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
from pathlib import Path
import json
from collections import defaultdict
import multiprocessing
import joblib

tax_order = ["kingdom", "phylum", "class","order", "family","genus","species"]
tax_order_rank = {"kingdom":6, "phylum":5, "class":4, "order":3, "family":2, "genus":1, "species":0}
tax_order_rank_r = {val:key for key,val in tax_order_rank.items()}

q_ranks_tid = pd.read_csv("./udance_rank_tid.tsv", sep="\t")
r_ranks_tid = pd.read_csv("./10k_ranks_tid.tsv", sep="\t")

consult_results_dir = Path("./results/bac_output_CONSULT_results/p_results/")

T=0.51

def calc_votes(result_path, r_ranks_tid, q_ranks_tid, tax_order_rank, tax_order_rank_r):
    all_votes = {}
    with open(result_path, 'r') as reader:
        lines = reader.readlines()
    genome = result_path.stem.split("_")[-1].strip()
    name = None
    for idx, line in enumerate(lines[:5]):
        if idx%3 != 0:
            votes = defaultdict(dict)
            each_vote = defaultdict(dict)
            for match in line.split(" ")[1:]:
                taxid, dist = list(map(lambda x: (x.strip()), match.split(":")))
                taxid, dist = int(taxid), int(dist)
                if taxid >=0:
                    ref_row_idx = np.nonzero((r_ranks_tid == taxid).any(axis=1).to_numpy())[0][0]
                    tax_order_s = r_ranks_tid.iloc[ref_row_idx].index[(r_ranks_tid == taxid).any()][0]
                    if taxid == 0:
                        tax_order_s == "species"
                    for order in range(tax_order_rank[tax_order_s], len(tax_order_rank)):
                        each_vote[tax_order_rank_r[order]][str(r_ranks_tid.at[ref_row_idx, tax_order_rank_r[order]])] = each_vote[tax_order_rank_r[order]].get(str(r_ranks_tid.at[ref_row_idx, tax_order_rank_r[order]]), 0) + (1-dist/32)**32
            all_votes[name + ":" + line.split(" ")[0]] = each_vote
        else:
            name = line.strip()
    return all_votes

def calc_votes_(x):
    return calc_votes(x, r_ranks_tid, q_ranks_tid, tax_order_rank, tax_order_rank_r)

pool = multiprocessing.Pool(120)
all_results = pool.map(calc_votes_, [x for x in consult_results_dir.iterdir() if not str(x.stem).startswith(".")])

for idx, path in enumerate([x for x in consult_results_dir.iterdir() if not str(x.stem).startswith(".")]):
    genome = path.stem.split("_")[-1]
    with open(f'results-CONSULT/{genome}.json', 'w') as f:
        json.dump(all_results[idx], f)
