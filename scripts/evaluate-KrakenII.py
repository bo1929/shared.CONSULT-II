import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing
from collections import defaultdict

NUM_THREADS = 96
taxa_order = list(
    ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
)
taxa_order_rank = {
    "superkingdom": 6,
    "phylum": 5,
    "class": 4,
    "order": 3,
    "family": 2,
    "genus": 1,
    "species": 0,
}

def evaluate_classification(result_file):
    evaluation_dict = defaultdict(list)

    # example filename : output_1000x_G000251165.txt
    genome = result_file.stem.split("_")[-1]
    true_taxa = query_ranks.loc[genome]

    with open(result_file) as f:
        all_reads = f.readlines()
    for read_line in all_reads:
        read_name = read_line.split("\t")[1]
        is_classified = read_line.split("\t")[0]
        pred_tID = int(read_line.split("\t")[2])

        pred_taxa = {}
        rank = None

        if is_classified != "U":
            parent, rank = (
                reference_taxonomy.at[pred_tID, 1],
                reference_taxonomy.at[pred_tID, 2],
            )
            pred_taxa[rank] = pred_tID

            while parent != 1:
                pred_tID = parent
                parent, rank = (
                    reference_taxonomy.at[parent, 1],
                    reference_taxonomy.at[parent, 2],
                )
                pred_taxa[rank] = pred_tID

        evaluation_dict["genome"].append(genome)
        evaluation_dict["read"].append(read_name)

        seen = False
        miss = False

        for idx, taxon in enumerate(taxa_order):
            if true_taxa[taxon] == 0:
                if taxon != "species":
                    evaluation_dict[taxon].append(evaluation_dict[taxa_order[idx - 1]][-1])
                else:
                    raise ValueError("Taxonomic ID is 0 at species level.")
            elif taxon in pred_taxa.keys():
                if pred_taxa[taxon] == true_taxa[taxon]:
                    for i in range(idx, len(taxa_order)):
                        evaluation_dict[taxa_order[i]].append("TP")
                    break
                else:
                    evaluation_dict[taxon].append("FP")
                seen = True
            elif not seen:
                if (not miss) and (true_taxa[taxon] != reference_ranks[taxon]).all():
                    evaluation_dict[taxon].append("TN")
                else:
                    evaluation_dict[taxon].append("FN")
                    miss = True
            else:
                if taxon != "species":
                    evaluation_dict[taxon].append(evaluation_dict[taxa_order[idx - 1]][-1])
                else:
                    raise ValueError("Some unkown condition or conflicting taxonomic ID.")
    return pd.DataFrame(evaluation_dict)


if __name__ == "__main__":
    query_ranks = pd.read_csv("../misc/uDance-ranks_tid.tsv", sep="\t")
    reference_ranks = pd.read_csv("../misc/10kBacteria-ranks_tid.tsv", sep="\t")
    reference_taxonomy = pd.read_csv(
        "../misc/ReferenceTaxonomy-nodes.dmp", sep="|", header=None
    )
    kraken_result_dir = Path("../results/KrakenII-bacteria/")
    # Note that we have KrakenII output files with the following with 
    # filenames as "output_1000x_<GENOME_NAME>.txt".

    # Make sure we use superkingdom consistently.
    reference_ranks = reference_ranks.rename(columns={"kingdom": "superkingdom"})
    query_ranks = query_ranks.rename(columns={"kingdom": "superkingdom"})

    # Clean whitespaces and newlines etc.
    reference_taxonomy[2] = reference_taxonomy[2].apply(lambda x: x.strip())

    # Set DataFrame indices to genome names.
    query_ranks.index = query_ranks["genome"]
    reference_ranks.index = reference_ranks["genome"]
    reference_ranks = reference_ranks.drop("genome", axis=1)
    query_ranks = query_ranks.drop("genome", axis=1)

    # Set DataFrame indices to taxonomic IDs.
    reference_taxonomy.index = reference_taxonomy[0]

    pool = multiprocessing.Pool(NUM_THREADS)
    all_evaluations = [
        pool.map(
            evaluate_classification,
            [
                file
                for file in list(kraken_result_dir.iterdir())
                if not file.stem.startswith(".")
            ],
        )
    ]
    pd.concat(all_evaluations[0]).to_csv("../results/KrakenII-bacteria-eval.csv")
