import pandas as pd
from pathlib import Path
from collections import defaultdict


PATH_RANKS = Path("./10kBacteria-ranks_tid.tsv")
DIR_LOOKUPS = Path("./taxonomy_lookup/")
DIR_LOOKUPS.mkdir(parents=True, exist_ok=True)

taxa_levels = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
taxa_level_ranks = {
    "kingdom": 7,
    "phylum": 6,
    "class": 5,
    "order": 4,
    "family": 3,
    "genus": 2,
    "species": 1,
}

if __name__ == "__main__":
    ranks_tid = pd.read_csv(PATH_RANKS, sep="\t")
    ancestor_lookup = defaultdict(list)
    genome_lookup = defaultdict(int)
    level_lookup = defaultdict(str)

    for i in range(ranks_tid.shape[0]):
        genome_lookup[ranks_tid.iloc[i]["genome"]] = ranks_tid.iloc[i]["species"]
        for ix, taxa1 in enumerate(taxa_levels[:]):
            tid1 = ranks_tid.iloc[i][taxa1]
            if tid1 not in ancestor_lookup.keys():
                for taxa2 in taxa_levels[ix:]:
                    ancestor_lookup[tid1].append(ranks_tid.iloc[i][taxa2])
                level_lookup[ranks_tid.iloc[i][taxa1]] = taxa1
            else:
                break

    with open(DIR_LOOKUPS / "ancestor_lookup", "w") as f:
        lookup_str = ""
        for taxa, ancestors in ancestor_lookup.items():
            lookup_str += f"{taxa} {','.join([str(a_taxa) for a_taxa in ancestors])}"
            lookup_str += "\n"
        f.write(lookup_str)

    with open(DIR_LOOKUPS / "genome_lookup", "w") as f:
        lookup_str = ""
        for key, val in genome_lookup.items():
            lookup_str += f"{key} {val}"
            lookup_str += "\n"
        f.write(lookup_str)

    with open(DIR_LOOKUPS / "level_lookup", "w") as f:
        lookup_str = ""
        for key, val in level_lookup.items():
            lookup_str += f"{key} {val}"
            lookup_str += "\n"
        f.write(lookup_str)
