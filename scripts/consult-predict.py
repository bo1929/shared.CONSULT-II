import json
from pathlib import Path

tax_order = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
results = Path("./results-CONSULT/")

for path in results.iterdir():
    with open(path) as f:
        d = json.load(f)
    genome = path.stem

    prev_pred = (None, None, None, -1, -1)
    for idx, (read, votes) in enumerate(d.items()):
        read = read.strip()
        with open(f"./CONSULT-predictions-07/{genome}.txt", 'a') as f:
            read_name, form = read.split(":")
            if len(votes)> 0:
                max_vote = max(votes["kingdom"].values())
                th = 0.5*max_vote/10
                if th > 0.01:
                    for lvl in tax_order:
                        for pred_name, val in votes.get(lvl, {}).items():
                            if val > th:
                                pred = (read_name, form, pred_name, val, th)
                                break
                        if pred[3] > th:
                            break
                else:
                    pred = (read_name, form, "NA", 0, 0)
            else:
                pred = (read_name, form, "NA", 0, 0)

            if idx%2  == 1:
                if pred[4] < prev_pred[4]:
                    pred = prev_pred
                f.write(f"{pred[0]}:{pred[1]} {pred[2]}:{pred[3]}\n")
                # print(f"{pred[0]}:{pred[1]} {pred[2]}:{pred[3]}\n")
                prev_pred = (read, None, None, -1, -1)
            else:
                prev_pred = pred
