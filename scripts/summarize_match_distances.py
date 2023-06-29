import sys
import collections as coll

BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.35)]


def get_bin(dist):
    bin_i = (0.35, 1)
    for bin_i in BINS:
        if (dist >= bin_i[0]) and (dist < bin_i[1]):
            break
    return bin_i

if __name__ == "__main__":
    results = coll.defaultdict(lambda: coll.defaultdict(lambda: coll.defaultdict(int)))

    with open(sys.argv[2], "r") as dist_file:
        dist_map = {genome[:-4]:dist for genome, dist in map(lambda x: (x[0], x[2]), map(lambda x: x.split(), dist_file.readlines()))}

    print(f"Bin\tRank\tCount\tHamming_distance\tMatch")

    with open(sys.argv[1], "r") as results_file:
        for line in results_file.readlines()[1: ]:
            ls = line.strip().split(",")
            dist_to_closest = float(dist_map[ls[1]])
            bin_genome = get_bin(dist_to_closest)
            print(f"{bin_genome}\t{ls[2]}\t{ls[3]}\t{ls[4]}\t{ls[5]}")
