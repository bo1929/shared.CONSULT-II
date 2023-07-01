import sys
import collections as coll

BINS = [(0.0, 0.001), (0.001, 0.02), (0.02, 0.06), (0.06, 0.12), (0.12, 0.16), (0.16,0.22), (0.22, 0.35)]
CONT = True


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

    if CONT:
        with open(sys.argv[1], "r") as results_file:
            for line in results_file.readlines()[1: ]:
                ls = line.strip().split(",")[1:]
                results[ls[0]]["species"][ls[2]] +=1
                results[ls[0]]["genus"][ls[3]] +=1
                results[ls[0]]["family"][ls[4]] +=1
                results[ls[0]]["order"][ls[5]] +=1
                results[ls[0]]["class"][ls[6]] +=1
                results[ls[0]]["phylum"][ls[7]] +=1
                results[ls[0]]["superkingdom"][ls[8]] +=1

        for genome, results_bin in results.items():
            for key_rank, results_rank in results_bin.items():
                for key_type, count in results_rank.items():
                    dist_to_closest = float(dist_map[genome])
                    print(f"{genome}\t{dist_to_closest}\t{key_rank}\t{key_type}\t{count}")
    else:
        with open(sys.argv[1], "r") as results_file:
            for line in results_file.readlines()[1: ]:
                ls = line.strip().split(",")[1:]
                dist_to_closest = float(dist_map[ls[0]])
                bin_genome = get_bin(dist_to_closest)
                results[bin_genome]["species"][ls[2]] +=1
                results[bin_genome]["genus"][ls[3]] +=1
                results[bin_genome]["family"][ls[4]] +=1
                results[bin_genome]["order"][ls[5]] +=1
                results[bin_genome]["class"][ls[6]] +=1
                results[bin_genome]["phylum"][ls[7]] +=1
                results[bin_genome]["superkingdom"][ls[8]] +=1

        for key_bin, results_bin in results.items():
            for key_rank, results_rank in results_bin.items():
                for key_type, count in results_rank.items():
                    print(f"{key_bin}\t{key_rank}\t{key_type}\t{count}")
