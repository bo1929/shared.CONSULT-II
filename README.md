# Shared files for CONSULT-II paper
## Public libraries
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (140 Gb - large but performant with defaults)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v020-WoL140G.tar.gz)
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (32 Gb - lighter-weight but still highly accurate)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v030-WoL32G.tar.gz)
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (18 Gb - lightweight and robust)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v030-WoL18G.tar.gz)
- [Taxonomy lookup tables for all WoL libraries](https://github.com/bo1929/shared.CONSULT-II/raw/master/misc/taxonomy_lookup-RefSeq2019.tar.gz)

## Queries used for read classification experiments
- Simulated short reads from bacterial genomes: [download](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/bacteria_queries_sampreads.tar.gz)
- Simulated short reads from archaeal genomes: [download](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/archaea_queries_sampreads.tar.gz)

## misc
- `./misc/ReferenceTaxonomy-nodes.dmp.gz`:
Full taxonomy we used.
- `./misc/10kBacteria-metadata.tsv`:
Metadata for reference genomes, including genome sizes that we used for profiling.
- `./misc/ReferenceTaxonomyRWoL-nodes.dmp.gz`
Reduced taxonomy, limited to species that are in the reference dataset, is used to create lookup tables so that 16bit is sufficient for taxon IDs.
- `./misc/dist-bacteria-to-closest.txt`:
Mash distances to the closest genomes (MinGND) for each bacterial query that we used, and the corresponding closest genome.
- `./misc/dist-archaea-to-closest.txt`:
Mash distances to the closest genomes (MinGND) for each archaeal query that we used, and the corresponding closest genome.
- `./misc/uDance-ranks_tid.tsv`:
Taxon ID mapping for reference genomes across different ranks (uDance dataset, includes bacterial queries).
- `./misc/10kBacteria-ranks_tid.tsv`:
Taxon ID mapping for reference genomes across different ranks (WoL dataset, includes references and archaeal queries).
- `./misc/taxonomy_lookup-WoLv1.tar.gz`:
Taxonomy lookup tables for library build using WoL-v1 dataset for read classification experiments. Needed for classification/profiling.
- `./misc/taxonomy_lookup-RefSeq2019.tar.gz`:
Taxonomy lookup tables for library build using RefSeq 2019 snapshot for CAMI-II challenge. Needed for classification/profiling.
## results
- `./results/profile-filtered_tv1000-norm_genome_sizes-marine.tar.gz`
CAMI-II marine dataset profiling results from all methods combined without genome size correction, OPAL output with `-n` option, `tv1000`     stands for the total vote threshold of 1000 for the entire sample.
- `./results/profile-filtered_tv1000-norm_free_baseline-strain_madness.tar.gz`
CAMI-II strain-madness dataset profiling results from all methods combined without genome size correction, OPAL output with `-n` option.
- `./results/profile-filtered_tv1000-norm_free_baseline-marine.tar.gz`
CAMI-II marine dataset profiling results from all methods combined with genome size correction, OPAL output with `-n` option.
- `./results/profile-filtered_tv1000-norm_genome_sizes-strain_madness.tar.gz`
CAMI-II strain-madness dataset profiling results from all methods combined with genome size correction, OPAL output with `-n` option.
- `./results/estimated_profiles-level_norm_tv1000.tar.gz`:
Results for CAMI-I profiling without genome size correction, and unclassified abundance is not allowed (CONSULT-II with `--force-unit` option), generated by OPAL using `-n` option.
- `./results/estimated_profiles-rootwgs_norm_tv1000.tar.gz`
Main results for CAMI CAMI-I high-complexity dataset, including different CONSULT-II approaches (with/without unclassified abundance) and other tools, generated by OPAL using `-n` option.
- `./results/estimated_profiles-levelwgs_norm_tv1000.tar.gz`
Results for CAMI-I profiling with genome size correction, and unclassified abundance is allowed, generated by OPAL with the `-n` option.
- `./results/estimated_profiles-root_norm_tv1000.tar.gz`
Results for CAMI-I profiling with genome size correction and unclassified abundance is allowed, generated by OPAL with the `-n` option.
- `./results/all_tools-profiling_evaluation-CAMI1_hc.tsv`:
CAMI-I profiling results all tools together for the figure.
- `./results/unclassified_abundances.tsv`
Unclassified abundance values across ranks for CAMI-I data.
- `./results/version_comparison-profiling_evaluation-CAMI1_hc.tsv`:
Comparison of profiling methods of CAMI-I results from RECOMB-CG paper and Bioinformatics paper, using Eq. 6 of the paper, without genome size correction and unclassified abundance.
- `./results/summary_scores_c.csv`
Comparison of different total vote thresholds ($0.0$, $0.01$, $0.03$) on bacterial query genomes for read classification, with MinGND bins.
- `./results/summary_scores_pu.csv`
Comparison of different (soft-)LCA methods ($w=2$, $w=4$, LCA) on bacterial queries for read classification, with distance to the closest (MinGND) bins.
- `./results/summary_p_scores_sizes-bacteria.csv`
Comparison of different library size of CONSULT-II (140Gb vs 32Gb) on bacterial queries for read classification, with distance to the closest (MinGND) values.
- `./results/summary_p_scores_methods-bacteria.csv`
All methods combined, classification results per bacterial query genome, with distance to the closest (MinGND) values.
- `./results/summary_scores_methods-bacteria.csv`
All methods combined, classification results per bacterial query genome, with bins for distance to the closest (MinGND).
- `./results/summary_p_scores_methods-archaea.csv`:
All classification results for archaeal queries, for Supplementary figure, with distance to the closest (MinGND) values.
- `./results/summary_scores_methods-archaea.csv`:
All classification results for archaeal queries, for Supplementary figure, with bins for distance to the closest (MinGND).
- `./results/resource_tool_comparison.tsv`:
Resource usages of all three tools for classification of smaller-sized query data.
- `./results/CONSULTII-evaluation-bacteria-w4s5_th05_c03.csv.xz`:
CONSULT-II classification results for $w=4$ and $s=5$, used in the empirical evaluation of soft LCA heuristic and final tool comparison, majority vote threshold (`th`) is $0.5$ and total vote threshold is $c=0.03)$.
- `./results/CONSULTII-evaluation-bacteria-h14_b10_c03.csv.xz`:
CONSULT-II classification results on bacterial queries for smaller library ($h=14$ and $b=10$) and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-evaluation-archaea-w4s5_th05_c03.csv.xz`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN values for each read with default values and the total vote threshold value of $0.03$.
- `./results/CONSULTII-evaluation-bacteria-w4s5_th05_c01.csv.xz`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN values for each read with default values and the total vote threshold value of $0.01$.
- `./results/CLARK-evaluation-archaea.csv.xz`
CLARK classification results on archaeal queries as TP/FP/TN/FN values for each read using default parameters.
- `./results/KrakenII-evaluation-archaea.csv.xz`
Kraken-II classification results on archaeal queries as TP/FP/TN/FN values for each read using default parameters.
- `./results/CONSULTII-evaluation-bacteria-h13_b16_c03.csv.xz`
CONSULT-II classification results on bacterial queries for smallest library ($h=13$ and $b=16$) and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-evaluation-bacteria-h14_b10_c00.csv.xz`
CONSULT-II classification results on bacterial queries for smaller library ($h=14$ and $b=10$) and the total vote threshold is $c=0.00$.
- `./results/CONSULTII-evaluation-archaea-w4s5_th05_c00.csv.xz`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN values for each read with default values and the total vote threshold value of $0.00$.
- `./results/CONSULTII-evaluation-bacteria-w4s5_th05_c00.csv.xz`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN values for each read with default values and the total vote threshold value of $0.00$.
- `./results/CONSULTII-evaluation-bacteria-wINFs5_th05_c00.csv.xz`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN values for each read with non-probabilistic LCA and the total vote threshold value of $0.00$.
- `./results/CONSULTII-evaluation-bacteria-w2s5_th05_c00.csv.xz`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN values for each read with $w=2$ for soft-LCA and the total vote threshold value of $0.00$.
- `./results/CLARK-evaluation-bacteria.csv.xz`
CLARK classification results on bacterial queries as TP/FP/TN/FN values for each read using default parameters.
- `./results/KrakenII-evaluation-bacteria.csv.xz`:
Kraken-II classification results on bacterial queries as TP/FP/TN/FN values for each read using default parameters.
- `./results/CONSULTII-summary_match_distances-bacteria-w4s5.tsv`
CONSULT-II k-mer match statistics (Hamming distance, matching rank, and counts) for each MinGND bin of bacterial queries using default parameter values.
- `./results/CONSULTII-summary_matches-bacteria-w4s5.csv`:
CONSULT-II k-mer match statistics (Hamming distance, matching rank, and counts) for each bacterial query genome bin using default parameter values.
- `./results/CONSULTII-version_comparison.txt`
Resource usage comparison with the original CONSULT and CONSULT-II using its heuristic to set parameters automatically.
- `./results/CLARK-summary_evaluation-archaea.tsv`
CLARK classification results on archaeal queries as TP/FP/TN/FN counts for each MinGND bin.
- `./results/CLARK-summary_evaluation-bacteria.tsv`
CLARK classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin.
- `./results/CLARK-summary_p_evaluation-archaea.tsv`
CLARK classification results on archaeal queries as TP/FP/TN/FN counts for each query genome with corresponding distance to the closest (MinGND) value.
- `./results/CLARK-summary_p_evaluation-bacteria.tsv`:
CLARK classification results on bacteria queries as TP/FP/TN/FN counts for each query genome with corresponding distance to the closest (MinGND) value.
- `./results/KrakenII-summary_p_evaluation-bacteria.tsv`
Kraken-II classification results on bacteria queries as TP/FP/TN/FN counts for each query genome with corresponding distance to the closest (MinGND) value.
- `./results/KrakenII-summary_evaluation-archaea.tsv`
Kraken-II classification results on archaeal queries as TP/FP/TN/FN counts for each MinGND bin.
- `./results/KrakenII-summary_p_evaluation-archaea.tsv`
Kraken-II classification results on archaeal queries as TP/FP/TN/FN counts for each query genome with corresponding distance to the closest (MinGND) value.
- `./results/KrakenII-summary_evaluation-bacteria.tsv`
Kraken-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin.

- `./results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c03.tsv`:
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using default parameters, and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c01.tsv`:
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using default parameters, and the total vote threshold is $c=0.01$.
- `./results/CONSULTII-summary_evaluation-bacteria-w4s5_th05_c00.tsv`:
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using default parameters, and the total vote threshold is $c=0.00$.
- `./results/CONSULTII-summary_evaluation-bacteria-wINFs5_th05_c00.tsv`:
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin using non-probabilistic LCA and the total vote threshold is $0.00$.
- `./results/CONSULTII-summary_evaluation-bacteria-w2s5_th05_c00.tsv`:
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin using $w=2$ for soft-LCA and the total vote threshold is $0.00$.
- `./results/CONSULTII-summary_evaluation-bacteria-h14_t2_l2_b10_c00.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using a smaller library ($h=14$ and $b=10$), the total vote threshold is $0.00$.
- `./results/CONSULTII-summary_evaluation-bacteria-h14_t2_l2_b10_c03.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using a smaller library ($h=14$ and $b=10$), the total vote threshold is $0.03$.
- `./results/CONSULTII-summary_evaluation-bacteria-h13_t2_l2_b16_c03.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each MinGND bin, using the light-weight library ($h=13$ and $b=16$), the total vote threshold is $0.03$.
- `./results/CONSULTII-summary_evaluation-archaea-w4s5_th05_c00.tsv`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN counts for each MinGND bin, using default parameters, and the total vote threshold is $c=0.00$.
- `./results/CONSULTII-summary_evaluation-archaea-w4s5_th05_c03.tsv`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN counts for each MinGND bin, using default parameters, and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-summary_p_evaluation-bacteria-w4s5_th05_c03.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each query genome, using default parameters, and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-summary_p_evaluation-bacteria-w4s5_th05_c00.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each query genome, using default parameters, and the total vote threshold is $c=0.00$.
- `./results/CONSULTII-summary_p_evaluation-bacteria-w4s5_th05_c01.tsv`
CONSULT-II classification results on bacterial queries as TP/FP/TN/FN counts for each query genome, using default parameters, and the total vote threshold is $c=0.01$.
- `./results/CONSULTII-summary_p_evaluation-archaea-w4s5_th05_c03.tsv`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN counts for each query genome, using default parameters, and the total vote threshold is $c=0.03$.
- `./results/CONSULTII-summary_p_evaluation-archaea-w4s5_th05_c00.tsv`
CONSULT-II classification results on archaeal queries as TP/FP/TN/FN counts for each query genome, using default parameters, and the total vote threshold is $c=0.00$.
## scripts
- `./scripts/match_evaluation_summary.R`
- `./scripts/evaluate_KrakenII.py`
- `./scripts/prepprocess_methods_psummary.py`
- `./scripts/profiling_size_correction.R`
- `./scripts/prepprocess_methods_summary.py`
- `./scripts/size_comparison.R`
- `./scripts/evaluate_CLARK.py`
- `./scripts/construct_taxonomy_lookup.py`
- `./scripts/calc_probabilities.R`
- `./scripts/prepprocess_c_summary.py`
- `./scripts/make_down_list.py`
- `./scripts/profiling_detailed_analysis.R`
- `./scripts/novelty_bins.R`
- `./scripts/postprocess_profiles.py`
- `./scripts/summarize_match_distances.py`
- `./scripts/profiling_tool_comparision.R`
- `./scripts/profiling_normalization_comparison.R`
- `./scripts/total_vote_impact.R`
- `./scripts/summarize_evaluations_KrakenII.py`
- `./scripts/profiling_unkown_taxa.R`
- `./scripts/classification_tool_comparison.R`
- `./scripts/profiling_cami2_analysis.R`
- `./scripts/evaluate_CONSULTII.py`
- `./scripts/prepprocess_pu_summary.py`
- `./scripts/pu_comparison_LCA.R`
- `./scripts/resource_tool_comparison.R`
- `./scripts/vote_threshold_comparison.R`
- `./scripts/summarize_evaluations_CLARK.py`
- `./scripts/summarize_evaluations_CONSULTII.py`
## figures
- `./figures/algorithm-illustration.pdf`
- `./figures/profiling-size_correction-CAMI2.pdf`
- `./figures/profiling-size_correction-CAMI1.pdf`
- `./figures/classification_tool_comparison-archaea-1.pdf`
- `./figures/resource_tool_comparison.pdf`
- `./figures/classification_size_comparison-bacteria-1.pdf`
- `./figures/novelty_bins-wolabels-bacteria.pdf`
- `./figures/classification_tool_comparison-archaea-2.pdf`
- `./figures/classification_size_comparison-bacteria-2.pdf`
- `./figures/profiling-marine-tool_comparison.pdf`
- `./figures/classification_tool_comparison-archaea-3.pdf`
- `./figures/novelty_bins-wlabels-archaea.pdf`
- `./figures/size_comparison.pdf`
- `./figures/profiling-legend-CAMI2.pdf`
- `./figures/total_nvote_impact.pdf`
- `./figures/novelty_bins-wlabels-bacteria.pdf`
- `./figures/algorithm-illustration.pages`
- `./figures/classification_tool_comparison-bacteria-1.pdf`
- `./figures/classification_tool_comparison-bacteria-3.pdf`
- `./figures/total_vote_impact.pdf`
- `./figures/expected_num_matches.pdf`
- `./figures/classification_tool_comparison-bacteria-2.pdf`
- `./figures/profiling_with_unclassified.pdf`
- `./figures/profiling_normalization_comparison-all_ranks_and_from_species.pdf`
- `./figures/profiling_normalization_comparison-all_ranks.pdf`
- `./figures/probability_LCA_update.pdf`
- `./figures/profiling_unclassified_comparison.pdf`
- `./figures/profiling_unclassified_taxon_proportions.pdf`
- `./figures/profiling-strain_madness-tool_comparison.pdf`
- `./figures/profiling-cami2_combined-tool_comparison.pdf`
- `./figures/vote_threshold_comparison-1.pdf`
- `./figures/pu_comparison_LCA.pdf`
- `./figures/profiling_tool_comparison.pdf`
- `./figures/vote_threshold_comparison-2.pdf`
- `./figures/number_of_matches.pdf`
