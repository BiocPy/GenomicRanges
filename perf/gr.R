library(rtracklayer)
library(GenomicRanges)
library(microbenchmark)

query <- rtracklayer::import("consensus_peaks_bicnn.bed")
subject <- rtracklayer::import("Astro.bw")

results_overlaps <- microbenchmark(
  findOverlaps(query, subject), times = 3
)
print(results_overlaps)

results_pretend_single_chrom <- microbenchmark(
  findOverlaps(ranges(query), ranges(subject)), times = 3
)
print(results_pretend_single_chrom)

results_nearest <- microbenchmark(
  nearest(query, subject), times = 3
)
print(results_nearest)

results_nearest_pretend_single_chrom <- microbenchmark(
  nearest(ranges(query), ranges(subject)), times = 3
)
print(results_nearest_pretend_single_chrom)
