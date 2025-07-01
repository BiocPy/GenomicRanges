library(rtracklayer)
library(GenomicRanges)
library(microbenchmark)

bed_gr <- rtracklayer::import("consensus_peaks_bicnn.bed")
bigwig_gr <- rtracklayer::import("Astro.bw")

overlaps <- function(query, subject) {
  findOverlaps(query, subject)
}

results_overlaps <- microbenchmark(
  overlaps(bigwig_gr, bed_gr), times = 3
)
print(results_overlaps)

pretend_single_chrom_overlaps <- function(query, subject) {
  findOverlaps(ranges(query), ranges(subject))
}

results_pretend_single_chrom <- microbenchmark(
  pretend_single_chrom_overlaps(bigwig_gr, bed_gr), times = 3
)
print(results_pretend_single_chrom)
