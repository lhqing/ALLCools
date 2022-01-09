# DMR Caller Benchmark

## Introduction

Here we perform DMR calling on the pseudo-bulk methylomes. These methylomes are merged from single-cells based on cell
cluster labels and additional experimental variables (if any).

Unlike bulk-WGBS, these pseudo-bulk methylomes usually has large difference on the coverage, which impact the power and
FDR of DMR calling. We want to test how existing tools perform at different genome coverage levels using both simulated
and real single-cell methylome input.

## Tools to test

1. methylpy
2. DSS
3. dmrseq

## Pseudo-bulk methylome coverage:

Using DG as an example. Merge 10, 20, 50, 100, 200, 500, 1000 cells, plot

1. Percentage of C covered
2. ave. C coverage
3. compare to bulk WGBS from ENCODE

## Simulation

1. choose a deep pseudo-bulk profile as reference
2. perform reads level Monte Carlo simulation by
    1. sampling individual reads (120bp region) in the genome based on coverage;
    2. in each read, simulate the CpG base call based on the reference methylation level;
    3. aggregate all the reads to generate pseudo-bulk profile.
3. two kinds of simulated samples:
    1. null simulations: simulated from the same reference, used to estimate FDR.
    2. known DMR simulations: simulated from the same reference, but changing reference level (using different effect
       sizes) at known DMR locus. These samples can be used to estimate specificity and sensitivity.

## Pairwise DMR calling - cluster A vs cluster B

### Power analysis

Simulation: Use known DMR simulations at different coverage level. Smaller effect size will need larger coverage to
reach the same power.

### Saturation analysis

Real data: Merge pseudo-bulk methylomes from different number of cells, and identify DMRs from each pair, then calculate
the % of DMRs (segregate by different effect size) identified comparing to the full-coverage profile.

### False positives and FDR

Simulation: using null simulations at different coverage level, any DMR identified are false discoveries Real data:
using pseudo-bulk methylomes from different number of **shuffled** cells, calculate FDR comparing to cluster DMR at
matched coverage.

### Test imbalanced coverage

The above test can be performed on

1) Balanced coverage: cluster A and cluster B has the same coverage;
2) imbalanced coverage: cluster A and cluster B has very different coverage.

## Multi-group DMR calling

The same as pairwise DMR calling, but including 5 different clusters and use multi-group comparison in the same model.
