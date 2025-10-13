# GBWT-Query

GBWT-Query packages the GBWT-compatible indexes and queries needed to detect long shared haplotype segments directly on pangenome graphs. The primary entry point is `ibdCaller`, a CLI that enumerates candidate identity-by-descent (IBD) segments between all haplotypes stored in a GBWT.

## Overview
- `ibdCaller` scans each forward-orientation haplotype, computes long shared matches, and writes a TSV report (`hap1_id  hap2_id  hap1_start  hap2_start  length_nodes`).
- Reusable headers in `include/` expose the supporting structures (`FastLCP`, `lf_gbwt::GBWT`, `CompText`, and query helpers) for integration in other tools.
- Toy panels and pre-built indexes under `toyData/` enable quick experimentation.

## Quick Start
```bash
# install sdsl-lite and gbwt one directory above the repos
git clone https://github.com/vgteam/sdsl-lite
git clone https://github.com/jltsiren/gbwt
cd sdsl-lite && ./install.sh .. && cd ../gbwt && ./install.sh ..

# build and run the caller on a toy GBWT
cd ../GBWT-Query/test
make ibdCaller
./ibdCaller ../toyData/indexes/bi_fig1/bi_fig1 4 --version 4 --output bi_fig1_ibd.tsv --verbose
```
For production runs keep the generated `*.fastlocate`, `*.fastlcp`, `*.lfgbwt`, and `*.comptext` files on disk to avoid rebuilding.

## Key Files
- `include/fast_lcp.h`, `include/lf_gbwt.h`, `include/compText.h`: succinct indexes layered over a GBWT.
- `include/querySupport.h`, `include/setMaximalMatchQuery.h`, `include/longMatchQuery.h`: long and set-maximal match implementations used by `ibdCaller`.
- `test/ibdCaller.cpp`: command-line driver for IBD discovery.
- `test/generateIndices.cpp`, `test/doTest.cpp`, `test/verifyStructuresAndQueries.cpp`: legacy helpers for regenerating GBWT-derived indexes and validating query correctness.

## Citation
The algorithms powering these tools are described in “Haplotype Matching with GBWT for Pangenome Graphs” by Sanaullah *et al.* (2024; see publication for DOI). Please cite that work when publishing results derived from this repository.
