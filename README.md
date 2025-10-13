# GBWT-Query

GBWT-Query provides the data structures and algorithms required to interrogate GBWT haplotype indexes, with a focus on identifying shared haplotypes directly on pangenome graphs. Besides the set maximal and long match query implementations from the accompanying paper, we bundle `ibdCaller`, a CLI that scans every haplotype in a GBWT and reports long shared substrings as candidate identity-by-descent (IBD) segments.

The query algorithms implemented here is from “Haplotype Matching with GBWT for Pangenome Graphs” by Sanaullah *et al.* (2024; see publication for DOI). 

## IBD exploration on a GBWT pangenome
- **What it does** `ibdCaller` walks every haplotype stored in the GBWT, computes long shared matches using the query engines in this repository, and emits them as putative IBD segments.
- **Input requirements** A GBWT file plus the auxiliary indexes (`FastLocate`, `FastLCP`, `lf_gbwt::GBWT`, `CompText`). If the auxiliary structures are missing, `ibdCaller` will build them in-memory from the GBWT so you can prototype immediately.
- **Output interpretation** Results are written as TSV (`hap1_id  hap2_id  hap1_start  hap2_start  length_nodes`). Coordinates are counts of GBWT nodes along the forward haplotypes; map them back to graph coordinates with the GBWT metadata utilities or the graph used to create the index.
- **Scope and caveats** Only forward-orientation threads are reported; reverse hits can be recovered by re-running `ibdCaller` on a reverse GBWT or by adapting the tool. Segment lengths are measured in GBWT nodes rather than nucleotides, so projecting to sequence space requires external tooling (GBZ, ODGI, VG, etc.). Versions 3, 4, and 2_4 assume a bidirectional GBWT.

## Quickstart
1. **Install dependencies** (installs will live one directory above the repos):
   ```bash
   git clone https://github.com/vgteam/sdsl-lite
   git clone https://github.com/jltsiren/gbwt
   cd sdsl-lite && ./install.sh .. && cd ../gbwt && ./install.sh ..
   ```
   This places headers in `../include`, static libraries in `../lib`, and the GBWT utilities in `../gbwt/bin`.
2. **Build GBWT-Query tools**:
   ```bash
   cd GBWT-Query/test
   make ibdCaller          # or just `make` to build every helper
   ```
3. **Call IBD segments on the sample dataset** (long matches ≥4 nodes in the bundled bidirectional example):
   ```bash
   ./ibdCaller ../toyData/indexes/bi_fig1/bi_fig1 4 --version 4 --output bi_fig1_ibd.tsv --verbose
   ```
   For small experiments the auxiliary indexes will be built in-memory. For large cohorts keep the generated `.ri`, `.flcp`, `.lfgbwt`, and `.comptext` files on disk to avoid rebuilding work.

# Implemented Indexes
1. `FastLCP`: Augmentation of FastLocate of that adds the locatePrev and LCP computation capabilities.
2. `lf_gbwt::GBWT`: Barebones GBWT that allows computation of LF in efficient worst case time complexity.
3. `CompText`: Compressed Text Data structure. Implementation of the $O\left(r \log\left(\frac{n}{r}\right)\right)$ space text access data structure by Gagie et al. ("Fully Functional Suffix Trees and Optimal Text Searching in BWT-Runs Bounded Space," https://doi.org/10.1145/3375890). Allows $O\left(\log\left(\frac{n}{r}\right) + l\right)$ time queries of $l$ consecutive characters of the text. This implementation is modified to use less space and allows point queries in $O\left(\log\left(\frac{n}{r}\right)\right)$ time.

# Implemented Haplotype Matching Algorithms
1. Set Maximal Match Query: Given a query path-string $Q$, outputs all matches to a path-string $p$ in the GBWT, $(i,j,k)$, s.t. $Q[i,i+k) = p[j,j+k)$ and there does not exist any path-string $p'$ and index $j'$ in the GBWT s.t. $Q[i-1,i+k) = p'[j'-1, j'+k)$ or $Q[i,i+k+1) = p'[j', j'+k+1)$.
   
3. Long Match Query: Given a query path-string $Q$, outputs all matches to a path-string $p$ in the GBWT, $(i,j,k)$, s.t. $k \geq L$, $Q[i,i+k) = p[j,j+k)$, and $Q[i-1,i+k) \neq p[j-1, j+k)$ and $Q[i,i+k+1) \neq p[j,j+k+1)$.

## Haplotype Matching Algorithm Versions
Each haplotype matching algorithm has various version implementations. Each version uses different indexes as input. Set Maximal Match Query has 6 numbered versions, Long Match Query has 4 (specifically, 2, 3, 4, and 2_4).The GBWT should contain no document array samples except for in version 0. There are also brute force versions for long and set maximal match query. Finally set maximal match query has a naive $O\left(\left\lvert Q\right\rvert^2t_{LF}+cs\right)$ time algorithm that requires just a GBWT (document array samples should be provided). Each numbered version uses the following indexes as input. 
* 0: `GBWT`. Document array samples should be contained in the GBWT, otherwise this implementation run very slowly.
* 1: `GBWT` and `FastLocate`.
* 2: `GBWT`, `FastLocate`, and `FastLCP`.
* 3: `lf_gbwt::GBWT`, `FastLocate`, and `FastLCP`.
* 4: `lf_gbwt::GBWT`, `FastLocate`, `FastLCP`, and `CompText`.
* 2_4: `GBWT`, `FastLocate`, `FastLCP`, and `CompText`.

# Repository contents
The `include/` headers expose the compressed indexes and query routines that underpin the IBD search workflow. They can be embedded in other projects when you need direct access to `FastLCP`, `lf_gbwt::GBWT`, or the long/set maximal match engines.

## Supporting headers
- `fast_lcp.h`, `lf_gbwt.h`, `compText.h`: space-efficient indexes layered on top of a GBWT.
- `querySupport.h`, `setMaximalMatchQuery.h`, `longMatchQuery.h`: building blocks for the haplotype-matching algorithms used by `ibdCaller`.
- `testing.h`, `ioHelp.h`: lightweight helpers for validating indexes and formatting results.

## Legacy utilities
Historical benchmarking and verification programs from the reference publication are kept in [`test/`](/test/). They remain useful for regenerating auxiliary indexes or validating changes:
- `generateIndices.cpp` removes a random subset of haplotypes from a GBWT, rebuilds the auxiliary indexes, and stores the removed paths for benchmarking.
- `doTest.cpp` replays the long/set maximal match benchmarks on trimmed panels and reports agreement between algorithmic variants.
- `verifyStructuresAndQueries.cpp` stress-tests the in-memory constructions of `FastLCP`, `lf_gbwt::GBWT`, and `CompText` after random path insertions.
- `make 1000GGenerate` / `make 1000GTest` recreate the chromosome 21 experiments from the published paper using the bundled 1000 Genomes subset (original dataset DOI: [https://doi.org/10.1038/nature15393](https://doi.org/10.1038/nature15393)).

## IBD workflow with `ibdCaller`
- **Prepare indexes** `ibdCaller` works with just the GBWT file; when auxiliary indexes are missing it builds them in-memory. For production runs keep the generated `*.fastlocate`, `*.fastlcp`, `*.lfgbwt`, and `*.comptext` artifacts on disk to avoid rebuilding.
- **Build and run** Inside `test/`, run `make ibdCaller` and invoke the binary with the shared basename and the minimum segment length (in GBWT node count):

```bash
cd test
make ibdCaller
./ibdCaller ../toyData/indexes/bi_fig1/bi_fig1 50 --version 4 --output bi_fig1_ibd.tsv
```
- **Inspect results** Output is tab-delimited with header `hap1_id  hap2_id  hap1_start  hap2_start  length_nodes`. The caller reports one line per forward-orientation segment of length ≥ `min_length`. Hap IDs refer to internal GBWT sequence identifiers; use the GBWT metadata utilities (`gbwt::Metadata`) or `build_gbwt --list-paths` to translate them to sample / contig labels when available.
- **Iterate** Tune `--max-haps` to restrict processing to the first N haplotypes, and switch `--version` between 2, 3, 4, or 24 depending on which index combination you have available. Versions 3/4/2_4 require a bidirectional GBWT and leverage `lf_gbwt` + `CompText` to accelerate long block extensions.
- **Post-process** Merge or filter segments with your favourite tooling (e.g. pandas, awk, Rust) and project them back to sequence coordinates using the graph that produced the GBWT. Reverse-orientation matches are skipped by default; adapt the code if bi-directional segments are required.
# Compilation
Compilation of code including the header files provided in this repository requires the use of the GBWT library (https://github.com/jltsiren/gbwt). The specific version this code was built on is available at https://github.com/jltsiren/gbwt/blob/0bfeb0723bdc71db075aacf99a77704769d56a55. Follow the instructions in the GBWT readme to compile the GBWT library. The GBWT library (and its dependency, [vgteam's fork of sdsl-lite](https://github.com/vgteam/sdsl-lite)) must be linked in order to compile code that uses header files from this repository. Finally note, the requirements are the same as that of the GBWT library: (C++14, OpenMP).

Sample datasets are provided in [/toyData/indexes](/toyData/indexes). The raw paths for these datasets are in [/toyData/raw/](/toyData/raw). 

Running `make` in the [/test/](/test/) directory builds `ibdCaller` alongside the legacy benchmarking utilities. Use `make ibdCaller` to build just the IBD caller. Running `make test` will run basic checks by verifying structures and replaying the bundled toy experiments. Note, `SDSL_DIR` in the Makefile should be set to the SDSL directory before compilation. `gbwtDir` in the Makefile should be set to the directory that contains the `build_gbwt` binary for `generateIndices` to work correctly. 

## Simple Compilation
1. In an empty directory `mkdir GBWT_AND_QUERY; cd GBWT_AND_QUERY`
2. Clone all repos.
   *  `git clone https://github.com/jltsiren/gbwt`
   *  `git clone https://github.com/ucfcbb/GBWT-Query`
   *  `git clone https://github.com/vgteam/sdsl-lite`
3. Installation of prereqs locally
   *  `cd sdsl-lite`
   *  `./install.sh ../`
   *  `cd ../gbwt`
   *  `./install.sh ../`
4.  Optional: reproduce the 1000G benchmarking setup from the paper
   *  `cd ../GBWT-Query/test`
   *  `make 1000GGenerate 1000GTest`
5.  Optional: call IBD segments on the generated indexes
   *  `cd ../GBWT-Query/test`
   *  `make ibdCaller`
   *  `./ibdCaller ../toyData/indexes/1000G/chr21Trimmed 5000 --version 4 --output chr21_ibd.tsv`
