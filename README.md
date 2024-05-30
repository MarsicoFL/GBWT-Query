# GBWT-Query

Provides haplotype matching capabilities on the GBWT. Implements long match and set maximal match query on the GBWT. Related paper: "On Haplotype Matching on the GBWT," by Ahsan Sanaullah, Seba Villalobos, Degui Zhi, and Shaojie Zhang. DOI: 

Contact Author: Ahsan Sanaullah (ahsan.sanaullah@ucf.edu)

# Implemented Indexes
1. `FastLCP`: Augmentation of FastLocate of that adds the locatePrev and LCP computation capabilities.
2. `lf_gbwt::GBWT`: Barebones GBWT that allows computation of LF in efficient worst case time complexity.
3. `CompText`: Compressed Text Data structure. Implementation of the $O(r \log(\frac{n}{r}))$ space text access data structure by Gagie et al. ("Fully Functional Suffix Trees and Optimal Text Searching in BWT-Runs Bounded Space," https://doi.org/10.1145/3375890). Allows $O(\log(\frac{n}{r}) + l)$ time queries of l consecutive characters of the text. This implementation is modified to use less space (by a constant factor) and only allows point queries in $O(log(\frac{n}{r}))$ time.

# Implemented Haplotype Matching Algorithms
1. Set Maximal Match Query: Given a query path-string $Q$, outputs all matches to a path-string $p$ in the GBWT, $(i,j,k)$, s.t. $Q[i,i+k) = p[j,j+k)$ and there does not exist any path-string $p'$ and index $j'$ in the GBWT s.t. $Q[i-1,i+k) = p'[j'-1, j'+k)$ or $Q[i,i+k+1) = p'[j', j'+k+1)$.
   
3. Long Match Query: Given a query path-string $Q$, outputs all matches to a path-string $p$ in the GBWT, $(i,j,k)$, s.t. $k \geq L$, $Q[i,i+k) = p[j,j+k)$, and $Q[i-1,i+k) \neq p[j-1, j+k)$ and $Q[i,i+k+1) \neq p[j,j+k+1)$.

## Haplotype Matching Algorithm Versions
Each haplotype matching algorithm has various version implementations. Each version uses different indexes as input. Set Maximal Match Query has 6 numbered versions, Long Match Query has 4 (specifically, 2, 3, 4, and 2_4).The GBWT should contain no document array samples except for in version 0. There are also brute force versions for long and set maximal match query. Finally set maximal match query has a naive $O(|Q|^2t_{LF}+cs)$ time algorithm that requires just a GBWT (document array samples should be provided). Each numbered version uses the following indexes as input. 
* 0: `GBWT`. Document array samples should be contained in the GBWT, otherwise this implementation run very slowly.
* 1: `GBWT` and `FastLocate`.
* 2: `GBWT`, `FastLocate`, and `FastLCP`.
* 3: `lf_gbwt::GBWT`, `FastLocate`, and `FastLCP`.
* 4: `lf_gbwt::GBWT`, `FastLocate`, `FastLCP`, and `CompText`.
* 2_4: `GBWT`, `FastLocate`, `FastLCP`, and `CompText`.

# Purpose of this repository
This repository has two primary purposes. First, it is provided for readers of the paper who wish to view the implementation. Second, the classes/indexes/queries implemented here may be reused in other code. 
## Archival
The code used to perform the experiments in the paper are provided in [/test/](/test/). 
* [generateIndices.cpp](/test/generateIndices.cpp) takes as input a base name and an integer $N$. It reads the GBWT with the base name, then removes $N$ random haplotypes. Finally, it builds the `GBWT`, the `FastLocate`, and all indexes implemented here for the trimmed dataset. Finally, it writes all removed paths and built indexes to the disk.
* [doTest.cpp](/test/doTest.cpp)
## Reuse of class and query implementation
# Compilation
Compilation of code including the header files provided in this repository requires the use of the GBWT library (https://github.com/jltsiren/gbwt). The specific version this code was built on is available at https://github.com/jltsiren/gbwt/blob/0bfeb0723bdc71db075aacf99a77704769d56a55. Follow the instructions in the GBWT readme to compile the GBWT library. The GBWT library (and its dependency, vgteam's fork of sdsl-lite) must be linked in order to compile. Finally note, the requirements are the same as that of the GBWT library: (C++14, OpenMP).
