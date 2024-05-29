#include<iostream>
#include<stdlib.h>
#include<ioHelp.h>
#include<setMaximalMatchQuery.h>
#include<longMatchQuery.h>
#include<testing.h>
#include<fast_lcp.h>
#include<gbwt/internal.h>
#include<lf_gbwt.h>
#include<compText.h>

const int ioWidth = 3;

int main(int argc, char* argv[]){
    if (argc != 4){
        std::cout << "There must be exactly three arguments passed, first the input gbwt, then its r-index, and then the number of incremental tests to perform. " 
            << argc-1 << " arguments passed.\n";
        exit(1);
    }
    
    std::ifstream gin(argv[1]), rin(argv[2]);
    if (!gin.is_open()){
        std::cout << "Input file '" << argv[1] << "' failed to open!\n";
        exit(1);
    }
    if (!rin.is_open()){
        std::cout << "Input file '" << argv[2] << "' failed to open!\n";
        exit(1);
    }
    
    int numTests = atoi(argv[3]);
    if (numTests < 0) {
        std::cout << "Number of tests needs to be positive! You provided " << argv[3] << ". Which is interpreted as num tests = " << numTests << ".\n";
        exit(1);
    } 

    gbwt::GBWT x;
    x.load(gin);
    std::cout << "Sanity check\n" << gbwt::indexType(x) << '\n';
    gbwt::printStatistics(x, argv[1], std::cout);
    
    gbwt::FastLocate r;
    r.load(rin);
    std::cout << "Sanity check\n" << gbwt::indexType(r) << '\n';
    gbwt::printStatistics(r, argv[2]);
    r.setGBWT(x);
    gbwt::printStatistics(r, argv[2]);

    FastLCP l(r);

    std::cout << "Finished construction of FastLCP" << std::endl;

    lf_gbwt::GBWT lfg(x);

    std::cout << "Finished construction of LF GBWT" << std::endl;
    
    CompText ct(l);
    //ct.buildFullMem(l);

    std::cout << "Finished construction of Compressed Text" << std::endl;


/*    std::cout << "Printing GBWT:" << '\n';
    for (gbwt::node_type i = 0; i < x.sigma(); ++i){
        std::cout << "Printing Node " << i << '\n';
        std::cout << "Size: " << x.nodeSize(i) << '\n';
        gbwt::CompressedRecord a = x.record(i);
        std::vector<gbwt::edge_type> outgoing = x.edges(i);

        std::cout << "BWT: ";
        for (unsigned j = 0; j < x.nodeSize(i); ++j)
            std::cout << a[j] << ' ';
        std::cout << "\nRun: ";
        for (gbwt::CompressedRecordIterator it(a); !(it.end()); ++it)
            std::cout << it->first << ',' << it->second << ' ';
        std::cout << '\n' << "Out: ";
        for (unsigned j = 0; j < outgoing.size(); ++j)
            std::cout << outgoing[j].first << ',' << outgoing[j].second << ' ';
        std::cout << '\n';

    }
 */   
    printGBWTandRindex(std::cout, x, r, l, lfg, ct, ioWidth);

    //gbwt::vector_type a{2,6,8,10,12,14,14,14,18};
    //gbwt::vector_type a{6,8,15,12,18};
    //gbwt::vector_type a{6,19,13,11,9,20,3};
    gbwt::vector_type a{17,15,15,15,15,2,4,8,14,14,14,14,14,14,14,18};
    //(start position in Q, length, path ID)
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>>                 matchWithoutPosition = {setMaximalMatchQueryN2(x,a), setMaximalMatchQuery0(x,a)};
    //(start position in Q, length, path ID, start position in path)
    std::cout << "Finished matches without positons computation" << std::endl;
    //std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matchWithPosition = {setMaximalMatchQueryBrute(x,a), setMaximalMatchQuery1(x, r, a), setMaximalMatchQuery2(x, r, l, a)};
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matchWithPosition = {setMaximalMatchQueryBrute(x,a), setMaximalMatchQuery1(x, r, a)};
    std::cout << "Finished first two matches with positons computation" << std::endl;
    matchWithPosition.push_back(setMaximalMatchQuery2(x, r, l, a));
    std::cout << "Finished first three matches with positons computation" << std::endl;
    matchWithPosition.push_back(setMaximalMatchQuery3(lfg, r, l, a));
    matchWithPosition.push_back(setMaximalMatchQuery4(lfg, r, l, ct, a));
    matchWithPosition.push_back(setMaximalMatchQuery2_4(x, r, l, ct, a));

    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> LongMatches;
    gbwt::size_type i = 1;
    do {
        LongMatches  = {longMatchQueryBrute(x, i, a), longMatchQuery2(x, r, l, i, a), longMatchQuery3(lfg, r, l, i, a), longMatchQuery4(lfg, r, l, ct, i, a), longMatchQuery2_4(x, r, l, ct, i, a)};
        std::cout << "BR: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[0] << "\n";
        std::cout << "L2: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[1] << "\n";
        std::cout << "L3: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[2] << "\n";
        std::cout << "L4: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[3] << "\n";
        std::cout << "L2_4: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[4] << "\n";
        for (auto& a : LongMatches)
            sort(a.begin(), a.end());
        std::cout << "\nAfter sorting:\n\n";
        std::cout << "BR: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[0] << "\n";
        std::cout << "L2: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[1] << "\n";
        std::cout << "L3: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[2] << "\n";
        std::cout << "L4: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[3] << "\n";
        std::cout << "L2_4: L = " << i << " Matches of " << a << ":\n"
            << LongMatches[4] << "\n";
        std::cout << ((matchesAgree(LongMatches))? "Long Matches agree" : "Long Matches don't agree!") << '\n';
        i *= 2;
    } while (LongMatches[0].size());

    std::cout << "Finished matches with positons computation" << std::endl;
    std::cout << "BR: Matches of " << a << ":\n"
        << matchWithPosition[0] << "\n";
    std::cout << "N2: Matches of " << a << ":\n" 
        << matchWithoutPosition[0] << "\n";
    std::cout << "V0: Matches of " << a << ":\n" 
        << matchWithoutPosition[1] << "\n";
    std::cout << "V1: Matches of " << a << ":\n" 
        << matchWithPosition[1] << "\n";
    std::cout << "V2: Matches of " << a << ":\n" 
        << matchWithPosition[2] << "\n";
    std::cout << "V3: Matches of " << a << ":\n" 
        << matchWithPosition[3] << "\n";
    std::cout << "V4: Matches of " << a << ":\n" 
        << matchWithPosition[4] << "\n";
    std::cout << "V2_4: Matches of " << a << ":\n" 
        << matchWithPosition[5] << "\n";
    for (auto& a : matchWithoutPosition)
        sort(a.begin(), a.end());
    for (auto& a : matchWithPosition)
        sort(a.begin(), a.end());
    std::cout << "\nAfter sorting:\n\n";
    std::cout << "BR: Matches of " << a << ":\n"
        << matchWithPosition[0] << "\n";
    std::cout << "N2: Matches of " << a << ":\n" 
        << matchWithoutPosition[0] << "\n";
    std::cout << "V0: Matches of " << a << ":\n" 
        << matchWithoutPosition[1] << "\n";
    std::cout << "V1: Matches of " << a << ":\n" 
        << matchWithPosition[1] << "\n";
    std::cout << "V2: Matches of " << a << ":\n" 
        << matchWithPosition[2] << "\n";
    std::cout << "V3: Matches of " << a << ":\n" 
        << matchWithPosition[3] << "\n";
    std::cout << "V4: Matches of " << a << ":\n" 
        << matchWithPosition[4] << "\n";
    std::cout << "V2_4: Matches of " << a << ":\n" 
        << matchWithPosition[5] << "\n";

    std::cout << ((matchesAgree(matchWithoutPosition))? "Matches without position agree" : "Matches without position don't agree!") << '\n';
    std::cout << ((matchesAgree(matchWithPosition))? "Matches with position agree" : "Matches with position don't agree!") << '\n';
    std::cout << ((matchesAgree(matchWithoutPosition[0], matchWithPosition[0]))? 
            "First matches with position and first matches without position agree" : 
            "First matches with position and first matches without position don't agree!") << '\n';
    //printGBWTandRindex(std::cout, x, r, ioWidth);
    std::cout << "testIncremental: " << ((testIncremental(x, r, l, lfg, ct, numTests))? "passed" : "failed") << '\n';
    //printGBWTandRindex(std::cout, x, r, ioWidth);
    return 0;
}
