/*
MIT License

Copyright (c) 2024 Ahsan Sanaullah
Copyright (c) 2024 S. Zhang Lab at UCF

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<string>
#include<algorithm>
#include<iomanip>
#include<vector>
#include<array>


#include<gbwt/gbwt.h>
#include<gbwt/fast_locate.h>
#include<fast_lcp.h>
#include<lf_gbwt.h>
#include<compText.h>

#include<ioHelp.h>
#include<testing.h>
#include<setMaximalMatchQuery.h>

template <class T>
T safeOpenFile(const std::string & name) {
    T file(name);
    if (!file.is_open()) {
        std::cout << "File '" << name << "' failed to open!\n";
        exit(1);
    }
    return file;
}

std::tuple<std::array<double, 7>,bool,double,gbwt::size_type> 
TestSetMaximalMatchQuery(const gbwt::GBWT& gwithDA, const gbwt::GBWT& g, const gbwt::FastLocate& r, const FastLCP& l, const lf_gbwt::GBWT& lfg, const CompText& ct, const gbwt::vector_type& Q) {
    std::array<double, 7> ans;
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matchWithoutPosition;
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matches(5);
    ans[0] = gbwt::readTimer();
    matchWithoutPosition = setMaximalMatchQuery0(gwithDA, Q);
    ans[1] = gbwt::readTimer();
    matches[0] = setMaximalMatchQuery1(g, r, Q);
    ans[2] = gbwt::readTimer();
    matches[1] = setMaximalMatchQuery2(g, r, l, Q);
    ans[3] = gbwt::readTimer();
    matches[2] = setMaximalMatchQuery3(lfg, r, l, Q);
    ans[4] = gbwt::readTimer();
    matches[3] = setMaximalMatchQuery4(lfg, r, l, ct, Q);
    ans[5] = gbwt::readTimer();
    matches[4] = setMaximalMatchQuery2_4(g, r, l, ct, Q);
    ans[6] = gbwt::readTimer();


    //verifying equal
    double start = gbwt::readTimer();
    std::sort(matchWithoutPosition.begin(), matchWithoutPosition.end());
    for (auto& match : matches) { std::sort(match.begin(), match.end()); }
    bool agree = matchesAgree(matches) && matchesAgree(matchWithoutPosition, matches[0]);
    return {ans, agree, gbwt::readTimer()-start, matches[0].size()};
}

std::tuple<std::array<double,5>,bool,double,gbwt::size_type> 
TestLongMatchQuery(const gbwt::GBWT& g, const gbwt::FastLocate& r, const FastLCP& l, const lf_gbwt::GBWT& lfg, const CompText& ct, const gbwt::size_type L, const gbwt::vector_type& Q) {
    std::array<double, 5> ans;
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matches(4);
    ans[0] = gbwt::readTimer();
    matches[0] = longMatchQuery2(g, r, l, L, Q);
    ans[1] = gbwt::readTimer();
    matches[1] = longMatchQuery3(lfg, r, l, L, Q);
    ans[2] = gbwt::readTimer();
    matches[2] = longMatchQuery4(lfg, r, l, ct, L, Q);
    ans[3] = gbwt::readTimer();
    matches[3] = longMatchQuery2_4(g, r, l, ct, L, Q);
    ans[4] = gbwt::readTimer();

    //verifying equal
    double start = gbwt::readTimer();
    for (auto& match : matches) { std::sort(match.begin(), match.end()); }
    bool agree = matchesAgree(matches);
    return {ans, agree, gbwt::readTimer()-start, matches[0].size()};
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "There must be exactly two arguments passed, first the base name for all the data structures, then L. " << argc - 1 << " arguments passed." << std::endl;
        for (int i = 0; i < argc; ++i)
            std::cout << "Argument " << i << ": '" << argv[i] << "'" << std::endl;
        exit(1);
    }
    bool validateStructures = false;
    gbwt::size_type L = atoi(argv[2]);
    std::cout << "L: " << L << std::endl;

    double start = gbwt::readTimer();

    //get query haplotype paths
    std::vector<gbwt::size_type> hapsToRemove;
    sdsl::int_vector<> queryPaths;
    std::ifstream queryIn = safeOpenFile<std::ifstream>(argv[1] + std::string("QueryPaths"));
    sdsl::load(hapsToRemove, queryIn);
    sdsl::load(queryPaths, queryIn);
    queryIn.close();
    std::vector<gbwt::size_type> pathStarts;
    pathStarts.push_back(0);
    for (gbwt::size_type i = 0; i < queryPaths.size(); ++i) { if (queryPaths[i] == gbwt::ENDMARKER) { pathStarts.push_back(i+1); } }
    std::cout << "Finished loading " << pathStarts.size()-1 << " sequences of total length " << queryPaths.size() << " to query in " << gbwt::readTimer() - start << " seconds" << std::endl;

    //open and load gbwts
    start = gbwt::readTimer();
    gbwt::GBWT g,gwithDA;
    std::ifstream gin = safeOpenFile<std::ifstream>(argv[1] + gbwt::GBWT::EXTENSION);
    g.load(gin);
    gin.close();
    std::ifstream gwithDAin = safeOpenFile<std::ifstream>(argv[1] + std::string("_withDA") + gbwt::GBWT::EXTENSION);
    gwithDA.load(gwithDAin);
    gwithDAin.close();
    std::cout << "Read GBWT in " << gbwt::readTimer() - start << " seconds" << std::endl;
    gbwt::printStatistics(g, argv[1] + gbwt::GBWT::EXTENSION, std::cout);
    gbwt::printStatistics(gwithDA, argv[1] + std::string("_withDA") + gbwt::GBWT::EXTENSION, std::cout);
    if (pathStarts.size()-1 != hapsToRemove.size()*(1 + int(g.bidirectional()))) {
        std::cout << "Error, pathStarts.size()-1 = " << pathStarts.size()-1 << ", hapsToRemove.size() = " << hapsToRemove.size() << " and g.bidirectional() = " << g.bidirectional() << std::endl;
        std::cout << "pathStarts.size()-1 should = hapsToRemove.size() if g isn't bidirectional and hapsToRemove.size() * 2 if g is bidirectional" << std::endl;
        exit(1);
    }

    //open and load FastLocate
    start = gbwt::readTimer();
    gbwt::FastLocate fl;
    std::ifstream flin = safeOpenFile<std::ifstream>(argv[1] + gbwt::FastLocate::EXTENSION);
    fl.load(flin);
    fl.index = &g;
    flin.close();
    std::cout << "Read FastLocate in " << gbwt::readTimer() - start << " seconds" << std::endl;
    gbwt::printStatistics(fl, argv[1] + gbwt::FastLocate::EXTENSION);

    //open and load FastLCP
    start = gbwt::readTimer();
    FastLCP flcp;
    std::ifstream flcpin = safeOpenFile<std::ifstream>(argv[1] + FastLCP::EXTENSION);
    flcp.load(flcpin);
    flcp.rindex = &fl;
    flcpin.close();
    std::cout << "Read FastLCP in " << gbwt::readTimer() - start << " seconds" << std::endl;
    //verify FastLCP
    if (validateStructures) {
        start = gbwt::readTimer();
        std::cout << "Verifying FastLCP locatePrev: " << std::flush << ((flcp.verifySuff())? "good" : "not good!" ) << std::endl;
        std::cout << "Verified in " << gbwt::readTimer() - start << " seconds" << std::endl;
    }

    //open and load lf gbwt
    start = gbwt::readTimer();
    lf_gbwt::GBWT lfg;
    std::ifstream lfgin = safeOpenFile<std::ifstream>(argv[1] + lf_gbwt::GBWT::EXTENSION);
    lfg.load(lfgin);
    lfgin.close();
    std::cout << "Read lf gbwt in " << gbwt::readTimer() - start << " seconds" << std::endl;
    //verify lf gbwt
    if (validateStructures) {
        start = gbwt::readTimer();
        std::cout << "Verifying lf gbwt: " << std::flush << ((lfg.verify(g))? "good" : "not good!" ) << std::endl;
        std::cout << "Verified in " << gbwt::readTimer() - start << " seconds" << std::endl;
    }

    //open and load comp text
    start = gbwt::readTimer();
    CompText ct;
    std::ifstream ctin = safeOpenFile<std::ifstream>(argv[1] + CompText::EXTENSION);
    ct.load(ctin);
    ct.source = &flcp;
    ctin.close();
    std::cout << "Read ctin in " << gbwt::readTimer() - start << " seconds" << std::endl;
    //verify comp text 
    if (validateStructures) {
        start = gbwt::readTimer();
        std::cout << "Verifying comp text: " << std::flush << ((ct.verifyText())? "good" : "not good!" ) << std::endl;
        std::cout << "Verified in " << gbwt::readTimer() - start << " seconds" << std::endl;
    }

    //test Queries
    std::cout << "Testing Queries" << std::endl;

    #pragma omp parallel for schedule(dynamic, 1)
    for (gbwt::size_type i = 0; i < pathStarts.size() - 1; ++i){
        double extractStart = gbwt::readTimer();
        //extract path
        gbwt::vector_type Q;
        Q.reserve(pathStarts[i+1]-1-pathStarts[i]);
        for (auto it = queryPaths.begin() + pathStarts[i]; it != queryPaths.begin() + pathStarts[i+1]-1; ++it)
            if (g.contains(*it) && g.nodeSize(*it))
                Q.push_back(*it);


        #pragma omp critical
        {
        std::cout << "Extracted Query " << i << " of length " << Q.size() << " in " << gbwt::readTimer() - extractStart << " seconds" << std::endl;
        }
        //set maximal match query
        auto output = TestSetMaximalMatchQuery(gwithDA, g, fl, flcp, lfg, ct, Q);
        auto& results = std::get<0>(output);
        #pragma omp critical
        {
            std::cout << "Set Maximal Match matches for path " << std::setw(4) << i << ((std::get<1>(output))? " agree. " : " do not agree! ") << "Checked in " << std::setw(12) << std::get<2>(output) << " seconds. There were a total of " << std::setw(10) << std::get<3>(output) << " matches." << std::endl;
            std::cout << "Set Maximal Match query times                 for path " << std::setw(4) << i << " :";
            std::cout << std::setw(10) << results[1] - results[0];
            for (unsigned j = 2; j < results.size(); ++j)
                std::cout << std::setw(10) << results[j] - results[j-1];
            std::cout << std::endl;
        }
        //long match query
        auto lutput = TestLongMatchQuery(g, fl, flcp, lfg, ct, L, Q);
        auto& longResults = std::get<0>(lutput);
        #pragma omp critical
        {
            std::cout << "Long Match  matches for       path " << std::setw(4) << i << ((std::get<1>(lutput))? " agree. " : " do not agree! ") << "Checked in " << std::setw(12) << std::get<2>(lutput) << " seconds. There were a total of " << std::setw(10) << std::get<3>(lutput) << " matches." << std::endl;
            std::cout << "Long Match        query times for L = " << std::setw(7) << L << " for path " << std::setw(4) << i << " :";
            std::cout << std::setw(30) << longResults[1] - longResults[0];
            for (unsigned j = 2; j < longResults.size(); ++j)
                std::cout << std::setw(10) << longResults[j] - longResults[j-1];
            std::cout << std::endl;
        }
    }

    //std::cout << "Testing Queries Done!" << std::endl;
    return 0;
}
