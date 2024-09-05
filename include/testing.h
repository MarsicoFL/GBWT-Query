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

#ifndef GBWT_QUERY_TESTING_H
#define GBWT_QUERY_TESTING_H

#include<gbwt/gbwt.h>
#include<gbwt/dynamic_gbwt.h>
#include<iostream>
#include<vector>
#include<random>
#include<fast_lcp.h>
#include<lf_gbwt.h>
#include<ioHelp.h>
#include<compText.h>
#include<cmath>
#include<setMaximalMatchQuery.h>
#include<longMatchQuery.h>
#include<cstdio>

//------------------------------------------------------------------------------
//Comparing matches

//Assumes vectors sorted
bool matchesAgree(std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches1,
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches2){
    return matches1 == matches2;
}

//Assumes vectors sorted
bool matchesAgree(std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches1,
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches2){
    if (matches1.size() != matches2.size())
        return false;
    for (unsigned i = 0; i < matches1.size(); ++i)
        if (std::get<0>(matches1[i]) != std::get<0>(matches2[i])
                || std::get<1>(matches1[i]) != std::get<1>(matches2[i])
                || std::get<2>(matches1[i]) != std::get<2>(matches2[i]))
            return false; 
    return true;
}

//Assumes vectors sorted
bool matchesAgree(std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches1,
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches2){
    return matches1 == matches2;
}

//Assumes vectors sorted
bool matchesAgree(std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matches){
    for (unsigned i = 0; i<matches.size()-1; ++i)
        if (!matchesAgree(matches[i], matches[i+1]))
            return false;
    return true;
}

//Assumes vectors sorted
bool matchesAgree(std::vector<std::vector<
       std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>
        >> matches){
    for (unsigned i = 0; i<matches.size()-1; ++i)
        if (!matchesAgree(matches[i], matches[i+1]))
            return false;
    return true;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Testing Queries through incremental GBWT building and querying

//generate haplotype randomly
//starts at the endmarker, generates next node through three options (if in node in GBWT): 
//  take an LF from the current node proportonal probability to numbers of edges (with prob 1/(1+mu+theta))
//  select uniform random node (with prob mu/(1+mu+theta))
//  add new node (with prob theta/(1+mu+theta))
//may generate empty haplotypes
gbwt::vector_type generateHaplotype(const gbwt::GBWT & x, double mu, double theta, std::mt19937_64 & gen, gbwt::node_type & alphabetSize){
    gbwt::size_type currentNodes = x.effective();
    gbwt::node_type currentNode = gbwt::ENDMARKER;
    gbwt::vector_type path;
    std::uniform_real_distribution<> dis1(0, 1+mu+theta), dis2(1, 1+mu+theta);
    std::uniform_int_distribution<> disNode(0,currentNodes-1);

    double p;
    do {
        p = (x.contains(currentNode) && !x.empty(currentNode))? dis1(gen) : dis2(gen);
        if (p < 1){
            std::uniform_int_distribution<> dis(0, x.nodeSize(currentNode)-1);
            currentNode = x.LF(currentNode, dis(gen)).first;
        }
        else if (p < 1+mu){
            currentNode = x.toNode(disNode(gen));
        }
        else{
            currentNode = currentNodes + x.header.offset;
            ++currentNodes;
            disNode = std::uniform_int_distribution<>(0, currentNodes-1);
        }
        path.push_back(currentNode);
    } while (currentNode != gbwt::ENDMARKER);
    path.pop_back();
    alphabetSize = currentNodes + x.header.offset;
    return path;
}

bool queriesEqual(const gbwt::GBWT & x, const gbwt::FastLocate & r, const FastLCP & l, const lf_gbwt::GBWT & lfg, const CompText & ct, const gbwt::vector_type Q){
    //(start position in Q, length, path ID)
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>>                 matchWithoutPosition = {setMaximalMatchQueryN2(x,Q), setMaximalMatchQuery0(x,Q)};
    //(start position in Q, length, path ID, start position in path)
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matchWithPosition = {setMaximalMatchQueryBrute(x,Q), setMaximalMatchQuery1(x, r, Q)};

    matchWithPosition.push_back(setMaximalMatchQuery2(x, r, l, Q));
    matchWithPosition.push_back(setMaximalMatchQuery3(lfg, r, l, Q));
    matchWithPosition.push_back(setMaximalMatchQuery4(lfg, r, l, ct, Q));
    matchWithPosition.push_back(setMaximalMatchQuery2_4(x, r, l, ct, Q));
    for (auto& a : matchWithoutPosition)
        sort(a.begin(), a.end());
    for (auto& a : matchWithPosition)
        sort(a.begin(), a.end());
    return matchesAgree(matchWithoutPosition) 
        && matchesAgree(matchWithPosition) 
        && matchesAgree(matchWithoutPosition[0], matchWithPosition[0]);
}

bool longMatchQueriesEqual(const gbwt::GBWT & x, const gbwt::FastLocate & r, const FastLCP & l, const lf_gbwt::GBWT & lfg, const CompText & ct, const gbwt::vector_type Q, std::mt19937_64& gen, const gbwt::size_type maxTests){
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> LongMatches;
    gbwt::size_type i = 1;
    bool agree, overall = true;
    do {
        LongMatches  = {longMatchQueryBrute(x, i, Q), longMatchQuery2(x, r, l, i, Q), longMatchQuery3(lfg, r, l, i, Q), longMatchQuery4(lfg, r, l, ct, i, Q), longMatchQuery2_4(x, r, l, ct, i, Q)};
        for (auto& a : LongMatches)
            sort(a.begin(), a.end());
        agree = matchesAgree(LongMatches);
        overall = overall && agree;
        i *= 2;
    } while (LongMatches[0].size());
    std::uniform_int_distribution<> dis(1, i), testdis(1,maxTests+1);
    gbwt::size_type num = testdis(gen);
    for (gbwt::size_type i = 0; i < num; ++i) {
        gbwt::size_type L = dis(gen);
        std::cout << "L = " << L << std::endl;
        LongMatches  = {longMatchQueryBrute(x, L, Q), longMatchQuery2(x, r, l, L, Q), longMatchQuery3(lfg, r, l, L, Q), longMatchQuery4(lfg, r, l, ct, L, Q), longMatchQuery2_4(x, r, l, ct, L, Q)};
        for (auto& a : LongMatches)
            sort(a.begin(), a.end());
        agree = matchesAgree(LongMatches);
        overall = overall && agree;
    }
    return overall;
}

template<class T, class U>
bool verifySerializeLoad(const std::string& tempFilename, const T& originalStructure, const U& compare) {
    std::ofstream out(tempFilename);
    if (!out.is_open()) {
        std::cerr << "File '" << tempFilename << "' failed to open for writing!\n";
        assert(false);
        return false;
    }
    originalStructure.serialize(out);
    out.close();
    std::ifstream in(tempFilename);
    if (!in.is_open()) { 
        std::cerr << "File '" << tempFilename << "' failed to open for reading!\n";
        assert(false);
        return false;
    }
    T newStructure;
    newStructure.load(in);
    in.close();
    return newStructure.verify(compare);
}
    

//Testing Queries through incremental GBWT building and querying
bool testIncremental(gbwt::GBWT & x, gbwt::FastLocate & r, FastLCP & l, lf_gbwt::GBWT & lfg, CompText & ct, unsigned n){
    bool overall = true, result, longResult, indexes = true;
    std::string tempFilename = std::tmpnam(nullptr);
    std::cout << "tempFilename: " << tempFilename << std::endl;
    if (!l.verifySuff()) { indexes = false; std::cout << "FastLCP not good!" << std::endl; }
    if (!lfg.verify(x)) {indexes = false; std::cout << "LF GBWT not good!" << std::endl; }
    if (!verifySerializeLoad(tempFilename, lfg, x)) {indexes = false; std::cout << "LF GBWT serialize/load not good!" << std::endl; }
    if (!ct.verifyText()) { indexes = false; std::cout << "CompText not good!" << std::endl; }

    std::random_device rd;
    unsigned seed = rd();
    //unsigned seed = 3126255437;
    //unsigned seed = 0x54459889;
    std::cout << "seed for testIncremental: " << seed << "\n";
    std::mt19937_64 gen(seed);
    gbwt::vector_type Q, hap;
    gbwt::node_type alphabetSize = x.sigma();

    gbwt::GBWTBuilder builder(64, 1);
    auto paths = extractPaths(x);
    for (unsigned i = 0; i < paths.size(); i += 2)
        builder.insert(paths[i], true);

    //auto logbase10 = [] (gbwt::size_type a) { double log = log10(a); return int(log) + (log > int(log)); };
    for (unsigned i = 0; i < n; ++i){
        Q = generateHaplotype(x, 0.05, 0.025, gen, alphabetSize);
        if (Q.size()){
            builder.insert(Q, true);
            builder.finish();

            x = builder.index;
            r = gbwt::FastLocate(x);
            l = FastLCP(r);
            if (!l.verifySuff()) { indexes = false; std::cout << "FastLCP not good!" << std::endl; }
            lfg = lf_gbwt::GBWT(x);
            if (!lfg.verify(x)) {indexes = false; std::cout << "LF GBWT not good!" << std::endl; }
            if (!verifySerializeLoad(tempFilename, lfg, x)) {indexes = false; std::cout << "LF GBWT serialize/load not good!" << std::endl; }
            ct.buildFullMem(l);
            if (!ct.verifyText()) { indexes = false; std::cout << "CompText not good!" << std::endl; }
        }
        //printGBWTandRindex(std::cout, x, r, l, lfg, ct, std::max(3, logbase10(std::max(x.sequences(), x.sigma())) + 1));
        result = queriesEqual(x, r, l, lfg, ct, Q = generateHaplotype(x, 0.05, 0.0, gen, alphabetSize));
        if (!result)
            std::cout << "Set Maximal Match Queries not equal! for Q = " << Q << "!" << std::endl;
        longResult = longMatchQueriesEqual(x, r, l, lfg, ct, Q = generateHaplotype(x, 0.05, 0.0, gen, alphabetSize), gen, 9);
        if (!longResult)
            std::cout << "Long Match Queries not equal! for Q = " << Q << "!" << std::endl;
        overall = overall && result && longResult && indexes;
    }
    return overall;
}

#endif //GBWT_QUERY_TESTING_H
