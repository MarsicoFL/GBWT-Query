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
    //return matches1 == matches2;
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
    //std::cout << "begGen" << std::flush;

    double p;
    do {
        p = (x.contains(currentNode) && !x.empty(currentNode))? dis1(gen) : dis2(gen);
        if (p < 1){
            std::uniform_int_distribution<> dis(0, x.nodeSize(currentNode)-1);
            currentNode = x.LF(currentNode, dis(gen)).first;
            //std::cout << "p," << currentNode << ',';
        }
        else if (p < 1+mu){
            currentNode = x.toNode(disNode(gen));
            //std::cout << "m," << currentNode << ',';
        }
        else{
            currentNode = currentNodes + x.header.offset;
            ++currentNodes;
            disNode = std::uniform_int_distribution<>(0, currentNodes-1);
            //std::cout << "t," << currentNode << ',';
        }
        path.push_back(currentNode);
    } while (currentNode != gbwt::ENDMARKER);
    //std::cout << "\n";
    path.pop_back();
    //std::cout << "endGen" << std::flush;
    alphabetSize = currentNodes + x.header.offset;
    return path;
}

bool queriesEqual(const gbwt::GBWT & x, const gbwt::FastLocate & r, const FastLCP & l, const lf_gbwt::GBWT & lfg, const CompText & ct, const gbwt::vector_type Q){
    //std::cout << "Q: " << Q << '\n';
    //(start position in Q, length, path ID)
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>>                 matchWithoutPosition = {setMaximalMatchQueryN2(x,Q), setMaximalMatchQuery0(x,Q)};
    //(start position in Q, length, path ID, start position in path)
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> matchWithPosition = {setMaximalMatchQueryBrute(x,Q), setMaximalMatchQuery1(x, r, Q)};

    matchWithPosition.push_back(setMaximalMatchQuery2(x, r, l, Q));
    matchWithPosition.push_back(setMaximalMatchQuery3(lfg, r, l, Q));
    matchWithPosition.push_back(setMaximalMatchQuery4(lfg, r, l, ct, Q));
    matchWithPosition.push_back(setMaximalMatchQuery2_4(x, r, l, ct, Q));
    //auto s = extractPaths(x);
    //gbwt::size_type sample = r.locateFirst(gbwt::ENDMARKER);
/*    for (unsigned i = 0; i < s.size(); ++i){
        std::cout << "Path " <<  i << ", Length: (" << s[i].size() << ",(" 
            << r.seqId(firstSampleFromPath(r, i)) << ',' << r.seqOffset(firstSampleFromPath(r, i)) << "),(" 
            << r.seqId(sample) << ',' << r.seqOffset(sample) << ")): " << s[i] << std::endl;
        if (s[i].size() != r.seqOffset(firstSampleFromPath(r, i)) 
                    || s[i].size() != r.seqOffset(sample))
                throw std::logic_error("Lengths don't match for path");
        sample = r.locateNext(sample);
    }
    std::cout << "ENDMARKER has " << r.comp_to_run[1] << " runs" << std::endl;
    std::cout << "Outputting r-index:" << std::endl;
    for (unsigned i = 0; i < r.samples.size(); ++i){
        std::cout << "(" << r.seqId(r.samples[i]) 
            << ',' << r.seqOffset(r.samples[i]) << ")," << std::flush;
    }
    std::cout << std::endl;
    */
    for (auto& a : matchWithoutPosition)
        sort(a.begin(), a.end());
    for (auto& a : matchWithPosition)
        sort(a.begin(), a.end());
    //for (auto& a : matchWithoutPosition)
        //std::cout << a << '\n';
    //for (auto& a : matchWithPosition)
        //std::cout << a << '\n';
    return matchesAgree(matchWithoutPosition) 
        && matchesAgree(matchWithPosition) 
        && matchesAgree(matchWithoutPosition[0], matchWithPosition[0]);
}

bool longMatchQueriesEqual(const gbwt::GBWT & x, const gbwt::FastLocate & r, const FastLCP & l, const lf_gbwt::GBWT & lfg, const CompText & ct, const gbwt::vector_type Q, std::mt19937_64& gen, const gbwt::size_type maxTests){
    std::vector<std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>> LongMatches;
    std::cout << "Q: " << Q << '\n';
    gbwt::size_type i = 1;
    bool agree, overall = true;
    do {
        std::cout << "L = " << i << std::endl;
        LongMatches  = {longMatchQueryBrute(x, i, Q), longMatchQuery2(x, r, l, i, Q), longMatchQuery3(lfg, r, l, i, Q), longMatchQuery4(lfg, r, l, ct, i, Q), longMatchQuery2_4(x, r, l, ct, i, Q)};
        //std::cout << "BR: L = " << i << " Matches of " << a << ":\n"
        //    << LongMatches[0] << "\n";
        //std::cout << "L2: L = " << i << " Matches of " << a << ":\n"
        //    << LongMatches[1] << "\n";
        for (auto& a : LongMatches)
            sort(a.begin(), a.end());
        //std::cout << "\nAfter sorting:\n\n";
        //std::cout << "BR: L = " << i << " Matches of " << a << ":\n"
        //    << LongMatches[0] << "\n";
        //std::cout << "L2: L = " << i << " Matches of " << a << ":\n"
        //    << LongMatches[1] << "\n";
        agree = matchesAgree(LongMatches);
        overall = overall && agree;
        //std::cout << ((agree)? "Long Matches agree" : "Long Matches don't agree!") << '\n';
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

//Testing Queries through incremental GBWT building and querying
bool testIncremental(gbwt::GBWT & x, gbwt::FastLocate & r, FastLCP & l, lf_gbwt::GBWT & lfg, CompText & ct, unsigned n){
    bool overall = true, result, longResult, indexes = true;
    std::random_device rd;
    //unsigned seed = rd();
    //unsigned seed = 3157635622;
    //unsigned seed = 3127984630;
    //unsigned seed = 2823761535;
    unsigned seed = 1376515986;
    std::cout << "seed for testIncremental: " << seed << "\n";
    std::mt19937_64 gen(seed);
    gbwt::vector_type Q, hap;
    gbwt::node_type alphabetSize = x.sigma();

    gbwt::GBWTBuilder builder(64, 1);
    auto paths = extractPaths(x);
    for (unsigned i = 0; i < paths.size(); i += 2)
        builder.insert(paths[i], true);
    //gbwt::DynamicGBWT& d = builder.index;

    //auto logbase10 = [] (gbwt::size_type a) { double log = log10(a); return int(log) + (log > int(log)); };
    for (unsigned i = 0; i < n; ++i){
        Q = generateHaplotype(x, 0.05, 0.025, gen, alphabetSize);
        if (Q.size()){
            //hap = Q;
            //hap.push_back(gbwt::ENDMARKER);
            //gbwt::reversePath(Q, hap);
            //hap.push_back(gbwt::ENDMARKER);
            //std::cout << "seq: " << d.header.sequences << " siz: " << d.header.size << " alp:" << d.header.alphabet_size << std::endl;
            builder.insert(Q, true);
            //std::cout << "seq: " << d.header.sequences << " siz: " << d.header.size << " alp:" << d.header.alphabet_size << std::endl;
            //d.header.sequences += 2;
            //d.header.size += hap.size();
            //d.header.alphabet_size = alphabetSize;
            builder.finish();
            x = builder.index;
            //std::cout << "seq: " << x.header.sequences << " siz: " << x.header.size << " alp:" << x.header.alphabet_size << std::endl;
            r = gbwt::FastLocate(x);
            l = FastLCP(r);
            if (!l.verifySuff()) { indexes = false; std::cout << "FastLCP not good!" << std::endl; }
            lfg = lf_gbwt::GBWT(x);
            if (!lfg.verify(x)) {indexes = false; std::cout << "LF GBWT not good!" << std::endl; }
            //CompText tempCT(l);
            //std::cout << "default build compressed text structure is " << std::flush << ((tempCT.verifyText())? "good" : "not good!") << std::endl;
            //tempCT.buildFullMem(l);
            //ct = tempCT;
            ct.buildFullMem(l);
            if (!ct.verifyText()) { indexes = false; std::cout << "CompText not good!" << std::endl; }
            //ct = CompText(l);
            //std::cout << " len: " << r.header.max_length << std::endl;
        }
        //printGBWTandRindex(std::cout, x, r, l, lfg, ct, std::max(3, logbase10(std::max(x.sequences(), x.sigma())) + 1));
        //double start = gbwt::readTimer();
        //std::cout << "compressed text structure is " << std::flush << ((ct.verifyText())? "good" : "not good!") << std::endl;
        //std::cout << "verifying compressed text took " << gbwt::readTimer()-start << " seconds" << std::endl;
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
