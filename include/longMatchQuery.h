#ifndef GBWT_QUERY_LONGMATCHQUERY_H
#define GBWT_QUERY_LONGMATCHQUERY_H

#include<querySupport.h>
#include<set>
#include<vector>
#include<map>
#include<algorithm>
#include<exception>
#include<ioHelp.h>

//------------------------------------------------------------------------------
//Brute Force Long Match Search 

//adds all long matches starting at i to Lmatches and matches given all path substring 
//P[j,k) s.t. k-j >= L and P[j-1,k) = Q[i-1, i+k-j) is 
//void LMatchesStartingAt(const gbwt::size_type i, const std::vector<gbwt::vector_type>& paths, const gbwt::size_type L,
//        const gbwt::vector_type& Q, std::set<std::tuple<gbwt::size_type,long long,gbwt::size_type>>& Lmatches,
//        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
//    using std::get;
//    auto it = Lmatches.begin();
//    for (gbwt::size_type j = 0; j<paths.size(); ++j){
//        for (long long k = -i; k+i+L <= paths[j].size(); ++k) {
//            //std::cout << "i " << i << " j " << j << " k " << k << std::endl;
//            while (it != Lmatches.end() && get<0>(*it) <= j && get<2>(*it) < i + k)
//                ++it;
//            //std::cout << "*it: (" << std::get<0>(*it) << ',' << std::get<1>(*it) << ',' << std::get<2>(*it) << ")," << std::endl;
//            if (it != Lmatches.end() &&  get<0>(*it) == j && i + get<1>(*it) <= i + k){
//                k = get<2>(*it) - i;
//                continue;
//            }
//            gbwt::size_type matchLen = 0;
//            while (i + matchLen < Q.size() &&  k+i+matchLen < paths[j].size() &&
//                    Q[i + matchLen] == paths[j][k + i + matchLen]) { ++matchLen; }
//            if (matchLen >= L) {
//                Lmatches.emplace(j, k, k + i + matchLen);
//                matches.emplace_back(i, matchLen, j, k + i);
//            }
//        }
//    }
//}
//
//void removeLMatchesNotContaining(const gbwt::size_type i, std::set<std::tuple<gbwt::size_type,long long,gbwt::size_type>>& Lmatches) {
//    using std::get;
//    for (auto it = Lmatches.begin(); it != Lmatches.end();) {
//        auto rem = it++;
//        if (get<2>(*rem) - get<1>(*rem) <= i)
//            Lmatches.erase(rem);
//    }
//}
//
////(start position in Q, length, path ID, start position in path)
//std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
//longMatchQueryBrute(const std::vector<gbwt::vector_type>& paths, const gbwt::size_type L, const gbwt::vector_type& Q){
//    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
//    if (!Q.size()) { return matches; }
//    //Lmatches represents all path substrings P[j,k) that match Q[i, i+k-j) and k-j >= l and Q[i+k-j] != P[k]
//    //a tuple (j,k,l) in Lmatches represents a match between Q and path_j where 
//    //Q[i,l-k) = path_j[i+k,l)
//    std::set<std::tuple<gbwt::size_type,long long,gbwt::size_type>> Lmatches;
//    LMatchesStartingAt(0, paths, L, Q, Lmatches, matches);
//    //std::cout << "Lmatches after compute 0: ";
//    //for (auto j : Lmatches)
//    //    std::cout << '(' << std::get<0>(j) << ',' << std::get<1>(j) << ',' << std::get<2>(j) << "),";
//    //std::cout << std::endl;
//    //std::cout << "matches after compute 0: ";
//    //for (auto j : matches)
//    //    std::cout << '(' << std::get<0>(j) << ',' << std::get<1>(j) << ',' << std::get<2>(j) << ',' << std::get<3>(j) << "),";
//    //std::cout << std::endl;
//    for (gbwt::size_type i = 1; i < Q.size()-1; ++i){
//        removeLMatchesNotContaining(i, Lmatches);
//        //std::cout << "Lmatches after remove not containing " << i << ": ";
//        //for (auto j : Lmatches)
//        //    std::cout << '(' << std::get<0>(j) << ',' << std::get<1>(j) << ',' << std::get<2>(j) << "),";
//        //std::cout << std::endl;
//        LMatchesStartingAt(i, paths, L, Q, Lmatches, matches);
//        //std::cout << "Lmatches after compute " << i << ": ";
//        //for (auto j : Lmatches)
//        //    std::cout << '(' << std::get<0>(j) << ',' << std::get<1>(j) << ',' << std::get<2>(j) << "),";
//        //std::cout << std::endl;
//        //std::cout << "matches after compute 0: ";
//        //for (auto j : matches)
//        //    std::cout << '(' << std::get<0>(j) << ',' << std::get<1>(j) << ',' << std::get<2>(j) << ',' << std::get<3>(j) << "),";
//        //std::cout << std::endl;
//    }
//    return matches;
//}

void LMatchesStartingAt(const gbwt::size_type i, const std::vector<gbwt::vector_type>& paths, const gbwt::size_type L, const gbwt::vector_type& Q, 
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    for (gbwt::size_type j = 0; j < paths.size(); ++j){
        for (gbwt::size_type k = 0; k+L <= paths[j].size(); ++k){
            if (i > 0 && k > 0 && paths[j][k-1] == Q[i-1])
                continue;
            gbwt::size_type matchLen = 0;
            while (k+matchLen < paths[j].size() && i+matchLen < Q.size() && paths[j][k+matchLen] == Q[i+matchLen])
                ++matchLen; 
            if (matchLen >= L)
                matches.emplace_back(i, matchLen, j, k);
        }
    }
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQueryBrute(const std::vector<gbwt::vector_type>& paths, const gbwt::size_type L, const gbwt::vector_type& Q) {
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!Q.size()) { return matches; }
    for (gbwt::size_type i = 0; i+L <= Q.size(); ++i)
        LMatchesStartingAt(i, paths, L, Q, matches);
    return matches;
}

//(start position in Q, length, path ID, start position in path)
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQueryBrute(const gbwt::GBWT& x, const gbwt::size_type L, const gbwt::vector_type& Q) {
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    return longMatchQueryBrute(extractPaths(x), L,Q);
}
//------------------------------------------------------------------------------

//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER Is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery2(const gbwt::GBWT& x, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::size_type L, const gbwt::vector_type& Q) {
    //std::cout << "L2" << std::endl;
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffGBWT(x, r, l, Qs);

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);

    //block is [f^L_i, g^L_i - 1]
    //assumes a.back() = 0 and Qs.back() = gbwt::ENDMARKER
    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i){
        //std::cout << "Iteration i " << i << " of main loop" << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches
        std::tie(block, topSuff, botSuff) = AddLongMatchesFastLCP(x, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        //std::cout << "AddLongMatchesFastLCP values computed" << " block " << "[" << block.first << ',' << block.second << "] topSuff " << r.unpack(topSuff) << " botSuff " << r.unpack(botSuff) << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));

        //std::cout << "Computing next block" << std::endl;
        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if (gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            //std::cout << "Entering empty block manually, before " << block << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], x.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
            //block = {a[i-1], a[i-1] - 1};
            //botSuff = topSuff = s[i-1];
            //if (lcpa[i-1] >= L) {
            //    --block.first;
            //    if (a[i-1] == x.nodeSize(Qs[i-1])) { topSuff = botSuff = l.locateLast(Qs[i-1]); }
            //    else { topSuff = l.locatePrev(topSuff); }
            //    std::cout << "Adding above" << r.unpack(topSuff) << std::endl;
            //    inBlock.emplace(topSuff + Qs.size() - (i-1), i-1 + L);
            //}
            //if (lcpb[i-1] >= L){
            //    ++block.second;
            //    std::cout << "Adding below" << r.unpack(botSuff) << std::endl;
            //    inBlock.emplace(botSuff + Qs.size() - (i-1), i-1 + L);
            //}
            //std::cout << "after " << block << std::endl;
        }
        assert(inBlock.size() == gbwt::Range::length(block));
        if (!gbwt::Range::empty(block)) {
            //std::cout << "Expanding non empty block manually" << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, x.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
            //while (l.LCP(topSuff) >= L) { 
            //    topSuff = l.locatePrev(topSuff); 
            //    std::cout << "Adding above" << r.unpack(topSuff) << std::endl;
            //    --block.first; 
            //    inBlock.emplace(topSuff + Qs.size() - (i-1), i-1 + L);
            //}
            //while (block.second + 1 < x.nodeSize(Qs[i-1]) && 
            //        l.LCP(tempSuff = r.locateNext(botSuff)) >= L) { 
            //    std::cout << "Adding below" << r.unpack(tempSuff) << std::endl;
            //    ++block.second;
            //    botSuff = tempSuff;
            //    inBlock.emplace(botSuff + Qs.size() - (i-1), i-1 + L);
            //}
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery3(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::size_type L, const gbwt::vector_type& Q) {
    //std::cout << "L3" << std::endl;
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowLFGBWT(lfg, a, Qs);

    //std::cout << "lcpa " << lcpa << std::endl;
    //std::cout << "lcpb " << lcpb << std::endl;
    
    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i) {
        //std::cout << "Iteration i " << i << " of main loop" << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches leaving last block and update block to [f^{L+1}_{i-1}, g^{L+1}_{i-1} - 1]
        std::tie(block, topSuff, botSuff) = AddLongMatchesLFGBWT(lfg, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        //std::cout << "AddLongMatchesFastLCP values computed" << " block " << "[" << block.first << ',' << block.second << "] topSuff " << r.unpack(topSuff) << " botSuff " << r.unpack(botSuff) << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));

        //std::cout << "Computing next block" << std::endl;
        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if(gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            //std::cout << "Entering empty block manually, before " << block << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], lfg.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        }
        assert(inBlock.size() == gbwt::Range::length(block));
        if (!gbwt::Range::empty(block)) {
            //std::cout << "Expanding non empty block manually" << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, lfg.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    return matches;
}


//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER Is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery4(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, const CompText& ct, const gbwt::size_type L, const gbwt::vector_type& Q) {
    //std::cout << "L4" << std::endl;
    for(gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(lfg, l, ct, a, s, Qs);
    
    //std::cout << "lcpa " << lcpa << std::endl;
    //std::cout << "lcpb " << lcpb << std::endl;

    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i) {
        //std::cout << "Iteration i " << i << " of main loop" << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches leaving last block and update block to [f^{L+1}_{i-1}, g^{L+1}_{i-1} - 1]
        std::tie(block, topSuff, botSuff) = AddLongMatchesLFGBWT(lfg, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        //std::cout << "AddLongMatchesFastLCP values computed" << " block " << "[" << block.first << ',' << block.second << "] topSuff " << r.unpack(topSuff) << " botSuff " << r.unpack(botSuff) << std::endl;
        //std::cout << "inBlock ";
        //for (auto it = inBlock.begin(); it != inBlock.end(); ++it) {
        //    gbwt::size_type t = it->first - Qs.size() + i;
        //    std::cout << "\t(" << r.seqId(t) << ',' << r.seqOffset(t) << ')';
        //}
        //std::cout << std::endl;

        assert(inBlock.size() == gbwt::Range::length(block));

        //std::cout << "Computing next block" << std::endl;
        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if(gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            //std::cout << "Entering empty block manually, before " << block << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], lfg.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        }
        assert(inBlock.size() == gbwt::Range::length(block));

        //std::cout << "Computing next block" << std::endl;
        if (!gbwt::Range::empty(block)) {
            //std::cout << "Expanding non empty block manually" << std::endl;
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, lfg.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
            //std::cout << "block " << block << " topSuff " << topSuff << " botSuff " << botSuff << std::endl;
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    return matches;
}


//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER Is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery2_4(const gbwt::GBWT& x, const gbwt::FastLocate& r, const FastLCP& l, const CompText& ct, const gbwt::size_type L, const gbwt::vector_type& Q) {
    //std::cout << "L2_4" << std::endl;
    for(gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffGBWT(x, r, l, Qs);

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(x, l, ct, a, s, Qs);

    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i){
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches
        std::tie(block, topSuff, botSuff) = AddLongMatchesFastLCP(x, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));

        assert(inBlock.size() == gbwt::Range::length(block));

        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if (gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], x.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        }
        assert(inBlock.size() == gbwt::Range::length(block));

        if (!gbwt::Range::empty(block)) {
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, x.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    return matches;
}
#endif //GBWT_QUERY_LONGMATCHQUERY_H
