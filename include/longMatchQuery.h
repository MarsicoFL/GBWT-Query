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
    std::cout << "a: " ;
    for (auto x : a)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "s: " ;
    for (auto x : s)
        std::cout << x << ' ';
    std::cout << std::endl;

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);
    std::cout << "lcpa: " ;
    for (auto x : lcpa)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "lcpb: " ;
    for (auto x : lcpb)
        std::cout << x << ' ';
    std::cout << std::endl;

    //block is [f^L_i, g^L_i - 1]
    //assumes a.back() = 0 and Qs.back() = gbwt::ENDMARKER
    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i){
        std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff: " << topSuff << " botSuff: " << botSuff << std::endl;
        std::cout << "inBlock: ";
        for (auto a: inBlock) 
            std::cout << "(" <<  a.first << ":" << a.second << ") ";
        std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches
        std::tie(block, topSuff, botSuff) = AddLongMatchesFastLCP(x, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        std::cout << "midblock: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff " << topSuff << " botSuff " << botSuff << std::endl;
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
    std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery3(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::size_type L, const gbwt::vector_type& Q) {
    std::cout << "In longMatchQuery3" << std::endl;
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);
    std::cout << "a: " ;
    for (auto x : a)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "s: " ;
    for (auto x : s)
        std::cout << x << ' ';
    std::cout << std::endl;

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowLFGBWT(lfg, a, Qs);
    std::cout << "lcpa: " ;
    for (auto x : lcpa)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "lcpb: " ;
    for (auto x : lcpb)
        std::cout << x << ' ';
    std::cout << std::endl;

    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i) {
        std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff: " << topSuff << " botSuff: " << botSuff << std::endl;
        std::cout << "inBlock: ";
        for (auto a: inBlock) 
            std::cout << "(" <<  a.first << ":" << a.second << ") ";
        std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches leaving last block and update block to [f^{L+1}_{i-1}, g^{L+1}_{i-1} - 1]
        std::tie(block, topSuff, botSuff) = AddLongMatchesLFGBWT(lfg, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        std::cout << "midblock: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff " << topSuff << " botSuff " << botSuff << std::endl;
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        assert(inBlock.size() == gbwt::Range::length(block));

        //std::cout << "Computing next block" << std::endl;
        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if(gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], lfg.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        }
        assert(inBlock.size() == gbwt::Range::length(block));
        if (!gbwt::Range::empty(block)) {
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, lfg.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    std::cout << "Leaving longMatchQuery3" << std::endl;
    return matches;
}


//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER Is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery4(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, const CompText& ct, const gbwt::size_type L, const gbwt::vector_type& Q) {
    std::cout << "In longMatchQuery4" << std::endl;
    for(gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);
    std::cout << "a: " ;
    for (auto x : a)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "s: " ;
    for (auto x : s)
        std::cout << x << ' ';
    std::cout << std::endl;

    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(lfg, l, ct, a, s, Qs);
    std::cout << "lcpa: " ;
    for (auto x : lcpa)
        std::cout << x << ' ';
    std::cout << std::endl;
    std::cout << "lcpb: " ;
    for (auto x : lcpb)
        std::cout << x << ' ';
    std::cout << std::endl;
    
    gbwt::range_type block = gbwt::Range::empty_range();
    gbwt::size_type topSuff = gbwt::invalid_offset(), botSuff = gbwt::invalid_offset();
    std::map<gbwt::size_type,gbwt::size_type> inBlock;
    for (gbwt::size_type i = Qs.size()-1; i != 0; --i) {
        std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff: " << topSuff << " botSuff: " << botSuff << std::endl;
        std::cout << "inBlock: ";
        for (auto a: inBlock) 
            std::cout << "(" <<  a.first << ":" << a.second << ") ";
        std::cout << std::endl;
        assert(inBlock.size() == gbwt::Range::length(block));
        //output matches leaving last block and update block to [f^{L+1}_{i-1}, g^{L+1}_{i-1} - 1]
        std::tie(block, topSuff, botSuff) = AddLongMatchesLFGBWT(lfg, r, l, inBlock, i, Qs, block, topSuff, botSuff, matches);
        std::cout << "midblock: [" << block.first << "," << block.second << "]" << std::endl;
        std::cout << "topSuff " << topSuff << " botSuff " << botSuff << std::endl;
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
        assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
        //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
        assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        assert(inBlock.size() == gbwt::Range::length(block));

        //compute block [f^L_{i-1}, g^L_{i-1} - 1]
        if(gbwt::Range::empty(block) && std::max(lcpa[i-1], lcpb[i-1]) >= L) {
            std::tie(block, topSuff, botSuff) = LongMatchesFastLCPInitializeEmptyBlock(r, l, inBlock, a[i-1], Qs[i-1], lfg.nodeSize(Qs[i-1]), s[i-1], lcpa[i-1], lcpb[i-1], Qs.size(), i - 1, L);
            //assert if block empty topSuff, botSuff = gbwt::invalid_offset()
            assert(!gbwt::Range::empty(block) || (botSuff == gbwt::invalid_offset() && topSuff == gbwt::invalid_offset()));
            //assert if block not empty, topsuff, botsuff != gbwt::invalid_offset()
            assert( gbwt::Range::empty(block) || (botSuff != gbwt::invalid_offset() && topSuff != gbwt::invalid_offset()));
        }
        assert(inBlock.size() == gbwt::Range::length(block));

        if (!gbwt::Range::empty(block)) {
            std::tie(block, topSuff, botSuff) = LongMatchesExpandNonEmptyBlock(r, l, inBlock, lfg.nodeSize(Qs[i-1]), block, topSuff, botSuff, Qs.size(), i - 1, L);
        }
        assert(inBlock.size() == gbwt::Range::length(block));
    }
    std::cout << "block: [" << block.first << "," << block.second << "]" << std::endl;
    AddLongMatchesWholeBlock(r, l, inBlock, 0, Qs, block, topSuff, matches);
    assert(inBlock.empty());
    std::cout << "Leaving longMatchQuery4" << std::endl;
    return matches;
}


//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER Is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
longMatchQuery2_4(const gbwt::GBWT& x, const gbwt::FastLocate& r, const FastLCP& l, const CompText& ct, const gbwt::size_type L, const gbwt::vector_type& Q) {
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
