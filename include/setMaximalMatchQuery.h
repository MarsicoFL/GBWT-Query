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

#ifndef GBWT_QUERY_SETMAXIMALMATCHQUERY_H
#define GBWT_QUERY_SETMAXIMALMATCHQUERY_H

#include<querySupport.h>
#include<vector>
#include<algorithm>
#include<exception>

//------------------------------------------------------------------------------
//Brute Force Set Maximal Match Search

//finds all longest matches to Q starting at position i in Q. matches less than len will not be found
void longestMatchStartingAt(const gbwt::size_type i, const std::vector<gbwt::vector_type>& paths, 
        const gbwt::vector_type& Q, std::set<std::pair<gbwt::size_type,long long>>& starts, gbwt::size_type& len) {
    for (gbwt::size_type j = 0; j<paths.size(); ++j){
        for (long long k = -i; k+i+len <= paths[j].size(); ++k){
            gbwt::size_type matchLen = 0;
            if (!starts.count({j, k})){
                while (i+matchLen < Q.size() && k+i+matchLen < paths[j].size() && 
                        Q[i+matchLen] == paths[j][k+i+matchLen]) ++matchLen;
                if (matchLen > len) { starts.clear(); len = matchLen; }
                if (len && matchLen == len)
                    starts.insert({j,k});
            }
        }
    }
}

//Takes as input haplotype panel paths and query Q. Also i, len, and starts where (j,k) in starts implies path_j[i+k,i+k+len) = Q[i,i+len). 
//All such (j,k) are contained in starts. (Unless len = 0, in which case starts is empty)
//increments i and updates starts and len accordingly. Sets len to largest value possible with nonempty starts
void incrementStart(gbwt::size_type& i, const std::vector<gbwt::vector_type>& paths, 
        const gbwt::vector_type & Q, std::set<std::pair<gbwt::size_type,long long>> & starts, gbwt::size_type & len){
    ++i;
    if (len) { --len; }
    //len = ((len == 0)? (--len) : 0);
    if (len == 0) { starts.clear(); }
    longestMatchStartingAt(i, paths, Q, starts, len);
}

//(start position in Q, length, path ID, start position in path)
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQueryBrute(const std::vector<gbwt::vector_type>& paths, const gbwt::vector_type & Q){
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!Q.size()) { return matches; }
    gbwt::size_type len = 0, oldLen = 0;
    std::set<std::pair<gbwt::size_type,long long>> starts;
    longestMatchStartingAt(0, paths, Q, starts, len);
    for (std::pair<gbwt::size_type,long long> a : starts)
        matches.emplace_back(0, len, a.first, a.second);

    for (gbwt::size_type i = 0; i < Q.size()-1;){
        oldLen = len;
        incrementStart(i, paths, Q, starts, len);
        if (len >= oldLen)
            for (std::pair<gbwt::size_type,long long> a : starts)
                matches.emplace_back(i, len, a.first, i+a.second);
    }
    return matches;
}

//(start position in Q, length, path ID, start position in path)
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQueryBrute(const gbwt::GBWT& x, const gbwt::vector_type & Q) {
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    return setMaximalMatchQueryBrute(extractPaths(x), Q);
}

//------------------------------------------------------------------------------

//(start position in Q, length, path ID)
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQueryN2(const gbwt::GBWT & x, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!Q.size()) { return matches; }

    //Q[i,prevEnd) is the previous longest match
    gbwt::size_type prevEnd = 0, end = 0;
    gbwt::SearchState s, oldS;
    for (gbwt::size_type i = 0; i<Q.size(); ++i){
        //find longest match starting at i, if end point is larger than previous longest match, append set maximal matches to matches
        end = i;
        s = x.find(Q[i]);
        while (!s.empty() && end < Q.size()){
            ++end;
            oldS = s;
            s = x.extend(s, Q[end]);
        }

        if (end > i && end > prevEnd){
            std::vector<gbwt::size_type> hapIds = locate(x, oldS);
            for (unsigned j = 0; j < hapIds.size(); ++j)
                matches.emplace_back(i, end-i, hapIds[j]);
        }
        prevEnd = end;
    }
    
    return matches;
}

//(start position in Q, length, path ID)
//assumes gbwt::ENDMARKER is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQuery0(const gbwt::GBWT & x, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    //compute virtual insertion positions
    std::vector<gbwt::size_type> a = virtualInsertionGBWT(x, Qs);

    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);

    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev)
            AddMatchesGBWT(x, Q, Q.size()-k-len, len, matches);
        prev = len;
    }
    return matches;
}

//(start position in Q, length, path ID, start position in path)
//assumes gbwt::ENDMARKER is not in Q
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQuery1(const gbwt::GBWT & x, const gbwt::FastLocate & r, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    //compute virtual insertion positions
    std::vector<gbwt::size_type> a = virtualInsertionGBWT(x, Qs);

    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);
    
    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev)
            AddMatchesFastLocate(x, r, Q, Q.size()-k-len, len, matches);
        prev = len;
    }
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery2(const gbwt::GBWT & x, const gbwt::FastLocate &r, const FastLCP & l, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffGBWT(x, r, l, Qs);

    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);

    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCP(x, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
        }
        prev = len;
    }
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery3(const lf_gbwt::GBWT & lfg, const gbwt::FastLocate &r, const FastLCP & l, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);

    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowLFGBWT(lfg, a, Qs);

    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCPLFGBWT(lfg, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
        }
        prev = len;
    }
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery4(const lf_gbwt::GBWT & lfg, const gbwt::FastLocate & r, const FastLCP & l, const CompText & ct, const gbwt::vector_type & Q){
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);

    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);

    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(lfg, l, ct, a, s, Qs);

    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCPLFGBWT(lfg, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
        }
        prev = len;
    }
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery2_4(const gbwt::GBWT & x, const gbwt::FastLocate& r, const FastLCP& l, const CompText & ct, const gbwt::vector_type & Q){
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
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(x, l, ct, a, s, Qs);

    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCP(x, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
        }
        prev = len;
    }
    return matches;
}
#endif //GBWT_QUERY_SETMAXIMALMATCHQUERY_H
