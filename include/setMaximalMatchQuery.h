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
    //std::cout << "SB" << std::endl;
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    return setMaximalMatchQueryBrute(extractPaths(x), Q);
}

//------------------------------------------------------------------------------

//(start position in Q, length, path ID)
std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> 
setMaximalMatchQueryN2(const gbwt::GBWT & x, const gbwt::vector_type & Q){
    //std::cout << "SN" << std::endl;
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
        //std::cout << i << ',' << end << ',' << s  << ',' << oldS << std::endl;
        while (!s.empty() && end < Q.size()){
            ++end;
            oldS = s;
            s = x.extend(s, Q[end]);
            //std::cout << i << ',' << end << ',' << s  << ',' << oldS << std::endl;
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
    //std::cout << "S0" << std::endl;
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);
    //gbwt::vector_type Qs = Q;
    ////reverse Q and append the endmarker
    //for (unsigned i = 0; i<Qs.size()/2; ++i)
    //    std::swap(Qs[i], Qs[Qs.size()-1-i]);
    //Qs.push_back(gbwt::ENDMARKER); 
    //std::cout << "Query path after reversal and endmarker appending: " << Qs << "\n";


    //std::cout << "Computing virtual insertion positions" << std::endl;
    //compute virtual insertion positions
    std::vector<gbwt::size_type> a = virtualInsertionGBWT(x, Qs);
    //std::vector<gbwt::size_type> a(Qs.size());
    //int i = Qs.size()-1;
    //a[i] = 0;
    //for (--i; i >= 0; --i)
    //    if ((a[i] = x.LF(Qs[i+1], a[i+1], Qs[i])) == gbwt::invalid_offset())
    //        a[i] = 0;
    //std::cout << "Offsets within nodes for Q*: " << a << '\n';

    //std::cout << "Computing lcps" << std::endl;
    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);
    //gbwt::size_type computed_above, computed_below;
    //std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    //for (gbwt::size_type k = 0; k < Qs.size(); ++k){
    //    if (k == 0 || computed_above < k){
    //        gbwt::edge_type above = {Qs[k], a[k]-1};
    //        gbwt::size_type lcp = 0;
    //        computed_above = k;
    //        while (a[k+lcp] != 0 
    //                && above.first == Qs[k+lcp] 
    //                && above.first != gbwt::ENDMARKER 
    //                && above.second >= 0){
    //            if (above.second+1 == a[k+lcp])
    //                computed_above = k+lcp;
    //            ++lcp;
    //            above = x.inverseLF(above); 
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_above; ++kp)
    //            lcpa[kp] = lcp--;
    //    }
    //    if (k == 0 || computed_below < k){
    //        gbwt::edge_type below = {Qs[k], a[k]};
    //        gbwt::size_type lcp = 0;
    //        computed_below = k;
    //        while (a[k+lcp] != x.nodeSize(Qs[k+lcp]) 
    //                && below.first == Qs[k+lcp] 
    //                && below.first != gbwt::ENDMARKER 
    //                && below.second < x.nodeSize(below.first)){
    //            if (below.second == a[k+lcp])
    //                computed_below = k+lcp;
    //            ++lcp;
    //            below = x.inverseLF(below);
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_below; ++kp)
    //            lcpb[kp] = lcp--;
    //    }
    //}
    //std::cout << "lcpa: " << lcpa << std::endl << "lcpb: " << lcpb << std::endl;

    //std::cout << "Computing set maximal matches" << std::endl;
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
    //std::cout << "S1" << std::endl;
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);
    //gbwt::vector_type Qs = Q;
    ////reverse Q and append the endmarker
    //for (unsigned i = 0; i<Qs.size()/2; ++i)
    //    std::swap(Qs[i], Qs[Qs.size()-1-i]);
    //Qs.push_back(gbwt::ENDMARKER); 
    //std::cout << "Query path after reversal and endmarker appending: " << Qs << "\n";


    //std::cout << "Computing virtual insertion positions" << std::endl;
    //compute virtual insertion positions
    std::vector<gbwt::size_type> a = virtualInsertionGBWT(x, Qs);
    //std::vector<gbwt::size_type> a(Qs.size());
    //int i = Qs.size()-1;
    //a[i] = 0;
    //for (--i; i >= 0; --i)
    //    if ((a[i] = x.LF(Qs[i+1], a[i+1], Qs[i])) == gbwt::invalid_offset())
    //        a[i] = 0;
    //std::cout << "Offsets within nodes for Q*: " << a << '\n';

    //std::cout << "Computing lcps" << std::endl;
    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);
    //gbwt::size_type computed_above, computed_below;
    //std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    //for (gbwt::size_type k = 0; k < Qs.size(); ++k){
    //    if (k == 0 || computed_above < k){
    //        gbwt::edge_type above = {Qs[k], a[k]-1};
    //        gbwt::size_type lcp = 0;
    //        computed_above = k;
    //        while (a[k+lcp] != 0 
    //                && above.first == Qs[k+lcp] 
    //                && above.first != gbwt::ENDMARKER 
    //                && above.second >= 0){
    //            if (above.second+1 == a[k+lcp])
    //                computed_above = k+lcp;
    //            ++lcp;
    //            above = x.inverseLF(above); 
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_above; ++kp)
    //            lcpa[kp] = lcp--;
    //    }
    //    if (k == 0 || computed_below < k){
    //        gbwt::edge_type below = {Qs[k], a[k]};
    //        gbwt::size_type lcp = 0;
    //        computed_below = k;
    //        while (a[k+lcp] != x.nodeSize(Qs[k+lcp]) 
    //                && below.first == Qs[k+lcp] 
    //                && below.first != gbwt::ENDMARKER 
    //                && below.second < x.nodeSize(below.first)){
    //            if (below.second == a[k+lcp])
    //                computed_below = k+lcp;
    //            ++lcp;
    //            below = x.inverseLF(below);
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_below; ++kp)
    //            lcpb[kp] = lcp--;
    //    }
    //}
    //std::cout << "lcpa: " << lcpa << std::endl << "lcpb: " << lcpb << std::endl;

    //std::cout << "Computing set maximal matches" << std::endl;
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

/*std::pair<gbwt::size_type,gbwt::size_type>
safeLF(const FastLCP & l, gbwt::size_type ind, gbwt::size_type suff, 
        gbwt::node_type from, gbwt::node_type to)
{
    //std::cout << "safeLF from " << from << " to " << to << " with index " << ind 
        //<< " and suff " << l.rindex->unpack(suff) << std::endl;
    gbwt::CompressedRecord rec = l.rindex->index->record(from);

    gbwt::size_type outrank = rec.edgeTo(to);
    //std::cout << "outrank " << outrank << std::flush;
    if(outrank >= rec.outdegree()) 
        return {0, l.getSample(to, 0)};
    //std::cout << " outdegree " << rec.outdegree() << std::endl;

    gbwt::CompressedRecordRankIterator iter(rec, outrank);

    gbwt::size_type last_run_before = gbwt::invalid_offset();

    if(iter->first == outrank)
        last_run_before = iter.runId();

    while(!(iter.end()) && iter.offset() < ind)
    {
        ++iter;
        if(iter->first == outrank)
            last_run_before = iter.runId();
    }
    //std::cout << "After while loop, *iter " << *iter << " iter.rank " << iter.rank() 
        //<< " iter.offset " << iter.offset() << " iter.runId " << iter.runId()
        //<< " last_run_before " << last_run_before << std::endl;
    gbwt::size_type temp = iter.rank();
    if (ind < iter.offset() && iter->first == outrank)
        return {temp - (iter.offset() - ind), --suff};

    if (last_run_before != gbwt::invalid_offset())
        return {temp, 
            ((temp == l.rindex->index->nodeSize(to))? 
                gbwt::invalid_offset() : 
                l.rindex->locateNext(l.getSampleBot(from, last_run_before)-1))};

    while (iter->first != outrank)
        ++iter;
    //std::cout << "After 2nd while loop, iter " << *iter << " iter.rank " << iter.rank() 
        //<< " iter.offset " << iter.offset() << " iter.runId " << iter.runId()
        //<< std::endl;

    return {temp, l.getSample(from, iter.runId())-1};
}*/

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery2(const gbwt::GBWT & x, const gbwt::FastLocate &r, const FastLCP & l, const gbwt::vector_type & Q){
    //std::cout << "S2" << std::endl;
    for (gbwt::node_type a : Q)
        if (!x.contains(a) || x.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!x.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);
    //gbwt::vector_type Qs = Q;
    ////reverse Q and append the endmarker
    //for (unsigned i = 0; i<Qs.size()/2; ++i)
    //    std::swap(Qs[i], Qs[Qs.size()-1-i]);
    //Qs.push_back(gbwt::ENDMARKER); 
    //std::cout << "Query path after reversal and endmarker appending: " << Qs << "\n";


    //std::cout << "Computing virtual insertion positions" << std::endl;
    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffGBWT(x, r, l, Qs);
    //std::vector<gbwt::size_type> a(Qs.size());
    //std::vector<gbwt::size_type> s(Qs.size());
    //int i = Qs.size()-1;
    //a[i] = 0;
    //s[i] = l.getSample(gbwt::ENDMARKER, 0);
    //gbwt::CompressedRecord rec;

    //for (--i; i >= 0; --i)
    //    std::tie(a[i], s[i]) = safeLF(l, a[i+1], s[i+1], Qs[i+1], Qs[i]);
    //std::cout << "Offsets within nodes for Q*: " << a << '\n';
    //std::cout << "Suffixes within nodes for Q*: ";
    //if (s.size()){
        //std::cout << r.unpack(*s.begin());
        //std::for_each(++s.begin(), s.end(), [&r](const auto & x){std::cout << ',' << r.unpack(x);});
    //}
    //std::cout << '\n';
    //std::cout << "Suffixes (manual) within nodes for Q*: ";
    //if (s.size()){
        //gbwt::size_type suff = (r.decompressSA(Qs[0]).size() > *a.begin())?
                //r.decompressSA(Qs[0])[*a.begin()]:
                //gbwt::invalid_offset();
        //std::cout << r.unpack(suff);
        //for (unsigned i = 1; i < a.size(); ++i){
            //suff = (r.decompressSA(Qs[i]).size() > a[i])?
                //r.decompressSA(Qs[i])[a[i]]:
                //gbwt::invalid_offset();
            //std::cout << ',' << r.unpack(suff);
        //}
    //}
    //std::cout << '\n';

    //std::cout << "Computing lcps" << std::endl;
    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowGBWT(x, a, Qs);
    //gbwt::size_type computed_above, computed_below;
    //std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    //for (gbwt::size_type k = 0; k < Qs.size(); ++k){
    //    if (k == 0 || computed_above < k){
    //        gbwt::edge_type above = {Qs[k], a[k]-1};
    //        gbwt::size_type lcp = 0;
    //        computed_above = k;
    //        while (a[k+lcp] != 0 
    //                && above.first == Qs[k+lcp] 
    //                && above.first != gbwt::ENDMARKER 
    //                && above.second >= 0){
    //            if (above.second+1 == a[k+lcp])
    //                computed_above = k+lcp;
    //            ++lcp;
    //            above = x.inverseLF(above); 
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_above; ++kp)
    //            lcpa[kp] = lcp--;
    //    }
    //    if (k == 0 || computed_below < k){
    //        gbwt::edge_type below = {Qs[k], a[k]};
    //        gbwt::size_type lcp = 0;
    //        computed_below = k;
    //        while (a[k+lcp] != x.nodeSize(Qs[k+lcp]) 
    //                && below.first == Qs[k+lcp] 
    //                && below.first != gbwt::ENDMARKER 
    //                && below.second < x.nodeSize(below.first)){
    //            if (below.second == a[k+lcp])
    //                computed_below = k+lcp;
    //            ++lcp;
    //            below = x.inverseLF(below);
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_below; ++kp)
    //            lcpb[kp] = lcp--;
    //    }
    //}
    //std::cout << "lcpa: " << lcpa << std::endl << "lcpb: " << lcpb << std::endl;

    //std::cout << "Computing set maximal matches" << std::endl;
    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCP(x, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
            //gbwt::size_type curr = lcpa[k], plen,
            //    curr_s = (a[k] == x.nodeSize(Qs[k]))? 
            //        gbwt::invalid_offset():
            //        s[k];

            //while (curr == len){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //    curr_s = (curr_s == gbwt::invalid_offset())?
            //        l.locateLast(Qs[k]):
            //        l.locatePrev(curr_s);
            //    //std::cout << "prev, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs " ;
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    curr = std::min(curr, l.LCP(curr_s));
            //}
            //curr = lcpb[k];
            //curr_s = s[k];
            //while (curr == len){
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs ";
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    ++a[k];
            //    if (a[k] < x.nodeSize(Qs[k])){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //        curr_s = r.locateNext(curr_s);
            //    //std::cout << "next, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //        curr = std::min(curr, l.LCP(curr_s));
            //    } 
            //    else
            //        break;
            //}
        }
        prev = len;
    }
    return matches;
}

/*std::pair<gbwt::size_type,gbwt::size_type>
safeLF(const lf_gbwt::GBWT & lfg, const FastLCP & l, gbwt::size_type ind, gbwt::size_type suff, 
        gbwt::node_type from, gbwt::node_type to)
{
    //std::cout << "safeLF from " << from << " to " << to << " with index " << ind 
        //<< " and suff " << l.rindex->unpack(suff) << std::endl;
    const lf_gbwt::CompressedRecord & rec = lfg.record(from);

    gbwt::size_type outrank = rec.edgeTo(lfg.toComp(to)),
        newind = rec.LF(ind, lfg.toComp(to));
    //std::cout << "outrank " << outrank << std::flush;
    if(outrank >= rec.outdegree()) 
        return {0, l.getSample(to, 0)};

    if (ind < rec.size() && rec[ind] == lfg.toComp(to))
        return {newind, --suff};

    auto nextRun = rec.firstByAlphabet.successor(outrank*rec.size() + ind);
    if (nextRun->second < (outrank+1)*rec.size()){
        //successor run exists
        //std::cout << "Using successor run" << std::endl;
        gbwt::size_type run_id = rec.logicalRunId(nextRun->second - outrank*rec.size());
        return {newind, l.getSample(from, run_id) - 1};
    }
    else if (newind < lfg.record(to).size()){
        //predecessor run exists
        //std::cout << "Using predecessor run" << std::endl;
        auto prevRun = --nextRun;
        //std::cout << prevRun->first << ',' << prevRun->second << std::endl;
        gbwt::size_type run_id = rec.logicalRunId(prevRun->second - outrank*rec.size());
        //std::cout << "run_id: " << run_id << std::endl;
        //std::cout << "l.getSampleBot(from, run_id) - 1" << l.rindex->unpack(l.getSampleBot(from, run_id) -1) << std::endl;
        return {newind, l.rindex->locateNext(l.getSampleBot(from, run_id)-1)};
    }
    //std::cout << "nothing found?" << std::endl;
    return {newind, gbwt::invalid_offset()};
}*/

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery3(const lf_gbwt::GBWT & lfg, const gbwt::FastLocate &r, const FastLCP & l, const gbwt::vector_type & Q){
    //std::cout << "S3" << std::endl;
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);
    //gbwt::vector_type Qs = Q;
    ////reverse Q and append the endmarker
    //for (unsigned i = 0; i<Qs.size()/2; ++i)
    //    std::swap(Qs[i], Qs[Qs.size()-1-i]);
    //Qs.push_back(gbwt::ENDMARKER); 
    //std::cout << "Query path after reversal and endmarker appending: " << Qs << "\n";


    //std::cout << "Computing virtual insertion positions" << std::endl;
    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);
    //std::vector<gbwt::size_type> a(Qs.size());
    //std::vector<gbwt::size_type> s(Qs.size());
    //int i = Qs.size()-1;
    //a[i] = 0;
    //s[i] = l.getSample(gbwt::ENDMARKER, 0);
    //gbwt::CompressedRecord rec;

    //for (--i; i >= 0; --i)
    //    std::tie(a[i], s[i]) = safeLF(lfg, l, a[i+1], s[i+1], Qs[i+1], Qs[i]);
    //std::cout << "Offsets within nodes for Q*: " << a << '\n';
    //std::cout << "Suffixes within nodes for Q*: ";
    //if (s.size()){
    //    std::cout << r.unpack(*s.begin());
    //    std::for_each(++s.begin(), s.end(), [&r](const auto & x){std::cout << ',' << r.unpack(x);});
    //}
    //std::cout << '\n';
    //std::cout << "Suffixes (manual) within nodes for Q*: ";
    //if (s.size()){
        //gbwt::size_type suff = (r.decompressSA(Qs[0]).size() > *a.begin())?
                //r.decompressSA(Qs[0])[*a.begin()]:
                //gbwt::invalid_offset();
        //std::cout << r.unpack(suff);
        //for (unsigned i = 1; i < a.size(); ++i){
            //suff = (r.decompressSA(Qs[i]).size() > a[i])?
                //r.decompressSA(Qs[i])[a[i]]:
                //gbwt::invalid_offset();
            //std::cout << ',' << r.unpack(suff);
        //}
    //}
    //std::cout << '\n';

    //std::cout << "Computing lcps" << std::endl;
    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowLFGBWT(lfg, a, Qs);
    //gbwt::size_type computed_above, computed_below;
    //std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    //for (gbwt::size_type k = 0; k < Qs.size(); ++k){
    //    if (k == 0 || computed_above < k){
    //        gbwt::edge_type above = {Qs[k], a[k]-1};
    //        gbwt::size_type lcp = 0;
    //        computed_above = k;
    //        while (a[k+lcp] != 0 
    //                && above.first == Qs[k+lcp] 
    //                && above.first != gbwt::ENDMARKER 
    //                && above.second >= 0){
    //            if (above.second+1 == a[k+lcp])
    //                computed_above = k+lcp;
    //            ++lcp;
    //            above = lfg.inverseLF(above); 
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_above; ++kp)
    //            lcpa[kp] = lcp--;
    //    }
    //    if (k == 0 || computed_below < k){
    //        gbwt::edge_type below = {Qs[k], a[k]};
    //        gbwt::size_type lcp = 0;
    //        computed_below = k;
    //        while (a[k+lcp] != lfg.nodeSize(Qs[k+lcp]) 
    //                && below.first == Qs[k+lcp] 
    //                && below.first != gbwt::ENDMARKER 
    //                && below.second < lfg.nodeSize(below.first)){
    //            if (below.second == a[k+lcp])
    //                computed_below = k+lcp;
    //            ++lcp;
    //            below = lfg.inverseLF(below);
    //        }
    //        for (gbwt::size_type kp = k; kp <= computed_below; ++kp)
    //            lcpb[kp] = lcp--;
    //    }
    //}
    //std::cout << "lcpa: " << lcpa << std::endl << "lcpb: " << lcpb << std::endl;

    //std::cout << "Computing set maximal matches" << std::endl;
    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCPLFGBWT(lfg, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
            //gbwt::size_type curr = lcpa[k], plen,
            //    curr_s = (a[k] == lfg.nodeSize(Qs[k]))? 
            //        gbwt::invalid_offset():
            //        s[k];

            //while (curr == len){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //    curr_s = (curr_s == gbwt::invalid_offset())?
            //        l.locateLast(Qs[k]):
            //        l.locatePrev(curr_s);
            //    //std::cout << "prev, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs " ;
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    curr = std::min(curr, l.LCP(curr_s));
            //}
            //curr = lcpb[k];
            //curr_s = s[k];
            //while (curr == len){
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs ";
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    ++a[k];
            //    if (a[k] < lfg.nodeSize(Qs[k])){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //        curr_s = r.locateNext(curr_s);
            //    //std::cout << "next, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //        curr = std::min(curr, l.LCP(curr_s));
            //    } 
            //    else
            //        break;
            //}
        }
        prev = len;
    }
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery4(const lf_gbwt::GBWT & lfg, const gbwt::FastLocate & r, const FastLCP & l, const CompText & ct, const gbwt::vector_type & Q){
    //std::cout << "S4" << std::endl;
    for (gbwt::node_type a : Q)
        if (!lfg.contains(a) || lfg.nodeSize(a) == 0)
            return {};
    std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>> matches;
    if (!lfg.bidirectional()) { throw std::invalid_argument("Input GBWT must be bidirectional!"); }
    if (!Q.size()) { return matches; }

    //std::cout << "Query path provided: " << Q << "\n";
    gbwt::vector_type Qs = reverseAndAppendENDMARKER(Q);
    //gbwt::vector_type Qs = Q;
    ////reverse Q and append the endmarker
    //for (unsigned i = 0; i<Qs.size()/2; ++i)
    //    std::swap(Qs[i], Qs[Qs.size()-1-i]);
    //Qs.push_back(gbwt::ENDMARKER); 
    //std::cout << "Query path after reversal and endmarker appending: " << Qs << "\n";


    //std::cout << "Computing virtual insertion positions" << std::endl;
    //compute virtual insertion positions
    std::vector<gbwt::size_type> a, s;
    std::tie(a, s) = virtualInsertionWithSuffLFGBWT(lfg, r, l, Qs);
    //std::vector<gbwt::size_type> a(Qs.size());
    //std::vector<gbwt::size_type> s(Qs.size());
    //int i = Qs.size()-1;
    //a[i] = 0;
    //s[i] = l.getSample(gbwt::ENDMARKER, 0);
    //gbwt::CompressedRecord rec;

    //for (--i; i >= 0; --i)
    //    std::tie(a[i], s[i]) = safeLF(lfg, l, a[i+1], s[i+1], Qs[i+1], Qs[i]);
    //std::cout << "Offsets within nodes for Q*: " << a << std::endl;
    //std::cout << "Suffixes within nodes for Q*:          ";
    //if (s.size()){
    //    std::cout << r.unpack(*s.begin());
    //    std::for_each(++s.begin(), s.end(), [&r](const auto & x){std::cout << ',' << r.unpack(x);});
    //}
    //std::cout << std::endl;
    //std::cout << "Suffixes (manual) within nodes for Q*: ";
    //if (s.size()){
    //    gbwt::size_type suff = (r.decompressSA(Qs[0]).size() > *a.begin())?
    //            r.decompressSA(Qs[0])[*a.begin()]:
    //            gbwt::invalid_offset();
    //    std::cout << r.unpack(suff);
    //    for (unsigned i = 1; i < a.size(); ++i){
    //        suff = (r.decompressSA(Qs[i]).size() > a[i])?
    //            r.decompressSA(Qs[i])[a[i]]:
    //            gbwt::invalid_offset();
    //        std::cout << ',' << r.unpack(suff);
    //    }
    //}
    //std::cout << std::endl;
    //std::cout << "Above Suffixes within nodes for Q*:          " << std::flush;
    //if (s.size()){
    //    gbwt::size_type suff = (0 != *s.begin())?
    //            l.locatePrev(*s.begin()):
    //            gbwt::invalid_offset();
    //    std::cout << r.unpack(suff) << std::flush;
    //    for (unsigned i = 1; i < s.size(); ++i){
    //        suff = (0 != a[i])? ((a[i] == lfg.nodeSize(Qs[i]))? l.locateLast(Qs[i]) : l.locatePrev(s[i])):
    //            gbwt::invalid_offset();
    //        std::cout << ',' << r.unpack(suff) << std::flush;
    //    }
    //}
    //std::cout << std::endl;
    //std::cout << "Above Suffixes (manual) within nodes for Q*: ";
    //if (s.size()){
    //    gbwt::size_type suff = (r.decompressSA(Qs[0]).size() > *a.begin())?
    //            r.decompressSA(Qs[0])[(*a.begin())-1]:
    //            gbwt::invalid_offset();
    //    std::cout << r.unpack(suff);
    //    for (unsigned i = 1; i < a.size(); ++i){
    //        suff = (0 != a[i])?
    //            r.decompressSA(Qs[i])[a[i]-1]:
    //            gbwt::invalid_offset();
    //        std::cout << ',' << r.unpack(suff);
    //    }
    //}
    //std::cout << std::endl;

    //std::cout << "Computing lcps" << std::endl;
    //compute lcps
    std::vector<gbwt::size_type> lcpa, lcpb;
    std::tie(lcpa, lcpb) = lcpAboveBelowCT(lfg, l, ct, a, s, Qs);
    //std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    //for (gbwt::size_type k = 0; k < Qs.size(); ++k){
    //    //can add small optimization of only evaluate if not equal to previous suff
    //    gbwt::size_type lcp = 0;
    //    if (a[k] != 0){
    //        lcp = (k == 0 || lcpa[k-1] == 0)? 0 : (lcpa[k-1] - 1);
    //        gbwt::size_type suff = (a[k] == lfg.nodeSize(Qs[k]))? l.locateLast(Qs[k]) : l.locatePrev(s[k]);
    //        suff = ct.FLsuffToTrueSuff(suff);
    //        while (ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER)
    //            ++lcp;
    //    }
    //    lcpa[k] = lcp;
    //    lcp = 0;
    //    if (a[k] != lfg.nodeSize(Qs[k])){
    //        lcp = (k == 0 || lcpb[k-1] == 0)? 0 : (lcpb[k-1] - 1);
    //        gbwt::size_type suff = ct.FLsuffToTrueSuff(s[k]);
    //        while (ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER) 
    //            ++lcp;
    //    }
    //    lcpb[k] = lcp;
    //}
    //std::cout << "lcpa: " << lcpa << std::endl << "lcpb: " << lcpb << std::endl;

    //std::cout << "Computing set maximal matches" << std::endl;
    //output set maximal matches
    gbwt::size_type prev = 0, len;
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        len = std::max(lcpa[k], lcpb[k]);
        if (len >= prev){
            AddMatchesFastLCPLFGBWT(lfg, r, l, Q, Q.size()-k-len, len, a[k], s[k], lcpa[k], lcpb[k], matches);
            //gbwt::size_type curr = lcpa[k], plen,
            //    curr_s = (a[k] == lfg.nodeSize(Qs[k]))? 
            //        gbwt::invalid_offset():
            //        s[k];

            //while (curr == len){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //    curr_s = (curr_s == gbwt::invalid_offset())?
            //        l.locateLast(Qs[k]):
            //        l.locatePrev(curr_s);
            //    //std::cout << "prev, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs " ;
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    curr = std::min(curr, l.LCP(curr_s));
            //}
            //curr = lcpb[k];
            //curr_s = s[k];
            //while (curr == len){
            //    //std::cout << "pathLength Start " << std::flush;
            //    plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
            //    //std::cout << plen << " vs ";
            //    plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(curr_s)));
            //    //std::cout << plen;
            //    //std::cout << " End" << std::endl;
            //    matches.emplace_back(Q.size()-k-len, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
            //    ++a[k];
            //    if (a[k] < lfg.nodeSize(Qs[k])){
            //    //std::cout << "curr_s: " << r.unpack(curr_s) << " Start " << std::flush;
            //        curr_s = r.locateNext(curr_s);
            //    //std::cout << "next, curr_s: " << r.unpack(curr_s) << " End " << std::endl;
            //        curr = std::min(curr, l.LCP(curr_s));
            //    } 
            //    else
            //        break;
            //}
        }
        prev = len;
    }
    //std::cout << "done with matches" << std::endl;
    return matches;
}

std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>
setMaximalMatchQuery2_4(const gbwt::GBWT & x, const gbwt::FastLocate& r, const FastLCP& l, const CompText & ct, const gbwt::vector_type & Q){
    //std::cout << "S2_4" << std::endl;
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
