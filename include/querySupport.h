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

//locate is copied from the GBWT
/*
Copyright (c) 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024 Jouni Siren
Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#ifndef GBWT_QUERY_QUERYSUPPORT_H
#define GBWT_QUERY_QUERYSUPPORT_H

#include<vector>
#include<map>

#include<gbwt/gbwt.h>
#include<gbwt/fast_locate.h>
#include"fast_lcp.h"
#include"lf_gbwt.h"
#include"compText.h"

std::vector<gbwt::vector_type> extractPaths(const gbwt::GBWT&);
std::vector<gbwt::size_type> locate(const gbwt::GBWT&, const gbwt::SearchState&);
gbwt::vector_type reverseAndAppendENDMARKER(const gbwt::vector_type&);
std::vector<gbwt::size_type> virtualInsertionGBWT(const gbwt::GBWT&, const gbwt::vector_type&);
std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> virtualInsertionWithSuffGBWT(const gbwt::GBWT&, const gbwt::FastLocate&, const FastLCP&, const gbwt::vector_type&);
std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> virtualInsertionWithSuffLFGBWT(const lf_gbwt::GBWT&, const gbwt::FastLocate&, const FastLCP &, const gbwt::vector_type&);
std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> lcpAboveBelowGBWT(const gbwt::GBWT&, const std::vector<gbwt::size_type>&, const gbwt::vector_type&);
std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> lcpAboveBelowLFGBWT(const lf_gbwt::GBWT&, const std::vector<gbwt::size_type>&, const gbwt::vector_type&);
std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> lcpAboveBelowCT(const lf_gbwt::GBWT&, const FastLCP&, const CompText&, const std::vector<gbwt::size_type>&, const std::vector<gbwt::size_type>&, const gbwt::vector_type&);
void AddMatchesGBWT(const gbwt::GBWT&, const gbwt::vector_type&, gbwt::size_type, gbwt::size_type, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>&);
gbwt::size_type endmarkerSampleFastLocate(const gbwt::FastLocate&, gbwt::size_type);
void AddMatchesFastLocate(const gbwt::GBWT&, const gbwt::FastLocate&, const gbwt::vector_type&, gbwt::size_type, gbwt::size_type, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>&);
gbwt::size_type endmarkerSampleFastLCP(const FastLCP&, gbwt::size_type);
void AddMatchesFastLCP(const gbwt::GBWT&, const gbwt::FastLocate&, const FastLCP&, const gbwt::vector_type&, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>&);
void AddMatchesFastLCPLFGBWT(const lf_gbwt::GBWT&, const gbwt::FastLocate&, const FastLCP&, const gbwt::vector_type&, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>&);

std::vector<gbwt::vector_type>
extractPaths(const gbwt::GBWT & x){
    std::vector<gbwt::vector_type> result;
    for (gbwt::size_type i = 0; i < x.sequences(); ++i)
        result.push_back(x.extract(i));
    return result;
}

//copied from GBWT::locate, duplicates not removed
std::vector<gbwt::size_type>
locate(const gbwt::GBWT & x, const gbwt::SearchState& state)
{
    using namespace gbwt;
    std::vector<size_type> result;
    if(!(x.contains(state))) { return result; }

    // Initialize BWT positions for each offset in the range.
    std::vector<edge_type> positions(state.size());
    for(size_type i = state.range.first; i <= state.range.second; i++)
    {
        positions[i - state.range.first] = edge_type(state.node, i);
    }

    // Continue with LF() until samples have been found for all sequences.
    while(!(positions.empty()))
    {
        size_type tail = 0;
        node_type curr = invalid_node();
        CompressedRecord current;
        sample_type sample;
        edge_type LF_result;
        range_type LF_range;

        for(size_type i = 0; i < positions.size(); i++)
        {
            if(positions[i].first != curr)              // Node changed.
            {
                curr = positions[i].first; current = x.record(curr);
                sample = x.da_samples.nextSample(x.toComp(curr), positions[i].second);
                LF_range.first = positions[i].second;
                LF_result = current.runLF(positions[i].second, LF_range.second);
            }
            if(sample.first < positions[i].second)      // Went past the sample.
            {
                sample = x.da_samples.nextSample(x.toComp(curr), positions[i].second);
            }
            if(sample.first > positions[i].second)      // Not sampled, also valid for invalid_sample().
            {
                if(positions[i].second > LF_range.second) // Went past the existing LF() result.
                {
                    LF_range.first = positions[i].second;
                    LF_result = current.runLF(positions[i].second, LF_range.second);
                }
                positions[tail] = edge_type(LF_result.first, LF_result.second + positions[i].second - LF_range.first);
                tail++;
            }
            else                                        // Found a sample.
            {
                result.push_back(sample.second);
            }
        }
        positions.resize(tail);
        sequentialSort(positions.begin(), positions.end());
    }

    //removeDuplicates(result, false);
    return result;
}

//returns rev(q)$
gbwt::vector_type reverseAndAppendENDMARKER(const gbwt::vector_type & Q) {
    gbwt::vector_type Qs = Q; 
    Qs.resize(Q.size()+1, gbwt::ENDMARKER);
    for (gbwt::size_type i = 0; i < Q.size()/2; ++i)
        std::swap(Qs[i], Qs[Q.size()-1-i]);
    Qs.shrink_to_fit();
    return Qs;
}

//returns vector a where (Q[i], a[i]) is the virtual insertion position of the i-th suffix of Qs (Qs[i,Qs.size()-1])
//assumes Qs.back() = gbwt::ENDMARKER
std::vector<gbwt::size_type> virtualInsertionGBWT(const gbwt::GBWT& g, const gbwt::vector_type& Qs) {
    std::vector<gbwt::size_type> a (Qs.size());
    a.shrink_to_fit();
    gbwt::size_type pos;
    a.back() = 0;
    for (gbwt::size_type i = Qs.size() - 1; i != 0; --i) {
        pos = g.LF(Qs[i], a[i], Qs[i-1]);
        a[i-1] = (pos == gbwt::invalid_offset())? 0 : pos;
    }
    return a;
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>>
virtualInsertionWithSuffGBWT(const gbwt::GBWT& g, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::vector_type& Qs) {
    auto LFWithSuff = 
        [&g, &r, &l] (const gbwt::size_type prevPos, const gbwt::size_type prevSuff, const gbwt::node_type from, const gbwt::node_type to) -> std::pair<gbwt::size_type,gbwt::size_type> {
            //std::cout << "In LFWithSuff( prevPos: " << prevPos << ", prevSuff: " << prevSuff << ", from: " << from << ", to: " << to << ")" << std::endl;
            gbwt::CompressedRecord rec = g.record(from);

            gbwt::size_type outrank = rec.edgeTo(to);
            if (outrank >= rec.outdegree())
                return {0, r.locateFirst(to)};

            gbwt::CompressedRecordRankIterator iter(rec, outrank); 
            gbwt::size_type last_run_before = gbwt::invalid_offset();

            if(iter->first == outrank) { last_run_before = iter.runId(); }
            while(!(iter.end()) && iter.offset() < prevPos)
            {
                ++iter;
                if(iter->first == outrank)
                    last_run_before = iter.runId();
            }

            gbwt::size_type temp = iter.rank();
            if (prevPos < iter.offset() && iter->first == outrank)
                return {temp - (iter.offset() - prevPos), prevSuff - 1};
            
            if (last_run_before != gbwt::invalid_offset())
                return {temp, 
                    ((temp == g.nodeSize(to))? 
                     gbwt::invalid_offset() : 
                     r.locateNext(l.getSampleBot(from, last_run_before)-1))};

            while (iter->first != outrank)
                ++iter;
            return {temp, l.getSample(from, iter.runId())-1};
        };
    std::vector<gbwt::size_type> a(Qs.size());
    std::vector<gbwt::size_type> s(Qs.size());
    a.shrink_to_fit();
    s.shrink_to_fit();
    a.back() = 0;
    s.back() = r.locateFirst(gbwt::ENDMARKER);
    for (gbwt::size_type i = Qs.size() - 1; i != 0; --i)
        std::tie(a[i-1], s[i-1]) = LFWithSuff(a[i], s[i], Qs[i], Qs[i-1]);
    return {a, s};
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>>
virtualInsertionWithSuffLFGBWT(const lf_gbwt::GBWT & lfg, const gbwt::FastLocate& r, const FastLCP & l, const gbwt::vector_type& Qs) {
    auto LFWithSuff = 
        [&lfg, &r, &l] (const gbwt::size_type prevPos, const gbwt::size_type prevSuff, const gbwt::node_type from, const gbwt::node_type to) -> std::pair<gbwt::size_type,gbwt::size_type> {
            //std::cout << "In LFWithSuffLFGBWT( prevPos: " << prevPos << ", prevSuff: " << prevSuff << ", from: " << from << ", to: " << to << ")" << std::endl;
            auto ind = lfg.isSmallAndIndex(lfg.toComp(from));
            gbwt::comp_type compTo= lfg.toComp(to);
            //std::cout << "compTo: " << compTo << " ind: (" << ind.first << ", " << ind.second << ")" << std::endl;
            if (ind.first) {
                auto t = lfg.smallRecords.emptyAndNonEmptyIndex(ind.second);
                assert(!t.first); //since for every node in the query, there must be a path that contains it in the GBWT
                gbwt::size_type outrank = lfg.smallRecords.edgeTo(ind.second, compTo),
                    newPos = lfg.smallRecords.LF(ind.second, prevPos, compTo);
                //std::cout << "outrank " << outrank << " newPos " << newPos << std::endl;
                if (outrank >= lfg.smallRecords.outdegree(ind.second))
                    return {0, r.locateFirst(to)};
                //std::cout << "outrank not out of bounds" << std::endl;

                if (prevPos < lfg.smallRecords.size(ind.second) && lfg.smallRecords.bwtAt(ind.second, prevPos) == compTo)
                    return {newPos, prevSuff - 1};
                //std::cout << "prevPos bwt not equal to compTo" << std::endl;

                gbwt::size_type prefixLength = lfg.smallRecords.prefixSum.select_iter(t.second + 1)->second;
                auto nextRun = lfg.smallRecords.firstByAlphabet.successor(prefixLength*lfg.smallRecords.maxOutdegree + outrank* lfg.smallRecords.size(ind.second) + prevPos);
                if (nextRun->second - prefixLength*lfg.smallRecords.maxOutdegree < (outrank+1)*lfg.smallRecords.size(ind.second)) {
                    //successor run exists
                    //std::cout << "successor run exists" << std::endl;
                    gbwt::size_type run_id = lfg.smallRecords.logicalRunId(ind.second, nextRun->second - prefixLength*lfg.smallRecords.maxOutdegree - outrank*lfg.smallRecords.size(ind.second));
                    return {newPos, l.getSample(from, run_id) - 1};
                }
                else if (newPos < lfg.nodeSize(to)) {
                    //predecessor run exists
                    auto prevRun = --nextRun;
                    gbwt::size_type run_id = lfg.smallRecords.logicalRunId(ind.second, prevRun->second - prefixLength*lfg.smallRecords.maxOutdegree - outrank*lfg.smallRecords.size(ind.second));
                    return {newPos, r.locateNext(l.getSampleBot(from, run_id) - 1)};
                }
                return {newPos, gbwt::invalid_offset()};
            }
            else {
                const lf_gbwt::CompressedRecord& rec = lfg.largeRecords[ind.second];

                gbwt::size_type outrank = rec.edgeTo(compTo),
                    newPos = rec.LF(prevPos, compTo);
                if (outrank >= rec.outdegree())
                    return {0, r.locateFirst(to)};

                if (prevPos < rec.size() && rec[prevPos] == compTo)
                    return {newPos, prevSuff - 1};

                auto nextRun = rec.firstByAlphabet.successor(outrank*rec.size() + prevPos);
                if (nextRun->second < (outrank+1)*rec.size()){
                    //successor run exists
                    gbwt::size_type run_id = rec.logicalRunId(nextRun->second-outrank*rec.size());
                    return {newPos, l.getSample(from, run_id) - 1};
                }
                else if (newPos < lfg.nodeSize(to)){
                    //predecessor run exists
                    auto prevRun = --nextRun;
                    gbwt::size_type run_id = rec.logicalRunId(prevRun->second - outrank*rec.size());
                    return {newPos, r.locateNext(l.getSampleBot(from, run_id) - 1)};
                }
                //predecessor run exists but newPos == lfg.record(to).size
                return {newPos, gbwt::invalid_offset()};
            }
        };
    std::vector<gbwt::size_type> a(Qs.size());
    std::vector<gbwt::size_type> s(Qs.size());
    a.shrink_to_fit();
    s.shrink_to_fit();
    a.back() = 0;
    s.back() = r.locateFirst(gbwt::ENDMARKER);
    for (gbwt::size_type i = Qs.size() - 1; i != 0; --i)
        std::tie(a[i-1], s[i-1]) = LFWithSuff(a[i], s[i], Qs[i], Qs[i-1]);
    return {a, s};
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>>
lcpAboveBelowGBWT(const gbwt::GBWT& g, const std::vector<gbwt::size_type>& a, const gbwt::vector_type& Qs) {
    //length of lcpa and lcpb that has been computed
    gbwt::size_type computed_above = 0, computed_below = 0;
    std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    lcpa.shrink_to_fit();
    lcpb.shrink_to_fit();
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        if (computed_above <= k) {
            //compute lcpa[k]
            gbwt::size_type lcp = 0;
            computed_above = k + 1;
            if (a[k] != 0) {
                gbwt::edge_type above = {Qs[k], a[k]-1};
                while (above.first == Qs[k+lcp] && above.first != gbwt::ENDMARKER) {
                    if (above.second + 1 == a[k+lcp]) { computed_above = k + lcp + 1; }
                    ++lcp;
                    above = g.inverseLF(above);
                }
            }
            for (gbwt::size_type kp = k; kp < computed_above; ++kp, --lcp)
                lcpa[kp] = lcp;
        }
        if (computed_below <= k) {
            //compute lcpb[k]
            gbwt::size_type lcp = 0;
            computed_below = k + 1;
            if (a[k] != g.nodeSize(Qs[k])) {
                gbwt::edge_type below = {Qs[k], a[k]};
                while (below.first == Qs[k+lcp] && below.first != gbwt::ENDMARKER) {
                    if (below.second == a[k+lcp]) { computed_below = k + lcp + 1; }
                    ++lcp;
                    below = g.inverseLF(below);
                }
            }
            for (gbwt::size_type kp = k; kp < computed_below; ++kp, --lcp)
                lcpb[kp] = lcp;
        }
    }
    return {lcpa, lcpb};
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>>
lcpAboveBelowLFGBWT(const lf_gbwt::GBWT& lfg, const std::vector<gbwt::size_type>& a, const gbwt::vector_type& Qs) {
    //length of lcpa and lcpb that has been computed
    gbwt::size_type computed_above = 0, computed_below = 0;
    std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    lcpa.shrink_to_fit();
    lcpb.shrink_to_fit();
    for (gbwt::size_type k = 0; k < Qs.size(); ++k){
        if (computed_above <= k) {
            //compute lcpa[k]
            gbwt::size_type lcp = 0;
            computed_above = k + 1;
            if (a[k] != 0) {
                gbwt::edge_type above = {Qs[k], a[k]-1};
                while (above.first == Qs[k+lcp] && above.first != gbwt::ENDMARKER) {
                    if (above.second + 1 == a[k+lcp]) { computed_above = k + lcp + 1; }
                    ++lcp;
                    above = lfg.inverseLF(above);
                }
            }
            for (gbwt::size_type kp = k; kp < computed_above; ++kp, --lcp)
                lcpa[kp] = lcp;
        }
        if (computed_below <= k) {
            //compute lcpb[k]
            gbwt::size_type lcp = 0;
            computed_below = k + 1;
            if (a[k] != lfg.nodeSize(Qs[k])) {
                gbwt::edge_type below = {Qs[k], a[k]};
                while (below.first == Qs[k+lcp] && below.first != gbwt::ENDMARKER) {
                    if (below.second == a[k+lcp]) { computed_below = k + lcp + 1; }
                    ++lcp;
                    below = lfg.inverseLF(below);
                }
            }
            for (gbwt::size_type kp = k; kp < computed_below; ++kp, --lcp)
                lcpb[kp] = lcp;
        }
    }
    return {lcpa, lcpb};
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> 
lcpAboveBelowCT(const lf_gbwt::GBWT& lfg, const FastLCP& l, const CompText& ct, const std::vector<gbwt::size_type>& a, const std::vector<gbwt::size_type>& s, const gbwt::vector_type& Qs) {
    std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    lcpa.shrink_to_fit();
    lcpb.shrink_to_fit();
    for (gbwt::size_type k = 0; k < Qs.size(); ++k) {
        gbwt::size_type lcp = 0;
        if (a[k] != 0) {
            lcp = (k == 0 || lcpa[k-1] == 0)? 0 : (lcpa[k-1] - 1);
            gbwt::size_type suff = (a[k] == lfg.nodeSize(Qs[k]))? l.locateLast(Qs[k]) : l.locatePrev(s[k]);
            suff = ct.FLsuffToTrueSuff(suff);
            while(ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER)
                ++lcp;
        }
        lcpa[k] = lcp;
        lcp = 0;
        if (a[k] != lfg.nodeSize(Qs[k])){
            lcp = (k == 0 || lcpb[k-1] == 0)? 0 : (lcpb[k-1] - 1);
            gbwt::size_type suff = ct.FLsuffToTrueSuff(s[k]);
            while (ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER)
                ++lcp;
        }
        lcpb[k] = lcp;
    }
    return {lcpa, lcpb};
}

std::pair<std::vector<gbwt::size_type>,std::vector<gbwt::size_type>> 
lcpAboveBelowCT(const gbwt::GBWT& g, const FastLCP& l, const CompText& ct, const std::vector<gbwt::size_type>& a, const std::vector<gbwt::size_type>& s, const gbwt::vector_type& Qs) {
    std::vector<gbwt::size_type> lcpa(Qs.size()), lcpb(Qs.size());
    lcpa.shrink_to_fit();
    lcpb.shrink_to_fit();
    for (gbwt::size_type k = 0; k < Qs.size(); ++k) {
        gbwt::size_type lcp = 0;
        if (a[k] != 0) {
            lcp = (k == 0 || lcpa[k-1] == 0)? 0 : (lcpa[k-1] - 1);
            gbwt::size_type suff = (a[k] == g.nodeSize(Qs[k]))? l.locateLast(Qs[k]) : l.locatePrev(s[k]);
            suff = ct.FLsuffToTrueSuff(suff);
            while(ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER)
                ++lcp;
        }
        lcpa[k] = lcp;
        lcp = 0;
        if (a[k] != g.nodeSize(Qs[k])){
            lcp = (k == 0 || lcpb[k-1] == 0)? 0 : (lcpb[k-1] - 1);
            gbwt::size_type suff = ct.FLsuffToTrueSuff(s[k]);
            while (ct.at(suff+lcp) == Qs[k+lcp] && Qs[k+lcp] != gbwt::ENDMARKER)
                ++lcp;
        }
        lcpb[k] = lcp;
    }
    return {lcpa, lcpb};
}

void AddMatchesGBWT(const gbwt::GBWT& x, const gbwt::vector_type& Q, gbwt::size_type k,
        gbwt::size_type len, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches){
    //std::cout << "AddMatchesGBWT called: " << k << ' ' << len << std::endl;
    std::vector<gbwt::size_type> s = locate(x, gbwt::find(x, Q.begin()+k, Q.begin()+k+len));
    for (unsigned i = 0; i<s.size(); ++i)
        matches.emplace_back(k, len, s[i]);
}

//outputs length of path i (excluding endmarker)
gbwt::size_type endmarkerSampleFastLocate(const gbwt::FastLocate& r, gbwt::size_type i) {
    if (r.comp_to_run.size()<=1) 
        throw std::invalid_argument("endmarkerSampleFastLocate called on r-index with no nodes (excluding the endmaker)"); //return 0;
    gbwt::size_type max = std::min(r.comp_to_run[1]-1, i);
    if (r.seqId(r.samples[max]) == i) //common case in large bidirectional gbwts
        return r.samples[max]; 

    gbwt::size_type first = *(--std::lower_bound(r.samples.begin(), r.samples.begin()+max+1, r.pack(i+1,0)));
    while (r.seqId(first) != i){
        first = r.locateNext(first);
    }
    return first;
}

void AddMatchesFastLocate(const gbwt::GBWT& x, const gbwt::FastLocate& r, const gbwt::vector_type& Q, 
        gbwt::size_type k, gbwt::size_type len, 
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches){
    gbwt::size_type first = 0;
    gbwt::SearchState ss = r.find(Q.begin()+k, Q.begin()+k+len, first);
    gbwt::size_type plen;
    for (unsigned i = 0; i<ss.size(); ++i){
        plen = r.seqOffset(endmarkerSampleFastLocate(r, r.seqId(first))); 
        matches.emplace_back(k, len, r.seqId(first), plen - r.seqOffset(first) - len);
        first = r.locateNext(first);
    }
}

gbwt::size_type endmarkerSampleFastLCP(const FastLCP& l, gbwt::size_type i){
    auto iter = l.first.predecessor(l.rindex->pack(i+1, 0)-1);
    return iter->second + l.samples_lcp[iter->first];
}

void AddMatchesFastLCP(const gbwt::GBWT& g, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::vector_type& Q,
        gbwt::size_type k, gbwt::size_type len, 
        gbwt::size_type pos, gbwt::size_type suff, gbwt::size_type lcpa, gbwt::size_type lcpb,
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    //do above matches
    gbwt::size_type curr_s = suff, plen;
    while (lcpa == len) {
        curr_s = (curr_s == suff && pos == g.nodeSize(Q[k+len-1]))? l.locateLast(Q[k+len-1]) : l.locatePrev(curr_s);
        plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
        matches.emplace_back(k, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
        lcpa = std::min(lcpa, l.LCP(curr_s));
    }
    //do below matches
    curr_s = suff;
    while (lcpb == len) {
        plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
        matches.emplace_back(k, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
        ++pos;
        if (pos == g.nodeSize(Q[k+len-1])) { break; }
        curr_s = r.locateNext(curr_s);
        lcpb = std::min(lcpb, l.LCP(curr_s));
    }
}

void AddMatchesFastLCPLFGBWT(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, const gbwt::vector_type& Q,
        gbwt::size_type k, gbwt::size_type len,
        gbwt::size_type pos, gbwt::size_type suff, gbwt::size_type lcpa, gbwt::size_type lcpb,
        std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    //std::cout << "In AddMatchesFastLCPLFGBWT" << std::endl;
    //do above matches
    gbwt::size_type curr_s = suff, plen;
    while (lcpa == len) {
        curr_s = (curr_s == suff && pos == lfg.nodeSize(Q[k+len-1]))? l.locateLast(Q[k+len-1]) : l.locatePrev(curr_s);
        plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
        matches.emplace_back(k, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
        lcpa = std::min(lcpa, l.LCP(curr_s));
    }
    //std::cout << "Done above matches" << std::endl;
    //do below matches
    curr_s = suff;
    while (lcpb == len) {
        plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(curr_s)));
        matches.emplace_back(k, len, r.seqId(curr_s), plen - r.seqOffset(curr_s) - len);
        ++pos;
        if (pos == lfg.nodeSize(Q[k+len-1])) { break; }
        curr_s = r.locateNext(curr_s);
        lcpb = std::min(lcpb, l.LCP(curr_s));
    }
    //std::cout << "Leaving AddMatchesFastLCPLFGBWT" << std::endl;
}


void AddLongMatchesWholeBlock(const gbwt::FastLocate& r, const FastLCP& l, std::map<gbwt::size_type, gbwt::size_type>& inBlock, const gbwt::size_type currQsInd, const gbwt::vector_type& Qs, 
        const gbwt::range_type block, const gbwt::size_type topSuff, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    //std::cout << "In AddLongMatchesWholeBlock(currQsInd " << currQsInd << ", block [" << block.first << ", " << block.second << "], topSuff " << topSuff << ")" << std::endl;
    gbwt::size_type suff;
    for (gbwt::size_type blockInd = block.first; blockInd != block.second + 1; ++blockInd){
        suff = (blockInd == block.first)? topSuff : r.locateNext(suff);
        auto it = inBlock.find(suff + Qs.size() - currQsInd);
        assert(it != inBlock.end());
        gbwt::size_type plen = r.seqOffset(endmarkerSampleFastLCP(l, r.seqId(suff)));
        matches.emplace_back(Qs.size() - 1 - it->second, it->second - currQsInd, r.seqId(suff), plen - r.seqOffset(suff) - (it->second - currQsInd));
        inBlock.erase(it);
    }
}
        

//output suffixes in block that are leaving the block from Qs[currQsInd] to Qs[currQsInd-1], assumes 0<currQsInd<currQsInd.size()
//returns block i-1 of length L+1 and top and bottom suffixes of this block (inclusive)
//if block is empty, top and bot suffix are invalid and returned block is arbitrary empty block
//assumes block is valid
std::tuple<gbwt::range_type,gbwt::size_type,gbwt::size_type>
AddLongMatchesFastLCP(const gbwt::GBWT& g, const gbwt::FastLocate& r, const FastLCP& l, std::map<gbwt::size_type, gbwt::size_type>& inBlock, const gbwt::size_type currQsInd, 
        const gbwt::vector_type& Qs, const gbwt::range_type block, const gbwt::size_type topSuff, const gbwt::size_type botSuff, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    if (gbwt::Range::empty(block)) { return {gbwt::Range::empty_range(), gbwt::invalid_offset(), gbwt::invalid_offset()}; }
    assert(currQsInd != 0);
    assert(inBlock.size() == gbwt::Range::length(block));
    assert(block.second < g.nodeSize(Qs[currQsInd]) && block.first <= block.second);
    assert(topSuff < r.pack(g.sequences(), 0) && botSuff < r.pack(g.sequences(), 0));

    gbwt::node_type from = Qs[currQsInd], to = Qs[currQsInd-1];
    gbwt::CompressedRecord rec = g.record(from);
    //returns the runId of the logical run pos is part of
    auto runIdPos = [&rec] (gbwt::CompressedRecordRankIterator& iter, const gbwt::size_type pos) -> gbwt::size_type {
        assert(pos < iter.offset() && pos >= iter.offset() - iter->second);
        if (iter->first != 0 || rec.outgoing[0].first != gbwt::ENDMARKER) { return iter.runId(); }
        return iter.runId() - (iter.offset() - 1 - pos);

    };
    //returns runId of run, if run of endmarkers returns runId of first run
    auto firstRunId = [&rec,&runIdPos] (gbwt::CompressedRecordRankIterator& iter) -> gbwt::size_type {
        return runIdPos(iter, iter.offset() - iter->second);
    };
    //advances iter to run that contains pos
    //returns number of outranks in [0,pos) and if bwt[pos] == outrank
    //pos < node.size()
    auto LFWithIter = [] (gbwt::CompressedRecordRankIterator& iter, const gbwt::size_type pos) -> std::pair<gbwt::size_type,bool> {
        //assert first position of iter is not past pos
        assert(iter.offset() - iter->second <= pos);
        iter.readPast(pos);
        assert(!iter.end());
        gbwt::size_type ans = iter.rank();
        if (iter->first == iter.rank_support.value) { return { ans - (iter.offset() - pos), true }; }
        return { ans, false };
    };

    gbwt::range_type newBlock = gbwt::Range::empty_range();
    gbwt::size_type outrank = rec.edgeTo(to);

    //compute new Block
    bool topSuffContinue, botSuffContinue;
    gbwt::CompressedRecordRankIterator topIter(rec, outrank), botIter(rec, outrank);
    if (outrank < rec.outdegree()) {
        std::tie(newBlock.first , topSuffContinue) = LFWithIter(topIter, block.first);
        std::tie(newBlock.second, botSuffContinue) = LFWithIter(botIter, block.second);
        newBlock.second -= !botSuffContinue;
    }
    if (gbwt::Range::empty(newBlock)) {
        AddLongMatchesWholeBlock(r, l, inBlock, currQsInd, Qs, block, topSuff, matches);
        return {gbwt::Range::empty_range(), gbwt::invalid_offset(), gbwt::invalid_offset()};
    }

    assert(topIter.runId() <= botIter.runId());

    //make start and end handle endmarkers
    gbwt::size_type firstOutrankRunId = gbwt::invalid_offset(), lastOutrankRunId = gbwt::invalid_offset();
    while (true) {
        if (topIter->first == outrank) {
            lastOutrankRunId = runIdPos(topIter, std::min(block.second, topIter.offset() - 1));
            if (firstOutrankRunId == gbwt::invalid_offset()) { firstOutrankRunId = firstRunId(topIter); }
        }
        else { 
            gbwt::range_type runBlock = {std::max(block.first, topIter.offset() - topIter->second), std::min(block.second, topIter.offset() - 1)};
            gbwt::size_type firstSuff = (block.first == runBlock.first)? topSuff : l.getSample(from, firstRunId(topIter));
            AddLongMatchesWholeBlock(r, l, inBlock, currQsInd, Qs, runBlock, firstSuff, matches); 
        }
        if (topIter.runId() == botIter.runId()) { break; }
        ++topIter;
    } 
    assert((firstOutrankRunId) != gbwt::invalid_offset() && lastOutrankRunId != gbwt::invalid_offset());

    return {newBlock, ((topSuffContinue)? topSuff - 1 : l.getSample(from, firstOutrankRunId) - 1), ((botSuffContinue)? botSuff - 1 : l.getSampleBot(from, lastOutrankRunId) - 1)};
}

std::tuple<gbwt::range_type,gbwt::size_type,gbwt::size_type>
LongMatchesFastLCPInitializeEmptyBlock(const gbwt::FastLocate& r, const FastLCP& l, std::map<gbwt::size_type,gbwt::size_type>& inBlock, 
        const gbwt::size_type nodePos, const gbwt::node_type node, const gbwt::size_type nodeSize, const gbwt::size_type bSuff, const gbwt::size_type lcpa, const gbwt::size_type lcpb, 
        const gbwt::size_type queryLen, const gbwt::size_type queryPos, const gbwt::size_type L) {
    //std::cout << "LongMatchesFastLCPInitializeEmptyBlock(nodePos " << nodePos << ", node " << node << ", nodeSize " << nodeSize << ", bSuff " << bSuff << ", lcpa " << lcpa << ", lcpb " << lcpb 
        //<< ", queryLen " << queryLen << ", " << queryPos << ", L " << L << ")" << std::endl;
    assert(nodePos <= nodeSize);
    assert((lcpa == 0) == (nodePos == 0));
    assert((lcpb == 0) == (nodePos == nodeSize));
    assert(lcpa == L || lcpb == L);
    assert(lcpa <= L && lcpb <= L);
    gbwt::range_type block;
    gbwt::size_type suff;
    if (lcpb == L) {
        assert(lcpb == L);
        assert(nodePos < nodeSize);
        block = {nodePos, nodePos};
        suff = bSuff;
    }
    else {
        assert(lcpa == L);
        block = {nodePos-1, nodePos-1};
        suff = (nodePos == nodeSize)? l.locateLast(node) : l.locatePrev(bSuff);
    }
    inBlock.emplace(suff + queryLen - queryPos, queryPos + L);
    return {block, suff, suff};
}

//given [f^{L+1}_{queryPos}, g^{L+1}_{queryPos} - 1] 
//if [f^{L+1}_{queryPos}, g^{L+1}_{queryPos} - 1] and output block is nonempty, provided block is [a,a]
//where the LCP of the suffix at a and Q[queryPos...] is length L and is a longest LCP
//compute nonempty block [f^L_{queryPos}, g^L_{queryPos} - 1]
std::tuple<gbwt::range_type,gbwt::size_type,gbwt::size_type>
LongMatchesExpandNonEmptyBlock(const gbwt::FastLocate& r, const FastLCP& l, std::map<gbwt::size_type,gbwt::size_type>& inBlock,
        const gbwt::size_type nodeSize, const gbwt::range_type block, const gbwt::size_type topSuff, const gbwt::size_type botSuff,
        const gbwt::size_type queryLen, const gbwt::size_type queryPos, const gbwt::size_type L) {
    //std::cout << "LongMatchesExpandNonEmptyBlock(nodeSize " << nodeSize << ", block [" << block.first << "," << block.second << "], topSuff " << topSuff << ", botSuff " << botSuff << ", queryLen " << queryLen
        //<< ", queryPos " << queryPos << ", L " << L << ")" << std::endl;
    assert(block.first < nodeSize && block.second < nodeSize && block.second >= block.first);
    gbwt::range_type newBlock = block;
    gbwt::size_type newTopSuff = topSuff, newBotSuff = botSuff, nextSuff;
    bool succeed;
    while (l.LCP(newTopSuff) >= L) {
        assert(newBlock.first != 0);
        newTopSuff = l.locatePrev(newTopSuff);
        --newBlock.first;
        succeed = inBlock.emplace(newTopSuff + queryLen - queryPos, queryPos + L).second;
        assert(succeed);
    }
    while (newBlock.second + 1 < nodeSize) {
        nextSuff = r.locateNext(newBotSuff);
        if (l.LCP(nextSuff) < L) { break; }
        ++newBlock.second;
        newBotSuff = nextSuff;
        succeed = inBlock.emplace(newBotSuff + queryLen - queryPos, queryPos + L).second;
        assert(succeed);
    }
    return {newBlock, newTopSuff, newBotSuff};
}

//output suffixes in block that are leaving the block from Qs[currQsInd] to Qs[currQsInd-1], assumes 0<currQsInd<currQsInd.size()
//returns block i-1 of length L+1 and top and bottom suffixes of this block (inclusive)
//if block is empty, top and bot suffix are invalid and returned block is arbitrary empty block
//assumes block is valid
std::tuple<gbwt::range_type,gbwt::size_type,gbwt::size_type>
AddLongMatchesLFGBWT(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r, const FastLCP& l, std::map<gbwt::size_type, gbwt::size_type>& inBlock, const gbwt::size_type currQsInd, 
        const gbwt::vector_type& Qs, const gbwt::range_type block, const gbwt::size_type topSuff, const gbwt::size_type botSuff, std::vector<std::tuple<gbwt::size_type,gbwt::size_type,gbwt::size_type,gbwt::size_type>>& matches) {
    if (gbwt::Range::empty(block)) { return {gbwt::Range::empty_range(), gbwt::invalid_offset(), gbwt::invalid_offset()}; }
    assert(currQsInd != 0);
    assert(inBlock.size() == gbwt::Range::length(block));
    assert(block.second < lfg.nodeSize(Qs[currQsInd]) && block.first <= block.second);
    assert(topSuff < r.pack(lfg.sequences(), 0) && botSuff < r.pack(lfg.sequences(), 0));

    gbwt::comp_type compTo = lfg.toComp(Qs[currQsInd-1]);
    auto ind = lfg.isSmallAndIndex(lfg.toComp(Qs[currQsInd]));
    auto t = (ind.first)?lfg.smallRecords.emptyAndNonEmptyIndex(ind.second) : std::pair<bool, gbwt::size_type>{lfg.largeRecords[ind.second].size() == 0, ind.second};
    assert(!t.first); //since for every node in the query, there must be a path in the GBWT that contains it
    //const lf_gbwt::CompressedRecord& rec = lfg.record(Qs[currQsInd]);

    gbwt::range_type newBlock = gbwt::Range::empty_range();
    gbwt::size_type outrank = (ind.first)? lfg.smallRecords.edgeTo(ind.second, compTo) : lfg.largeRecords[ind.second].edgeTo(compTo);
    //std::cout << "outrank " << outrank << std::endl;
    if (outrank != gbwt::invalid_offset()) {
        newBlock.first  = (ind.first)? lfg.smallRecords.LF(ind.second, block.first, compTo) : lfg.largeRecords[ind.second].LF(block.first, compTo);
        newBlock.second = ((ind.first)? lfg.smallRecords.LF(ind.second, block.second + 1, compTo) : lfg.largeRecords[ind.second].LF(block.second + 1, compTo)) - 1;
    }
    //std::cout << "newBlock: [" << newBlock.first << "," << newBlock.second << "]" << std::endl;
    if (gbwt::Range::empty(newBlock)) {
        AddLongMatchesWholeBlock(r, l, inBlock, currQsInd, Qs, block, topSuff, matches);
        return {gbwt::Range::empty_range(), gbwt::invalid_offset(), gbwt::invalid_offset()};
    }

    gbwt::size_type prefixSum = (ind.first)? lfg.smallRecords.prefixSum.select_iter(t.second + 1)->second : 0;
    sdsl::sd_vector<>::one_iterator start = ((ind.first)? lfg.smallRecords.first.predecessor(prefixSum + block.first) : lfg.largeRecords[ind.second].first.predecessor(block.first)), 
        end = ((ind.first)? lfg.smallRecords.first.successor(prefixSum + block.second + 1) : lfg.largeRecords[ind.second].first.successor(block.second + 1));
    //std::cout << "start->second - prefixSum" << start->second - prefixSum 
        //<< "end->second - prefixSum" << end->second - prefixSum << std::endl;
    for (auto next = start; start != end; ++start) {
        //std::cout << "checking run" << std::endl;
        ++next; 
        gbwt::comp_type runVal = (ind.first)? lfg.smallRecords.alphabetByRun[start->first] : lfg.largeRecords[ind.second].alphabetByRun[start->first];
        if (runVal == outrank) { continue; }
        gbwt::range_type runBlock = {std::max(block.first, start->second - prefixSum), std::min(block.second, next->second - prefixSum - 1)};
        //std::cout << "adding run [" << runBlock.first << ", " << runBlock.second << "]" << std::endl;
        gbwt::size_type firstSuff = (block.first == runBlock.first)? topSuff : 
            l.getSample(Qs[currQsInd], 
                    ((ind.first)? lfg.smallRecords.logicalRunId(ind.second, start->second - prefixSum): lfg.largeRecords[ind.second].logicalRunId(start->second))
                    );
        //if (ind.first)
            //std::cout << "lfg.smallRecords.logicalRunId(ind.second, start->second - prefixSum) " << lfg.smallRecords.logicalRunId(ind.second, start->second - prefixSum) << std::endl;
        AddLongMatchesWholeBlock(r, l, inBlock, currQsInd, Qs, runBlock, firstSuff, matches);
    }

    gbwt::size_type newTopSuff, newBotSuff;
    gbwt::comp_type bwtValTop = (ind.first)? lfg.smallRecords.compAlphabetAt(ind.second, block.first) : lfg.largeRecords[ind.second].compAlphabetAt(block.first);
    gbwt::comp_type bwtValBot = (ind.first)? lfg.smallRecords.compAlphabetAt(ind.second, block.second): lfg.largeRecords[ind.second].compAlphabetAt(block.second);
    if (bwtValTop == outrank) { newTopSuff = topSuff - 1; }
    else { 
        gbwt::size_type nextOutrankRunStart = (ind.first)
            ? lfg.smallRecords.firstByAlphabet.successor(prefixSum*lfg.smallRecords.maxOutdegree + outrank*lfg.smallRecords.size(ind.second) + block.first)->second - prefixSum*lfg.smallRecords.maxOutdegree - outrank*lfg.smallRecords.size(ind.second)
            : lfg.largeRecords[ind.second].firstByAlphabet.successor(outrank*lfg.largeRecords[ind.second].size() + block.first )->second - outrank*lfg.largeRecords[ind.second].size();
        assert(nextOutrankRunStart <= block.second);
        newTopSuff = l.getSample(Qs[currQsInd], 
                ((ind.first)? lfg.smallRecords.logicalRunId(ind.second, nextOutrankRunStart): lfg.largeRecords[ind.second].logicalRunId(nextOutrankRunStart))
                ) - 1;
    }
    if (bwtValBot == outrank) {newBotSuff = botSuff - 1; }
    else {
        gbwt::size_type lastOutrankRunStart = (ind.first)
            ? lfg.smallRecords.firstByAlphabet.predecessor(prefixSum*lfg.smallRecords.maxOutdegree + outrank*lfg.smallRecords.size(ind.second) + block.second)->second - prefixSum*lfg.smallRecords.maxOutdegree - outrank*lfg.smallRecords.size(ind.second)
            : lfg.largeRecords[ind.second].firstByAlphabet.predecessor(outrank*lfg.largeRecords[ind.second].size() + block.second)->second - outrank*lfg.largeRecords[ind.second].size();
        gbwt::size_type lastOutrankRunEnd = (ind.first)
            ? lfg.smallRecords.first.successor(prefixSum + lastOutrankRunStart + 1)->second - prefixSum - 1
            : lfg.largeRecords[ind.second].first.successor(lastOutrankRunStart + 1)->second - 1;
        newBotSuff = l.getSampleBot(Qs[currQsInd], 
                ((ind.first)? lfg.smallRecords.logicalRunId(ind.second, lastOutrankRunEnd) : lfg.largeRecords[ind.second].logicalRunId(lastOutrankRunEnd))
                ) - 1;
    }

    return {newBlock, newTopSuff, newBotSuff};
}
#endif //GBWT_QUERY_QUERYSUPPORT_H
