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

#ifndef GBWT_QUERY_IOHELP_H
#define GBWT_QUERY_IOHELP_H

#include<iostream>
#include<gbwt/gbwt.h>
#include<gbwt/fast_locate.h>
#include<gbwt/internal.h>
#include<vector>
#include<iomanip>
#include"fast_lcp.h"
#include"lf_gbwt.h"
#include<sstream>
#include"compText.h"

template<class T1, class T2, class T3>
std::ostream& operator<<(std::ostream& out, const std::tuple<T1,T2,T3> & a){
    out << '(' << std::get<0>(a) << ',' << std::get<1>(a) << ',' << std::get<2>(a) << ')';
    return out;
}

template<class T1, class T2, class T3, class T4>
std::ostream& operator<<(std::ostream& out, const std::tuple<T1,T2,T3,T4> & a){
    out << '(' << std::get<0>(a) << ',' << std::get<1>(a) << ',' << std::get<2>(a) << ',' << std::get<3>(a) << ')';
    return out;
}

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T> & a){
    if (a.size()){
        out << *a.begin();
        std::for_each(++a.begin(), a.end(), [&out](const T & x){out << ',' << x;});
    }
    return out;
}

template<class T1, class T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1,T2> & p){
    out << '(' << p.first << ',' << p.second << ')';
    return out;
}

template<class T1, class T2>
std::string areEqual(const T1 & a, const T2 & b){
    std::stringstream out;
    out << a;
    if (a != b)
        out << " not equal! " << b;
    return out.str();
}

void printGBWTandRindex(std::ostream& out, const gbwt::GBWT & g, const gbwt::FastLocate & r, const FastLCP & l, const lf_gbwt::GBWT & lfg, const CompText & ct, int width){
    using std::setw;
    int labelWidth = 3*width+1;
    out << "Printing GBWT, lf GBWT, and r-index:" << '\n';
    for (gbwt::node_type i = 0; i < g.effective(); ++i){
        gbwt::node_type node = g.toNode(i);
        gbwt::size_type size = g.nodeSize(node);
        gbwt::CompressedRecord a = g.record(node); 
        const lf_gbwt::CompressedRecord & lfg_a = lfg.record(node);
        out << "Printing Node " << setw(width) << areEqual(node, lfg.toNode(i))
            << setw(labelWidth) << "Comp:" << setw(width) << i 
            << setw(labelWidth) << "Size:" << setw(width) << areEqual(size, lfg.nodeSize(node))
            << setw(labelWidth) << "ConcRuns:" << setw(width) << areEqual(a.runs().first, lfg_a.runs().first)
            << setw(labelWidth) << "LogiRuns:" << setw(width) << areEqual(a.runs().second, lfg_a.runs().second) <<'\n';
        out << std::flush;
        std::vector<gbwt::edge_type> outgoing = g.edges(node);



        std::vector<gbwt::size_type> suffixA, suffixB, LCPA, LCPB;
        suffixA.resize(size);
        suffixB.resize(size);
        LCPA.resize(size);
        LCPB.resize(size);
        suffixA[0] = r.locateFirst(node);
        suffixB.back() = l.locateLast(node);
        LCPA[0] = l.LCP(suffixA[0]);
        LCPB.back() = l.LCP(suffixB.back());


        for (unsigned j = 1; j < size; ++j){
            suffixA[j] = r.locateNext(suffixA[j-1]);
            suffixB[size-1-j] = l.locatePrev(suffixB[size-j]);
            LCPA[j] = l.LCP(suffixA[j]);
            LCPB[size-1-j] = l.LCP(suffixB[size-1-j]);
        }

        out << setw(labelWidth) << "suff (FLoc):";
        for (unsigned j = 0; j < size; ++j){
            out << setw(width) << areEqual(r.unpack(suffixA[j]), r.unpack(suffixB[j])) << ' ';
        }
        out << std::flush;
        out << "\n" << setw(labelWidth) << " LCP A (FLoc):";
        for (unsigned j = 0; j < size; ++j){
            out << setw(width) << areEqual(LCPA[j], LCPB[j]) << ' ';
            //verify LCP
            /*
            gbwt::edge_type curr = {g.toNode(i), j}, above = {g.toNode(i), j-1};
            unsigned lcp = 0;
            while (j != 0 && curr.first == above.first && curr.first != gbwt::ENDMARKER){
                ++lcp;
                curr = g.inverseLF(curr);
                above = g.inverseLF(above);
            }
            if (lcp != LCPA[j])
                out << " brute and FastLCP disagree, brute: " << setw(width) << lcp << " ";
            */
        }
        out << std::flush;
        out << "\n" << setw(labelWidth) << "BWT:";
        for (unsigned j = 0; j < size; ++j){
            out << setw(width) << areEqual(a[j], lfg.toNode(lfg_a[j])) << ' ';
        }
        out << std::flush;
        out << "\n" << setw(labelWidth) << "LF[i]";
        for (unsigned j = 0; j < size; ++j){
            gbwt::edge_type pos(node, j);
            out << setw(width) << areEqual(g.LF(pos), lfg.LF(pos)) << ' ';
        }

        out << std::flush;
        out << "\n" << setw(labelWidth) << "invLF[i]";
        for (unsigned j = 0; j < size; ++j){
            gbwt::edge_type pos(node, j);
            out << setw(width) << areEqual(g.inverseLF(pos), lfg.inverseLF(pos)) << ' ';
        }
        out << "\n" << setw(labelWidth) << "Run:";
        {
            auto lfgit = lfg_a.first.one_begin();
            gbwt::CompressedRecordIterator it(a);
            gbwt::size_type start = 0, newStart, run;
            for (;  lfgit != lfg_a.first.one_end() && !(it.end()); ++it){
                run = lfgit->first;
                ++lfgit;
                newStart = lfgit->second;
                out << '(' << setw(width) << areEqual(it->first, lfg_a.alphabetByRun[run])
                    << ',' << setw(width) << areEqual(it->second, newStart - start) << "),";
                start = newStart;
            }
            out << std::flush;
            out << " Finished runs: " << areEqual((it.end()), lfgit == lfg_a.first.one_end());
            out << std::flush;
        }
        out << '\n' << setw(labelWidth) << "Out:";
        for (unsigned j = 0; j < outgoing.size(); ++j){
            out << '(' << setw(width) << areEqual(outgoing[j].first, g.toNode(lfg_a.alphabet.select_iter(j+1)->second))
                << ',' << setw(width) << areEqual(outgoing[j].second, lfg_a.outgoing[j]) << "),";
        }
        out << std::flush;
        out << '\n';
    }
    out << "\nPrinting paths in GBWT:" << '\n';

    out << std::flush;
    for (unsigned i = 0; i < g.nodeSize(gbwt::ENDMARKER); ++i){
        out << setw(2*width) << "Path" << setw(width) << i << ':';
        gbwt::edge_type curr, init = {gbwt::ENDMARKER, i};
        curr = init;
        gbwt::vector_type Q;
        bool unequal = false;
        do {
            curr = g.LF(curr);
            Q.push_back(curr.first);
        } while (curr.first != init.first);
        Q.pop_back();
        if (Q != g.extract(i)){
            unequal = true;
            out << " Extracted by LF on GBWT and extract() on GBWT NOT EQUAL!\n"
                << setw(labelWidth) << "extract on GBWT:";
            for (const auto& a: g.extract(i))
                out << setw(width) << a;
        }
        if (Q != lfg.extract(i)){
            unequal = true;
            out << " Extracted by LF on GBWT and extract() on LF GBWT NOT EQUAL!\n"
                << setw(labelWidth) << "extract on LF GBWT:";
            for (const auto& a: lfg.extract(i))
                out << setw(width) << a;
        }
        if (unequal)
            out << setw(labelWidth) << "LF on GBWT:";
        for (const auto& a: Q)
            out << setw(width) << a;
        out << '\n';
        out << std::flush;
    }

    out << "\nPrinting full text by general LF: "; 

    auto fullLF = [&l] (gbwt::edge_type & position, gbwt::size_type & sequence){ 
        position = l.rindex->index->LF(position);
        if (position.first == gbwt::ENDMARKER) {
            if (sequence == 0)
                return false;
            --sequence;
            position = {gbwt::ENDMARKER, sequence};
        }
        return true;
    };

    gbwt::size_type sequence = g.sequences()-1;
    gbwt::edge_type position = {gbwt::ENDMARKER, sequence};
    gbwt::vector_type fullTextGBWT;
    fullTextGBWT.push_back(position.first);
    while (fullLF(position, sequence))
        fullTextGBWT.push_back(position.first);
    for (unsigned i = 0; i<fullTextGBWT.size()/2; ++i)
        std::swap(fullTextGBWT[i], fullTextGBWT[fullTextGBWT.size()-1-i]);
    out << fullTextGBWT;

    out << "\nPrinting full text by compr text: " << std::flush;
    gbwt::vector_type fullTextct;
    fullTextct.resize(ct.textLength());
    for (gbwt::size_type i = 0; i < fullTextct.size(); ++i)
        fullTextct[i] = ct.at(i);
    out << fullTextct << '\n' << std::flush;

    out << "Full texts are " << ((fullTextGBWT == fullTextct)? "equal" : "not equal!") << std::endl;
    out << "Full texts are " << ((ct.verifyText())? "equal" : "not equal!") << std::endl;
}
#endif //GBWT_QUERY_IOHELP_H
