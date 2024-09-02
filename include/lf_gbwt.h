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

//Classes in this header are derived from the GBWT and CompressedRecord, with the following license:
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

#ifndef GBWT_QUERY_LF_GBWT_H
#define GBWT_QUERY_LF_GBWT_H

//heavily derived/modified from FastLocate

#include<gbwt/gbwt.h>
#include<gbwt/internal.h>

namespace lf_gbwt{
    typedef std::uint8_t  byte_type;

    class GBWT;

    //all inputs and outputs to CompressedRecord are in the compressed alphabet (toComp)
    struct CompressedRecord
    {
        typedef gbwt::size_type size_type;

        //outgoing[i] stores the value BWT.rank(v,i) where v is the current node
        sdsl::int_vector<0> outgoing;

        //bit vector of length n where n is the length of the current node 
        //bit i is set if a concrete run starts at i
        sdsl::sd_vector<> first;

        //bit vector of length n*\sigma where \sigma is the size of outgoing 
        //bit i is set if a concrete run of value \floor(i/n) starts at i%n
        sdsl::sd_vector<> firstByAlphabet;

        //bit vector of length n where the ith set bit corresponds to the ith set bit in firstByAlphabet
        //the number of bits in [ith set bit, i+1th set bit) is the length of the concrete run corresponding to the
        //ith set bit
        sdsl::sd_vector<> firstByAlphComp;

        //bit vector of length gbwt.effective(), bit i is set if node i is contained in the BWT of this record
        sdsl::sd_vector<> alphabet;

        //int i stores the mapped alphabet value of the BWT run
        sdsl::int_vector<0> alphabetByRun;

        CompressedRecord() = default;
        CompressedRecord(const gbwt::CompressedRecord &, const GBWT *);

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
        void load(std::istream& in);

        size_type size() const { return this->first.size(); }
        bool empty() const { return (this->size() == 0); }
        std::pair<size_type, size_type> runs() const; // (concrete, logical)
        size_type outdegree() const { return this->outgoing.size(); }

        // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
        gbwt::edge_type LF(size_type i) const;

        // Returns `offset` such that `LF(offset) == (to, i)`, or `invalid_offset()`
        // if there is no such offset.
        // This can be used for computing inverse LF in a bidirectional GBWT.
        size_type offsetTo(gbwt::comp_type to, size_type i) const;

        // As above, but also reports the closed offset range ('run') and the identifier
        // ('run_id') of the logical run used for computing LF().
        gbwt::edge_type LF(size_type i, gbwt::range_type& run, size_type& run_id) const;

        // As above, but also sets 'run_end' to the last offset of the current logical run.
        gbwt::edge_type runLF(size_type i, size_type& run_end) const;

        // Returns invalid_offset() if there is no edge to the destination.
        size_type LF(size_type i, gbwt::comp_type to) const;

        // Returns Range::empty_range() if the range is empty or the destination is invalid.
        gbwt::range_type LF(gbwt::range_type range, gbwt::comp_type to) const;

        // the following two functions assume valid i, 0 <= i < this->size()
        size_type compAlphabetAt(size_type i) const { return this->alphabetByRun[this->first.predecessor(i)->first]; }
        // Returns BWT[i] within the record
        gbwt::comp_type operator[](size_type i) const { return this->alphabet.select_iter(this->compAlphabetAt(i)+1)->second; }

        bool hasEdge(gbwt::comp_type to) const { return this->edgeTo(to) != gbwt::invalid_offset(); }

        // Maps successor nodes to outranks.
        gbwt::rank_type edgeTo(gbwt::comp_type to) const { 
            auto iter = this->alphabet.predecessor(to);
            return (iter != this->alphabet.one_end() && iter->second == to)? iter->first : gbwt::invalid_offset(); 
        };

        // These assume that 'outrank' is a valid outgoing edge.
        gbwt::comp_type successor(gbwt::rank_type outrank) const { return this->alphabet.select_iter(outrank+1)->second; }
        //gbwt::node_type successor(gbwt::rank_type outrank) const { return this->alphabet.successor(outrank)->second; }
        
        //what is this? when is it used?
        //size_type offset(gbwt::rank_type outrank) const { return this->outgoing[outrank]; }

        //returns the run id by logical runs of an offset
        size_type logicalRunId(size_type i) const {
            if (i >= this->size())
                return gbwt::invalid_offset();
            size_type concRunId = this->first.predecessor(i)->first;
            if (!this->hasEdge(gbwt::ENDMARKER))
                return concRunId;
            size_type num0RunsBefore = this->firstByAlphabet.successor(i)->first;
            size_type num0Before = this->firstByAlphComp.select_iter(num0RunsBefore+1)->second;
            //error if in the middle of an endmarker run? counts endmarkers past current
            //if ((*this)[i] == gbwt::ENDMARKER) {
            //  num0Before -= (++this->first.predecessor(i))->second - i;
            return concRunId - num0RunsBefore + num0Before;
        }
    };

    //The GBWT is split into two subsets, one of small records and one of large records
    //small records are those with outdegree <  maxOutDegree
    struct SmallRecordArray{
        typedef gbwt::size_type size_type;

        size_type maxOutdegree = 0;
        size_type records = 0;
        size_type effective = 0;

        sdsl::sd_vector<> outDegreePrefixSum;

        sdsl::sd_vector<> emptyRecords;

        //bit vector, if the c-th set bit (0-indexed) is at position i in prefixSum,
        //then, i is c + the sum of the lengths of all nodes < c 
        //has \sigma set bits
        //note by GBWT assumption, \sigma < 2^32 and sum of their lengths <= (2^32 - 1)^2
        //therefore the max length of this array is 2^64  - 2^33 + 1 + 2^32 - 1 = 2^64 - 2^32
        sdsl::sd_vector<> prefixSum;

        //outgoing[i*maxOutdegree + j] stores the value BWT.rank(v,j) where v is the i-th node in the small record array
        //above wrong, use prefix sum
        sdsl::int_vector<0> outgoing;

        //stores start positions of runs of each node concatenated, starting at prefixSum
        sdsl::sd_vector<> first;

        //bit vector of length prefixSum.length() * maxOutdegree
        //call prefixSum(i) the position of the i-th set bit in prefix sum
        //call the i-th window (0-indexed) [maxOutdegree*prefixSum(i), maxOutdegree*prefixSum(i) + (prefixSum(i+1) - prefixSum(i+) - 1)*maxOutdegree)
        //the has maxOutdegree sections, each of length prefixSum(i+1) - prefixSum(i) - 1 (the length of the node)
        //the i-th window is equivalent to firstByAlphabet if this node was a CompressedRecord
        sdsl::sd_vector<>firstByAlphabet;


        //bit vector of length prefixSum.length()
        //call the i-th window (0-indexed) [prefixSum(i), prefixSum(i+1))
        //this windodw is equivalent to firstByAlphComp if node i was a CompressedRecord
        sdsl::sd_vector<> firstByAlphComp; 

        //bit vector of length \sigma^2
        //\sigma is the size the number of nodes in the GBWT ( = gbwt.effective())
        //note by GBWT assumption, \sigma <= 2^32
        //the i-th bit is set if the i/\sigma-th node has an edge to i%\sigma
        //
        //emptynodes windows of length gbwt.effective()
        sdsl::sd_vector<> alphabet;

        //stores the mapped alphabet value of the BWT run
        //concatenated for all smallRecords
        sdsl::int_vector<0> alphabetByRun;

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
        void load(std::istream& in);

        //these assume node is a valid index in smallRecordArray
        size_type size(const size_type node) const;
        bool empty(const size_type node) const;
        std::pair<bool, size_type> emptyAndNonEmptyIndex(const size_type node) const;
        std::pair<size_type, size_type> runs(const size_type node) const;
        size_type outdegree(const size_type node) const;
        size_type nonEmptyRecords() const { return records - emptyRecords.ones(); }

        gbwt::edge_type LF(const size_type node, const size_type pos) const;

        size_type offsetTo(const size_type node, const gbwt::comp_type to, size_type i) const;

        gbwt::edge_type LF(const size_type node, const size_type i, gbwt::range_type& run, size_type& run_id) const;

        gbwt::edge_type runLF(const size_type node, const size_type i, size_type& run_end) const;

        size_type LF(const size_type node, const size_type i, const gbwt::comp_type to) const;

        gbwt::range_type LF(const size_type node, const gbwt::range_type range, const gbwt::comp_type to) const;

        size_type compAlphabetAt(const size_type node, const size_type i) const;

        gbwt::comp_type bwtAt(const size_type node, const size_type i) const;

        bool hasEdge(const size_type node, const gbwt::comp_type to) const;

        gbwt::rank_type edgeTo(const gbwt::size_type node, const gbwt::comp_type to) const;

        gbwt::comp_type successor(size_type node, gbwt::rank_type outrank) const;

        size_type offset(size_type node, gbwt::rank_type outrank) const;

        size_type logicalRunId(const size_type node, const size_type i) const;

    };

    class GBWT
    {
        public:
            typedef gbwt::size_type size_type;

            //------------------------------------------------------------------------------

            GBWT() = default;
            explicit GBWT(const gbwt::GBWT& source);

            size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
            void load(std::istream& in);
            
            const static std::string EXTENSION; // .lfgbwt

            //------------------------------------------------------------------------------

            /*
               Low-level interface: Statistics.
            */

            size_type size() const { return this->header.size; }
            bool empty() const { return (this->size() == 0); }
            size_type sequences() const { return this->header.sequences; }
            size_type sigma() const { return this->header.alphabet_size; }
            size_type effective() const { return this->header.alphabet_size - this->header.offset; }

            std::pair<size_type, size_type> runs() const;

            bool bidirectional() const { return this->header.get(gbwt::GBWTHeader::FLAG_BIDIRECTIONAL); }

            //------------------------------------------------------------------------------

            /*
               Metadata interface.
            */

            bool hasMetadata() const { return this->header.get(gbwt::GBWTHeader::FLAG_METADATA); }
            void addMetadata() { this->header.set(gbwt::GBWTHeader::FLAG_METADATA); }
            void clearMetadata() { this->metadata.clear(); this->header.unset(gbwt::GBWTHeader::FLAG_METADATA); };

            //------------------------------------------------------------------------------

            /*
               High-level interface. The queries check that the parameters are valid. Iterators
               must be InputIterators. On error or failed search, the return values will be the
               following:

               find     empty search state
               prefix   empty search state
               extend   empty search state
               locate   invalid_sequence() or empty vector
               extract  empty vector
            */

            gbwt::vector_type extract(size_type sequence) const;

            bool contains(gbwt::node_type node) const
            {
                return ((node < this->sigma() && node > this->header.offset) || node == gbwt::ENDMARKER);
            }

            gbwt::node_type firstNode() const { return this->header.offset + 1; }
            gbwt::comp_type toComp(gbwt::node_type node) const { return (node == 0 ? node : node - this->header.offset); }
            gbwt::node_type toNode(gbwt::comp_type comp) const { return (comp == 0 ? comp : comp + this->header.offset); }

            size_type nodeSize(gbwt::node_type node) const { auto ind = this->isSmallAndIndex(this->toComp(node)); return (ind.first)? smallRecords.size(ind.second) : this->largeRecords[ind.second].size(); }
            bool empty(gbwt::node_type node) const { return this->nodeSize(node) == 0; }

            //------------------------------------------------------------------------------

            /*
               Low-level interface: Navigation and searching. The interface assumes that node
               identifiers are valid. This can be checked with contains().
            */

            //returns pair first value is true if node is small, second value is index in smallRecords if small, otherwise index in largeRecords
            std::pair<bool, gbwt::size_type> isSmallAndIndex(gbwt::comp_type comp) const {
                assert(comp < this->effective());
                auto nextSmall = isSmall.successor(comp);
                return (nextSmall->second == comp)? std::pair<bool,gbwt::size_type>{true, nextSmall->first} : std::pair<bool,gbwt::size_type>{false, comp - nextSmall->first};
            }

            // On error: invalid_edge().
            gbwt::edge_type LF(gbwt::edge_type position) const
            {
                auto ind = this->isSmallAndIndex(this->toComp(position.first));
                gbwt::edge_type ans = (ind.first)? this->smallRecords.LF(ind.second, position.second) : this->largeRecords[ind.second].LF(position.second);
                ans.first = this->toNode(ans.first);
                return ans;
            }
            
            gbwt::edge_type LF(gbwt::edge_type position, gbwt::range_type& run, gbwt::size_type& run_id) const {
                auto ind = this->isSmallAndIndex(this->toComp(position.first));
                gbwt::edge_type next = (ind.first)? this->smallRecords.LF(ind.second, position.second, run, run_id) : this->largeRecords[ind.second].LF(position.second, run, run_id);
                next.first = this->toNode(next.first);
                return next;
            }

            // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
            // On error: invalid_edge().
            gbwt::edge_type inverseLF(gbwt::node_type from, size_type i) const;

            // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
            // On error: invalid_edge().
            gbwt::edge_type inverseLF(gbwt::edge_type position) const
            {
                return this->inverseLF(position.first, position.second);
            }
            
            gbwt::node_type predecessorAt(gbwt::node_type, size_type i) const;
            //------------------------------------------------------------------------------

            /*
               Low-level interface: Sequences. The interface assumes that node identifiers are
               valid. This can be checked with contains().
            */

            //------------------------------------------------------------------------------

            gbwt::GBWTHeader  header;
            gbwt::Tags        tags;
            sdsl::sd_vector<> isSmall;
            SmallRecordArray smallRecords;
            std::vector<CompressedRecord> largeRecords;
            gbwt::Metadata    metadata;

            // Decompress and cache the endmarker, because decompressing it is expensive.
            //CompressedRecord *endmarker_record;

            //------------------------------------------------------------------------------

            /*
               Internal interface. Do not use.
            */

            /*const CompressedRecord& record(gbwt::node_type node) const { 
                if (!this->contains(node))
                    throw std::invalid_argument("node given to lf_gbwt.record invalid!");
                return this->bwt[this->toComp(node)]; 
            }*/

            bool verify(const gbwt::GBWT& g) {
                bool equal = (this->size() == g.size()) && (this->sequences() == g.sequences()) && (this->sigma() == g.sigma()) && (this->effective() == g.effective()) && (this->runs() == g.runs()) && (this->bidirectional() == g.bidirectional());
                equal = equal && this->verifyBWT(g);
                equal = equal && this->verifyLF(g);
                equal = equal && (!this->bidirectional() || this->verifyInverseLF(g));
                return equal;
            }

            //assumes number of nodes and node ids are equivalent
            bool verifyBWT(const gbwt::GBWT& g) {
                bool equal = true;
                #pragma omp parallel for schedule(dynamic, 1)
                for (gbwt::comp_type i = 0; i < this->effective(); ++i) {
                    auto ind = this->isSmallAndIndex(i);
                    bool nodeEqual = this->toNode(i) == g.toNode(i);
                    gbwt::CompressedRecord rec = g.record(g.toNode(i));
                    nodeEqual = nodeEqual && this->nodeSize(this->toNode(i)) == rec.size();
                    if (ind.first) {
                        for (gbwt::size_type j = 0; nodeEqual && j < this->nodeSize(this->toNode(i)); ++j)
                            nodeEqual = nodeEqual && rec[j] == this->toNode(this->smallRecords.bwtAt(ind.second, j));
                    }
                    else {
                        const CompressedRecord &lfrec = this->largeRecords[ind.second];
                        for (gbwt::size_type j = 0; nodeEqual && j < rec.size(); ++j)
                            nodeEqual = nodeEqual && rec[j] == this->toNode(lfrec[j]);
                    }
                    #pragma omp critical 
                    {
                       equal = equal && nodeEqual;
                    }
                }
                return equal;
            }

            bool verifyLF(const gbwt::GBWT& g) {
                bool equal = true;
                #pragma omp parallel for schedule(dynamic, 1)
                for (gbwt::comp_type i = 0; i < this->effective(); ++i) {
                    gbwt::node_type node = this->toNode(i);
                    bool nodeEqual = node == g.toNode(i);
                    gbwt::size_type nodeSize = g.nodeSize(node);
                    nodeEqual = nodeEqual && nodeSize == this->nodeSize(node);
                    for (gbwt::size_type j = 0; nodeEqual && j < nodeSize; ++j)
                        nodeEqual = nodeEqual && this->LF({node, j}) == g.LF({node, j});
                    #pragma omp critical 
                    {
                       equal = equal && nodeEqual;
                    }
                }
                return equal;
            }

            bool verifyInverseLF(const gbwt::GBWT& g) {
                bool equal = this->bidirectional() && g.bidirectional();
                #pragma omp parallel for schedule(dynamic, 1)
                for (gbwt::comp_type i = 0; i < this->effective(); ++i) {
                    gbwt::node_type node = this->toNode(i);
                    bool nodeEqual = node == g.toNode(i);
                    gbwt::size_type nodeSize = g.nodeSize(node);
                    nodeEqual = nodeEqual && nodeSize == this->nodeSize(node);
                    for (gbwt::size_type j = 0; nodeEqual && j < nodeSize; ++j)
                        nodeEqual = nodeEqual && this->inverseLF({node, j}) == g.inverseLF({node, j});
                    #pragma omp critical 
                    {
                       equal = equal && nodeEqual;
                    }
                }
                return equal;
            }
    }; // class GBWT

    const std::string GBWT::EXTENSION = ".lfgbwt"; // .lfgbwt

    //CompressedRecord member functions

    CompressedRecord::CompressedRecord(const gbwt::CompressedRecord & rec, const GBWT *source){
        size_type n = rec.size(), sigma = rec.outgoing.size(), runs = rec.runs().first;
        if (n == 0)
            return;

        //outgoing
        this->outgoing.resize(sigma);
        for (unsigned i = 0; i < sigma; ++i)
            this->outgoing[i] = rec.outgoing[i].second;
        sdsl::util::bit_compress(this->outgoing);

        //alphabet
        sdsl::sd_vector_builder alphabetBuilder(source->effective(), sigma);
        for (unsigned i = 0; i < sigma; ++i)
            alphabetBuilder.set(source->toComp(rec.outgoing[i].first));
        this->alphabet = sdsl::sd_vector<>(alphabetBuilder);

        //first, firstByAlphabet, firstByAlphComp, and alphabetByRun;
        sdsl::sd_vector_builder 
            firstBuilder(n, runs),
            firstByAlphabetBuilder(n*sigma, runs),
            firstByAlphCompBuilder(n, runs);
        this->alphabetByRun.resize(runs);

        //start location of run in firstByAlphabetAssist and length of run
        std::vector<std::vector<std::pair<size_type,size_type>>> firstByAlphabetAssist;
        firstByAlphabetAssist.resize(sigma);

        gbwt::CompressedRecordFullIterator iter(rec);
        size_type start, runInd = 0;
        while (!iter.end()){
            start = iter.offset() - iter.run.second;
            firstBuilder.set(start);
            firstByAlphabetAssist[iter.run.first].emplace_back(iter.run.first*n + start, iter.run.second);
            this->alphabetByRun[runInd] = iter.run.first;
            ++iter;
            ++runInd;
        }
        this->first = sdsl::sd_vector<>(firstBuilder);
        sdsl::util::bit_compress(this->alphabetByRun);

        size_type currInd = 0;
        for (const auto & alphArr : firstByAlphabetAssist){
            for (const auto & a : alphArr){
                firstByAlphabetBuilder.set(a.first);
                firstByAlphCompBuilder.set(currInd);
                currInd += a.second;
            }
        }
        this->firstByAlphabet = sdsl::sd_vector<>(firstByAlphabetBuilder);
        this->firstByAlphComp = sdsl::sd_vector<>(firstByAlphCompBuilder);
    }

    CompressedRecord::size_type CompressedRecord::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += sdsl::serialize(this->outgoing, out, child, "outgoing");
        written_bytes += sdsl::serialize(this->first, out, child, "first");
        written_bytes += sdsl::serialize(this->firstByAlphabet, out, child, "firstByAlphabet");
        written_bytes += sdsl::serialize(this->firstByAlphComp, out, child, "firstByAlphComp");
        written_bytes += sdsl::serialize(this->alphabet, out, child, "alphabet");
        written_bytes += sdsl::serialize(this->alphabetByRun, out, child, "alphabetByRun");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void CompressedRecord::load(std::istream& in) {
        sdsl::load(this->outgoing, in);
        sdsl::load(this->first, in);
        sdsl::load(this->firstByAlphabet, in);
        sdsl::load(this->firstByAlphComp, in);
        sdsl::load(this->alphabet, in);
        sdsl::load(this->alphabetByRun, in);
    }

    std::pair<CompressedRecord::size_type, CompressedRecord::size_type> CompressedRecord::runs() const {
        size_type totRuns = this->first.ones();
        if (!this->hasEdge(gbwt::ENDMARKER)) 
            return {totRuns, totRuns};
        if (this->outgoing.size() == 1)
            return {totRuns, this->size()};
        size_type numENDMARKERruns = this->firstByAlphabet.predecessor(this->size()-1)->first + 1;
        size_type numENDMARKERS = this->firstByAlphComp.select_iter(numENDMARKERruns+1)->second;
        return {totRuns, totRuns - numENDMARKERruns + numENDMARKERS};
    }; // (concrete, logical)

    // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
    gbwt::edge_type CompressedRecord::LF(size_type i) const {
        if (i >= this->size())
            return gbwt::invalid_edge();
        gbwt::comp_type next = (*this)[i];
        return {next, this->LF(i, next)};
    }

    //SmallRecordArray member functions

    SmallRecordArray::size_type SmallRecordArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += sdsl::serialize(this->maxOutdegree, out, child, "maxOutdegree");
        written_bytes += sdsl::serialize(this->records, out, child, "records");
        written_bytes += sdsl::serialize(this->outDegreePrefixSum, out, child, "outDegreePrefixSum");
        written_bytes += sdsl::serialize(this->emptyRecords, out, child, "emptyRecords");
        written_bytes += sdsl::serialize(this->prefixSum, out, child, "prefixSum");
        written_bytes += sdsl::serialize(this->outgoing, out, child, "outgoing");
        written_bytes += sdsl::serialize(this->first, out, child, "first");
        written_bytes += sdsl::serialize(this->firstByAlphabet, out, child, "firstByAlphabet");
        written_bytes += sdsl::serialize(this->firstByAlphComp, out, child, "firstByAlphComp");
        written_bytes += sdsl::serialize(this->alphabet, out, child, "alphabet");
        written_bytes += sdsl::serialize(this->alphabetByRun, out, child, "alphabetByRun");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void SmallRecordArray::load(std::istream& in) {
        sdsl::load(this->maxOutdegree, in);
        sdsl::load(this->records, in);
        sdsl::load(this->outDegreePrefixSum, in);
        sdsl::load(this->emptyRecords, in);
        sdsl::load(this->prefixSum, in);
        sdsl::load(this->outgoing, in);
        sdsl::load(this->first, in);
        sdsl::load(this->firstByAlphabet, in);
        sdsl::load(this->firstByAlphComp, in);
        sdsl::load(this->alphabet, in);
        sdsl::load(this->alphabetByRun, in);
    }

    SmallRecordArray::size_type SmallRecordArray::size(const SmallRecordArray::size_type node) const {
        assert(node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return 0;
        auto it = prefixSum.select_iter(p.second + 1);
        size_type beg = (it++)->second;
        return it->second - beg;
    }

    bool SmallRecordArray::empty(const size_type node) const {
        assert(node < records);
        return emptyRecords[node];
    }

    std::pair<bool, SmallRecordArray::size_type> SmallRecordArray::emptyAndNonEmptyIndex(const SmallRecordArray::size_type node) const {
        assert(node < records);
        auto it = emptyRecords.successor(node);
        return {it->second == node, node - it->first};
    }

    std::pair<SmallRecordArray::size_type, SmallRecordArray::size_type> SmallRecordArray::runs(const SmallRecordArray::size_type node) const {
        assert (node < records);
        auto pos = emptyAndNonEmptyIndex(node);
        if (pos.first) return {0, 0};
        auto end = prefixSum.select_iter(pos.second + 1);
        auto beg = end++; 
        auto runBeg = first.successor(beg->second);
        auto runEnd = first.successor(end->second);
        assert(runBeg->second == beg->second);
        assert(runEnd->second == end->second);
        size_type concRuns = runEnd->first - runBeg->first;
        //check if this has an edge to the endmarker
        if (!alphabet[pos.second*effective])
            return {concRuns, concRuns};
        auto begAlphRun = firstByAlphabet.successor(maxOutdegree*beg->second);
        auto endAlphRun = firstByAlphabet.successor(maxOutdegree*beg->second + (end->second - beg->second));
        size_type numENDMARKERruns = endAlphRun->first - begAlphRun->first;
        assert(numENDMARKERruns > 0);
        size_type numENDMARKERS = firstByAlphComp.select_iter(endAlphRun->first)->second - firstByAlphComp.select_iter(begAlphRun->first)->second;
        assert(numENDMARKERS > 0);
        return {concRuns, concRuns - numENDMARKERruns + numENDMARKERS};
    }

    SmallRecordArray::size_type SmallRecordArray::outdegree(const SmallRecordArray::size_type node) const {
        assert(node < records);
        auto pos = emptyAndNonEmptyIndex(node);
        if (pos.first) return 0;
        size_type n = nonEmptyRecords(); 
        size_type ans = alphabet.successor(n*pos.second + n)->first - alphabet.successor(effective*pos.second)->first;
        assert(ans != 0);
        return ans;
    }

    gbwt::edge_type SmallRecordArray::LF(const SmallRecordArray::size_type node, const SmallRecordArray::size_type pos) const {
        assert (node < records);
        if (pos >= size(node))
            return gbwt::invalid_edge();
        gbwt::comp_type next = bwtAt(node, pos);
        return {next, LF(node, next, pos)};
    }
    
    SmallRecordArray::size_type SmallRecordArray::offsetTo(const SmallRecordArray::size_type node, const gbwt::comp_type to, SmallRecordArray::size_type i) const {
        assert(node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return gbwt::invalid_offset();
        size_type outrank = edgeTo(node, to);
        if (outrank == gbwt::invalid_offset()) return gbwt::invalid_offset();
        size_type n = nonEmptyRecords();
        size_type outgoingPrefixSum = alphabet.successor(p.second * effective)->first;
        i -= outgoing[outgoingPrefixSum + outrank];
        auto prefixEnd = prefixSum.select_iter(p.second + 1);
        auto prefixBeg = prefixEnd++;
        size_type lengthPrefixSum = prefixBeg->second;
        size_type nodeLength = prefixEnd->second - prefixBeg->second;

        auto firstToRun = firstByAlphabet.successor((lengthPrefixSum * maxOutdegree) + nodeLength*outrank);
        auto firstToCompRun = firstByAlphComp.select_iter(firstToRun->first+1);

        auto runComp = firstByAlphComp.predecessor(firstToCompRun->second + i);
        i -= (runComp->second - firstToCompRun->second);
        auto run = firstByAlphabet.select_iter(runComp->first+1);
        if ((run->second - (lengthPrefixSum * maxOutdegree))/nodeLength != outrank)
            return gbwt::invalid_offset();
        return (run->second - (lengthPrefixSum*maxOutdegree) - (outrank*nodeLength)) + i;
    }

    gbwt::edge_type SmallRecordArray::LF(const SmallRecordArray::size_type node, const SmallRecordArray::size_type i, gbwt::range_type& run, SmallRecordArray::size_type& run_id) const {
        assert(node < records);
        size_type nodeLength = size(node);
        if (i >= nodeLength){
            run.first = run.second = gbwt::invalid_offset();
            run_id = gbwt::invalid_offset();
            return gbwt::invalid_edge();
        }
        size_type nonEmptyIndex = emptyAndNonEmptyIndex(node).second;
        run_id = logicalRunId(node, i);
        auto prefixBeg = prefixSum.select_iter(nonEmptyIndex + 1);
        auto it = first.predecessor(prefixBeg->second + i);
        run.first = it->second - prefixBeg->second;
        ++it;
        run.second = it->second - prefixBeg->second - 1;

        return LF(node, i);
    }

    SmallRecordArray::size_type SmallRecordArray::LF(const SmallRecordArray::size_type node, const SmallRecordArray::size_type i, const gbwt::comp_type to) const {
        assert(node < records);
        size_type nodeLength = size(node);
        bool isEmpty;
        size_type nonEmptyIndex;
        std::tie(isEmpty, nonEmptyIndex) = emptyAndNonEmptyIndex(node);
        //why > instead of >=?
        if (isEmpty || i > nodeLength || !hasEdge(node, to))
            return gbwt::invalid_offset();
        size_type outrank = edgeTo(node, to);
        if (outrank == gbwt::invalid_offset()) return gbwt::invalid_offset();
        auto prefixBeg = prefixSum.select_iter(nonEmptyIndex + 1);
        size_type lengthPrefixSum = prefixBeg->second;
        auto nextRun = firstByAlphabet.successor((lengthPrefixSum*maxOutdegree) + (outrank*nodeLength) + i);
        auto firstOutrankRun = firstByAlphabet.successor((lengthPrefixSum*maxOutdegree) + (outrank*nodeLength));
        auto nextRunComp = firstByAlphComp.select_iter(nextRun->first+1);
        auto firstOutrankComp = firstByAlphComp.select_iter(firstOutrankRun->first+1);
        size_type numOutrank = nextRunComp->second - firstOutrankComp->second;
        if (bwtAt(node, i) == to) {
            size_type afteriInRun = this->first.successor(lengthPrefixSum + i)->second - lengthPrefixSum - i;
            numOutrank -= afteriInRun;
        }
        size_type outgoingPrefixSum = alphabet.successor(nonEmptyIndex*effective)->first;

        return outgoing[outgoingPrefixSum + outrank] + numOutrank;
    }

    SmallRecordArray::size_type SmallRecordArray::compAlphabetAt(const SmallRecordArray::size_type node, const SmallRecordArray::size_type i) const {
        assert (node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return gbwt::invalid_offset();
        if (i >= size(node)) return gbwt::invalid_offset();
        size_type lengthPrefixSum = prefixSum.select_iter(p.second + 1)->second;
        size_type runInd = first.predecessor(lengthPrefixSum + i)->first;
        return alphabetByRun[runInd];
    }

    SmallRecordArray::size_type SmallRecordArray::bwtAt(const SmallRecordArray::size_type node, const SmallRecordArray::size_type i) const {
        assert (node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return gbwt::invalid_offset();
        if (i >= size(node)) return gbwt::invalid_offset();
        size_type outrank = compAlphabetAt(node, i);
        auto it = alphabet.select_iter(outrank + alphabet.successor(p.second*effective)->first + 1);
        return it->second - p.second*nonEmptyRecords();
    }

    bool SmallRecordArray::hasEdge(const SmallRecordArray::size_type node, const gbwt::comp_type to) const {
        return edgeTo(node, to) != gbwt::invalid_offset();
    }

    gbwt::rank_type SmallRecordArray::edgeTo(const gbwt::size_type node, const gbwt::comp_type to) const {
        assert(node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return gbwt::invalid_offset();
        auto it = alphabet.successor(effective*p.second + to);
        return (it->second == effective*p.second + to)? it->first - alphabet.successor(effective*p.second)->first : gbwt::invalid_offset();
    }

    gbwt::comp_type SmallRecordArray::successor(const SmallRecordArray::size_type node, const gbwt::rank_type outrank) const {
        assert(node < records);
        auto p = emptyAndNonEmptyIndex(node);
        if (p.first) return gbwt::invalid_offset();
        auto it = alphabet.select_iter(alphabet.successor(effective*p.second)->first + outrank + 1);
        if (it->second >= effective*(p.second+1))
            return gbwt::invalid_offset();
        return it->second - (effective*p.second);
    }

    SmallRecordArray::size_type SmallRecordArray::logicalRunId(const SmallRecordArray::size_type node, const SmallRecordArray::size_type i) const {
        assert(node < records);
        size_type nodeLength = size(node);
        if (i >= nodeLength)
            return gbwt::invalid_offset();
        size_type nonEmptyIndex = emptyAndNonEmptyIndex(node).second;
        auto prefixBeg = prefixSum.select_iter(nonEmptyIndex + 1);
        auto it = first.predecessor(prefixBeg->second + i);
        assert(first[prefixBeg->second] && first.predecessor(prefixBeg->second)->second == prefixBeg->second);
        size_type concRunId = it->first - first.predecessor(prefixBeg->second)->first;
        if (!hasEdge(node, gbwt::ENDMARKER))
            return concRunId;
        size_type num0RunsBefore = firstByAlphabet.successor(nonEmptyIndex*effective+i)->first;
        size_type num0Before = firstByAlphComp.select_iter(num0RunsBefore + 1)->second;
        auto alphPref = firstByAlphabet.successor(nonEmptyIndex*effective);
        num0RunsBefore -= alphPref->first;
        num0Before -= this->firstByAlphComp.select_iter(alphPref->first + 1)->second;
        if (bwtAt(node, i) == gbwt::ENDMARKER) {
            //remove 0s after i from num0before
            ++it;
            num0Before -= it->second - i - prefixBeg->second;
        }
        return concRunId - num0RunsBefore + num0Before;
    }

    //GBWT member functions

    gbwt::node_type GBWT::predecessorAt(gbwt::node_type revFrom, size_type i) const {
        auto revNode = [this](gbwt::comp_type x) { 
            if (x == gbwt::ENDMARKER) 
                return x;
            return gbwt::Node::reverse(this->toNode(x));
        };
        auto ind = this->isSmallAndIndex(this->toComp(revFrom));
        if (i >= this->nodeSize(revFrom))
            return gbwt::invalid_node();

        if (ind.first) {
            size_type prefixLength = this->smallRecords.prefixSum.select_iter(ind.second + 1)->second;
            size_type predoutrank = (this->smallRecords.firstByAlphabet.select_iter(
                    this->smallRecords.firstByAlphComp.predecessor(prefixLength + i)->first + 1
                    )->second - (prefixLength*this->smallRecords.maxOutdegree))/this->smallRecords.size(ind.second);

            //check if before predoutrank is its reverse
            auto alphIter = this->smallRecords.alphabet.select_iter(
                    this->smallRecords.alphabet.successor(ind.second * this->effective())->first + 1 + predoutrank
                    );
            if (predoutrank > 0) {
                auto prevAlphIter = alphIter;
                --prevAlphIter;
                if (this->toNode(prevAlphIter->second - ind.second * this->effective()) == revNode(alphIter->second - ind.second * this->effective())) {
                    size_type beforePrevAlph = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank-1)*this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    size_type prevAlphSize   = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank)  *this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    size_type AlphSize       = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank+1)*this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    AlphSize -= prevAlphSize;
                    prevAlphSize -= beforePrevAlph;
                    beforePrevAlph -= prefixLength;
                    if (i < beforePrevAlph + AlphSize)
                        return revNode(alphIter->second - ind.second*this->effective());
                    return this->toNode(alphIter->second - ind.second*this->effective());
                }
            }
            
            //check if after predoutrank is its reverse
            auto afterAlphIter = alphIter;
            ++afterAlphIter;
            if (afterAlphIter->second/this->effective() == ind.second) {
                if (this->toNode(afterAlphIter->second - ind.second * this->effective()) == revNode(alphIter->second - ind.second * this->effective())) {
                    size_type beforeAlph    = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank)  *this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    size_type AlphSize      = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank+1)*this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    size_type AfterAlphSize = this->smallRecords.firstByAlphComp.select_iter(1+
                            this->smallRecords.firstByAlphabet.successor((predoutrank+2)*this->smallRecords.size(ind.second) + prefixLength*this->smallRecords.maxOutdegree)->first
                            )->second;
                    AfterAlphSize -= AlphSize;
                    AlphSize -= beforeAlph;
                    beforeAlph -= prefixLength;
                    if (i < beforeAlph + AfterAlphSize)
                        return revNode(afterAlphIter->second - ind.second*this->effective());
                    return this->toNode(afterAlphIter->second - ind.second*this->effective());

                }
            }

            return revNode(alphIter->second - ind.second*this->effective());
        }
        else {
            const CompressedRecord & rev = this->largeRecords[ind.second];
            size_type predoutrank = rev.firstByAlphabet.select_iter(
                    rev.firstByAlphComp.predecessor(i)->first+1
                    )->second/rev.size();

            //check if before predoutrank is its reverse
            auto alphIter = rev.alphabet.select_iter(predoutrank+1);
            if (predoutrank > 0){
                auto prevAlphIter = --alphIter;
                ++alphIter;
                if (this->toNode(prevAlphIter->second) == revNode(alphIter->second)){
                    size_type beforePrevAlph = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank-1)*rev.size())->first
                            )->second;
                    size_type prevAlphSize = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank)*rev.size())->first
                            )->second;
                    size_type AlphSize = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank+1)*rev.size())->first
                            )->second;
                    AlphSize -= prevAlphSize;
                    prevAlphSize -= beforePrevAlph;

                    if (i < beforePrevAlph + AlphSize)
                        return revNode(alphIter->second);
                    return this->toNode(alphIter->second);
                }
            }

            //check if after predoutrank is its reverse 
            if (predoutrank < rev.outgoing.size()-1){
                auto afterAlphIter = ++alphIter;
                --alphIter;
                if (this->toNode(afterAlphIter->second) == revNode(alphIter->second)){
                    size_type beforeAlph = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank)*rev.size())->first
                            )->second;
                    size_type AlphSize = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank+1)*rev.size())->first
                            )->second;
                    size_type AfterAlphSize = rev.firstByAlphComp.select_iter(1+
                            rev.firstByAlphabet.successor((predoutrank+2)*rev.size())->first
                            )->second;
                    AfterAlphSize -= AlphSize;
                    AlphSize -= beforeAlph;

                    if (i < beforeAlph + AfterAlphSize)
                        return revNode(afterAlphIter->second);
                    return this->toNode(afterAlphIter->second);
                }
            }

            return revNode(alphIter->second);
        }
    }
        

    // As above, but also reports the closed offset range ('run') and the identifier
    // ('run_id') of the logical run used for computing LF().
    gbwt::edge_type CompressedRecord::LF(size_type i, gbwt::range_type& run, size_type& run_id) const {
        if (i >= this->size()){
            run.first = run.second = gbwt::invalid_offset();
            run_id = gbwt::invalid_offset();
            return gbwt::invalid_edge();
        }
        run_id = this->logicalRunId(i);
        auto it = this->first.predecessor(i);
        run.first = it->second;
        ++it;
        run.second = it->second-1;

        return this->LF(i);
    }

    // Returns invalid_offset() if there is no edge to the destination.
    CompressedRecord::size_type CompressedRecord::LF(size_type i, gbwt::comp_type to) const {
        if (i > this->size() || !this->hasEdge(to))
            return gbwt::invalid_offset();
        size_type outrank = this->edgeTo(to);
        auto nextRun = this->firstByAlphabet.successor(outrank*this->size() + i);
        auto firstOutrankRun = this->firstByAlphabet.successor(outrank*this->size());
        auto nextRunComp = this->firstByAlphComp.select_iter(nextRun->first+1);
        auto firstOutrankRunComp = this->firstByAlphComp.select_iter(firstOutrankRun->first+1);
        size_type numOutrank = nextRunComp->second 
            - firstOutrankRunComp->second;
        if ((*this)[i] == to){
            size_type afteriInRun = this->first.successor(i)->second - i;
            numOutrank -= afteriInRun;
        }
        
        return this->outgoing[outrank] + numOutrank;
    }


    GBWT::GBWT(const gbwt::GBWT& source) {
        double start = gbwt::readTimer();
        this->header = source.header;
        this->tags = source.tags;
        this->metadata = source.metadata;

        double bwtStart = gbwt::readTimer();
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            std::cerr << "lf_GBWT::GBWT::GBWT(): Constructing " << source.effective() << " nodes of total length " << this->size() << std::endl;
        }
        std::vector<CompressedRecord> bwt;
        bwt.resize(source.effective());
        source.bwt.forEach([this, &bwt] (size_type a, const gbwt::CompressedRecord& rec) {bwt[a] = CompressedRecord(rec, this);});
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            double seconds = gbwt::readTimer() - bwtStart;
            std::cerr << "lf_GBWT::GBWT::GBWT(): Constructed " << source.effective() << " nodes of total length " << this->size() << " in " << seconds << " seconds" << std::endl;
        }

        double smallRecordsStart = gbwt::readTimer();
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            std::cerr << "lf_GBWT::GBWT::GBWT(): Computing maxOutdegree for SmallRecordArray " << std::endl;
        }
        size_type maxOutdegreeFound = 0;
        //maps from outdegree a tuple (a,b,c) where 
        //a is the number of bytes of compressed records with that outdegree, 
        //b is the number of compressed records with that outdegree, and 
        //c is the total length of all nodes of this size
        std::unordered_map<size_type, std::tuple<size_type, size_type, size_type>> outdegreeCounts;
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_type i = 0; i < bwt.size(); ++i) {
            size_type outdegree = bwt[i].outdegree();
            size_type bytes = sdsl::size_in_bytes(bwt[i]);
            size_type length = bwt[i].size();

            #pragma omp critical 
            {
                std::tuple<size_type, size_type, size_type> &val = outdegreeCounts[outdegree];
                maxOutdegreeFound = std::max(maxOutdegreeFound, outdegree);
                //if C++17 supported, replace with structured binding declaration
                //at the moment minimum is C++14
                std::get<0>(val) += bytes;
                std::get<1>(val)++;
                std::get<2>(val) += length;
            }
        }
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            std::cerr << "lf_GBWT::GBWT::GBWT(): Finished getting sizes, counts, and lengths by outdegree. There are " 
                << outdegreeCounts.size() << " unique outdegrees and the maximum is " << maxOutdegreeFound << "." << std::endl;
        }
        //initialize this for if maxOutdegree = 0
        {
            sdsl::sd_vector_builder temp(bwt.size(), 0);
            this->isSmall = sdsl::sd_vector<>(temp);
        }
        this->largeRecords = bwt;

        SmallRecordArray tempSmallRecords;
        sdsl::sd_vector<> tempIsSmall;
        std::vector<CompressedRecord> tempLargeRecords;
        assert(maxOutdegreeFound < size_type(1) << size_type(32));
        for (size_type nextOutdegree = 0; nextOutdegree < maxOutdegreeFound + 1; ++nextOutdegree) {
            //find nextOutdegree to add
            if (outdegreeCounts.count(nextOutdegree) == 0)
                continue;
            double thisOutdegreeStart = gbwt::readTimer();
            if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
            {
                std::cerr << "lf_GBWT::GBWT::GBWT(): Computing if nodes with outdegree " << nextOutdegree << " should be stored in the smallRecordArray or in largeRecords."
                    << " There are " << std::get<1>(outdegreeCounts[nextOutdegree]) << " nodes with this outdegree. They have a total length of " << std::get<2>(outdegreeCounts[nextOutdegree])
                    << ". They take a total of " << std::get<0>(outdegreeCounts[nextOutdegree]) << " bytes." << std::endl;
            }
            if (nextOutdegree == 0) {
                //remember to add check for if max outdegree == 1, then every node in small records is empty,
                //no need to check emptyRecords array (it will not be initialized correctly)
                tempSmallRecords.maxOutdegree = 1;
                tempSmallRecords.records = std::get<1>(outdegreeCounts[nextOutdegree]);
            }
            else {
                //build smallRecordArray tempSmallRecords with maxOutdegree nextOutdegree+1
                tempSmallRecords.maxOutdegree = nextOutdegree + 1;

                size_type recordNum = 0, prevOutdegreePrefixSum = 0, prevLengthPrefixSum = 0, 
                          nonemptyRecordNum = 0;
                std::vector<size_type> outdegreePrefixSumAssist, alphabetAssist,
                    emptyRecordsAssist, lengthPrefixSumAssist, outgoingAssist,
                    firstAssist, firstByAlphabetAssist, alphabetByRunAssist,
                    firstByAlphCompAssist;
                for (size_type i = 0; i < this->effective(); ++i) {
                    const CompressedRecord& rec = bwt[i];
                    if (rec.outdegree() >= tempSmallRecords.maxOutdegree) {
                        largeRecords.push_back(rec);
                        continue;
                    }
                    if (rec.outdegree() == 0) {
                        emptyRecordsAssist.push_back(recordNum);
                    }
                    else {
                        gbwt::CompressedRecord gRec = source.record(source.toNode(i));
                        std::vector<std::vector<std::pair<size_type,size_type>>> localFirstByAlphabetAssist;
                        localFirstByAlphabetAssist.resize(tempSmallRecords.maxOutdegree);

                        size_type recLength = rec.size(), lengthPrefixSum = prevLengthPrefixSum;
                        outdegreePrefixSumAssist.push_back(prevOutdegreePrefixSum);
                        prevOutdegreePrefixSum += rec.outdegree();
                        lengthPrefixSumAssist.push_back(prevLengthPrefixSum);
                        prevLengthPrefixSum += recLength;
                        for (const auto& edge : rec.outgoing)
                            outgoingAssist.push_back(edge);
                        for (auto it = rec.alphabet.one_begin(); it != rec.alphabet.one_end(); ++it)
                            alphabetAssist.push_back((nonemptyRecordNum * this->effective()) + it->second);

                        gbwt::CompressedRecordFullIterator iter(gRec);
                        size_type start;
                        while(!iter.end()) {
                            start = iter.offset() - iter.run.second;
                            firstAssist.push_back(start + lengthPrefixSumAssist.back());
                            localFirstByAlphabetAssist[iter.run.first].emplace_back((tempSmallRecords.maxOutdegree * lengthPrefixSum) + (iter.run.first*recLength) + start, iter.run.second);
                            alphabetByRunAssist.push_back(iter.run.first);
                            ++iter;
                        }

                        size_type currInd = lengthPrefixSum;
                        for (const auto& alphArr : localFirstByAlphabetAssist) {
                            for (const auto& a : alphArr) {
                                firstByAlphabetAssist.push_back(a.first);
                                firstByAlphCompAssist.push_back(currInd);
                                currInd += a.second;
                            }
                        }
                        assert(currInd == prevLengthPrefixSum);
                        ++nonemptyRecordNum;
                    }
                    ++recordNum;
                }
                assert(recordNum - nonemptyRecordNum == emptyRecordsAssist.size());
                assert(this->effective() < (size_type(1) << size_type(32)));

                //build all structures using assist vectors
                auto buildSD = [] (sdsl::sd_vector<>& toBuild, const std::vector<size_type>& data, const size_type size, const size_type setBits) {
                    //assumes data sorted
                    //assert(std::is_sorted(data.begin(), data.end());
                    assert(size >= setBits);
                    assert(setBits == data.size());
                    sdsl::sd_vector_builder builder(size, setBits);
                    for (const size_type& x : data)
                        builder.set(x);
                    toBuild = sdsl::sd_vector<>(builder);
                };
                auto buildIntVec = [] (sdsl::int_vector<0>& toBuild, const std::vector<size_type>& data) {
                    toBuild.resize(data.size());
                    for (size_type i = 0; i < data.size(); ++i)
                        toBuild[i] = data[i];
                    sdsl::util::bit_compress(toBuild);
                };
                buildSD(tempSmallRecords.outDegreePrefixSum, outdegreePrefixSumAssist, prevOutdegreePrefixSum, nonemptyRecordNum);
                buildSD(tempSmallRecords.emptyRecords, emptyRecordsAssist, recordNum, recordNum - nonemptyRecordNum);
                buildSD(tempSmallRecords.prefixSum, lengthPrefixSumAssist, prevLengthPrefixSum, nonemptyRecordNum);
                buildIntVec(tempSmallRecords.outgoing, outgoingAssist);
                buildSD(tempSmallRecords.first, firstAssist, prevLengthPrefixSum, firstAssist.size());
                buildSD(tempSmallRecords.firstByAlphabet, firstByAlphabetAssist, prevLengthPrefixSum*tempSmallRecords.maxOutdegree, firstAssist.size());
                buildSD(tempSmallRecords.firstByAlphComp, firstByAlphCompAssist, prevLengthPrefixSum, firstAssist.size());
                buildSD(tempSmallRecords.alphabet, alphabetAssist, this->effective()*this->effective(), prevOutdegreePrefixSum);
                buildIntVec(tempSmallRecords.alphabetByRun, alphabetByRunAssist);
            }

            size_type prevIsSmallSize = sdsl::size_in_bytes(this->isSmall),
                      prevSmallRecordsSize = sdsl::size_in_bytes(this->smallRecords),
                      prevLargeRecordsSize = sdsl::size_in_bytes(this->largeRecords),
                      prevSize = prevIsSmallSize + prevSmallRecordsSize + prevLargeRecordsSize,
                      tempIsSmallSize = sdsl::size_in_bytes(tempIsSmall),
                      tempSmallRecordsSize = sdsl::size_in_bytes(tempSmallRecords),
                      tempLargeRecordsSize = sdsl::size_in_bytes(tempLargeRecords),
                      tempSize = tempIsSmallSize + tempSmallRecordsSize + tempLargeRecordsSize;

            if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
            {
                std::cerr << "lf_GBWT::GBWT::GBWT(): Finished splitting nodes in " << gbwt::readTimer() - thisOutdegreeStart << " seconds." << std::endl
                    << "Previously, the three structures (isSmall, smallRecords, largeRecords) took " << prevSize << " bytes in total, respectively ( " << prevIsSmallSize << ", " << prevSmallRecordsSize << ", " << prevLargeRecordsSize << ")." << std::endl
                    << "Now, the new structures take                                                " << tempSize << " bytes in total, respectively ( " << tempIsSmallSize << ", " << tempSmallRecordsSize << ", " << tempLargeRecordsSize << ")." << std::endl;
            }

            if (tempSize >= prevSize) {
                break;
            }

            this->isSmall = std::move(tempIsSmall);
            this->smallRecords = std::move(tempSmallRecords);
            this->largeRecords = std::move(tempLargeRecords);
            tempIsSmall = sdsl::sd_vector<>();
            tempSmallRecords = SmallRecordArray();
            tempLargeRecords = std::vector<CompressedRecord>();
        }
        
        

        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            double smallSeconds = gbwt::readTimer() - smallRecordsStart;
            std::cerr << "lf_GBWT::GBWT::GBWT(): Computed small and large record partition in " << smallSeconds << " seconds" << std::endl;

            double seconds = gbwt::readTimer() - start;
            std::cerr << "lf_GBWT::GBWT::GBWT(): Processed " << source.effective() << " nodes of total length " << this->size() << " in " << seconds << " seconds" << std::endl;
        }
    }

    GBWT::size_type GBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += sdsl::serialize(this->header, out, child, "header");
        written_bytes += sdsl::serialize(this->tags, out, child, "tags");
        written_bytes += sdsl::serialize(this->isSmall, out, child, "isSmall");
        written_bytes += sdsl::serialize(this->smallRecords, out, child, "smallRecords");
        written_bytes += sdsl::serialize(this->largeRecords, out, child, "largeRecords");
        if(this->hasMetadata())
        {
            written_bytes += sdsl::serialize(this->metadata, out, child, "metadata");
        }
        
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void GBWT::load(std::istream& in){
        sdsl::load(this->header, in);
        sdsl::load(this->tags, in);
        sdsl::load(this->isSmall, in);
        sdsl::load(this->smallRecords, in);
        sdsl::load(this->largeRecords, in);
        if(this->hasMetadata())
        {
            sdsl::load(this->metadata, in);
        }
    }

    std::pair<gbwt::size_type, gbwt::size_type> GBWT::runs() const {
        std::pair<size_type, size_type> result(0,0);
        for (const CompressedRecord& rec: this->largeRecords) {
            std::pair<size_type, size_type> temp = rec.runs();
            result.first  += temp.first;
            result.second += temp.second;
        }
        for (size_type i = 0; i < isSmall.ones(); ++i) {
            std::pair<size_type, size_type> temp = this->smallRecords.runs(i);
            result.first += temp.first;
            result.second += temp.second;
        }
        return result;
    }

    gbwt::vector_type GBWT::extract(size_type sequence) const {
        gbwt::vector_type ans;
        if (sequence >= this->sequences())
            return ans;
        gbwt::edge_type position = {gbwt::ENDMARKER, sequence};
        position = this->LF(position);
        while (position.first != gbwt::ENDMARKER){
            ans.push_back(position.first);
            position = this->LF(position);
        }
        ans.shrink_to_fit();
        return ans;
    }

    gbwt::edge_type GBWT::inverseLF(gbwt::node_type from, size_type i) const {
        if (!this->bidirectional() || from == gbwt::ENDMARKER) { return gbwt::invalid_edge(); }

        //find the predecessor node id
        gbwt::node_type revFrom = gbwt::Node::reverse(from);
        gbwt::node_type pred = this->predecessorAt(revFrom, i);
        if (pred == gbwt::invalid_node()) { return gbwt::invalid_edge(); }

        //determine the offset 
        auto ind = this->isSmallAndIndex(pred);
        size_type offset = (ind.first)? this->smallRecords.offsetTo(ind.second, this->toComp(from), i)
            : this->largeRecords[ind.second].offsetTo(this->toComp(from), i);
        if (offset == gbwt::invalid_offset()) { return {pred, offset}; return gbwt::invalid_edge(); }

        return { pred, offset};
    }

    CompressedRecord::size_type CompressedRecord::offsetTo(gbwt::comp_type to, size_type i) const {
        size_type outrank = this->edgeTo(to);
        if (outrank >= this->outdegree() || i < outgoing[outrank]) { return gbwt::invalid_offset(); }
        i -= outgoing[outrank];
        auto firstToRun = this->firstByAlphabet.successor(outrank*this->size());
        auto firstToCompRun = this->firstByAlphComp.select_iter(firstToRun->first+1);
        auto runComp = this->firstByAlphComp.predecessor(firstToCompRun->second + i);
        i -= (runComp->second - firstToCompRun->second);
        auto run = this->firstByAlphabet.select_iter(runComp->first+1);
        if (run->second/this->size() != outrank)
            return gbwt::invalid_offset();
        return (run->second-outrank*this->size()) + i;
    }
    
} // namespace lf_gbwt
#endif //GBWT_QUERY_LF_GBWT_H

