#ifndef GBWT_QUERY_LF_GBWT_H
#define GBWT_QUERY_LF_GBWT_H

//heavily derived/modified from FastLocate

#include<gbwt/gbwt.h>
#include<gbwt/internal.h>

namespace lf_gbwt{
    typedef std::uint8_t  byte_type;

    class GBWT;

    //all inputs and outputs to these functions are in the compressed alphabet (toComp)
    struct CompressedRecord
    {
        typedef gbwt::size_type size_type;

        //std::vector<gbwt::edge_type> outgoing;
        //outgoing[i] stores the value BWT.rank(v,i) where v is the current node
        sdsl::int_vector<0> outgoing;

        //const byte_type*       body;
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
        //sdsl::sd_vector<>::select_1_type firstByAlphCompSelect;

        //size_type              data_size;
        //bit vector of length gbwt.effective(), bit i is set if node i is contained in the BWT of this record
        sdsl::sd_vector<> alphabet;
        //sdsl::sd_vector<>::select_1_type alphabetSelect;

        //int i stores the mapped alphabet value of the BWT run
        sdsl::int_vector<0> alphabetByRun;

        CompressedRecord() = default;
        CompressedRecord(const gbwt::CompressedRecord &, const GBWT *);

        // Returns the size of the record corresponding to the given semiopen interval.
        //static size_type recordSize(const std::vector<byte_type>& source, size_type start, size_type limit);

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
        void load(std::istream& in);

        // Checks whether the record starting at the given position is empty.
        //static bool emptyRecord(const std::vector<byte_type>& source, size_type start);

        size_type size() const { return this->first.size(); }
        bool empty() const { return (this->size() == 0); }
        std::pair<size_type, size_type> runs() const; // (concrete, logical)
        size_type outdegree() const { return this->outgoing.size(); }

        // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
        gbwt::edge_type LF(size_type i) const;

        // Returns the predecessor node for the sequene at offset `i` in the other orientation
        // of this node, or `invalid_node()` if there is no such node.
        // This can be used for computing inverse LF in a bidirectional GBWT.
        //gbwt::node_type predecessorAt(size_type i) const;

        // Returns `offset` such that `LF(offset) == (to, i)`, or `invalid_offset()`
        // if there is no such offset.
        // This can be used for computing inverse LF in a bidirectional GBWT.
        size_type offsetTo(gbwt::node_type to, size_type i) const;

        // As above, but also reports the closed offset range ('run') and the identifier
        // ('run_id') of the logical run used for computing LF().
        gbwt::edge_type LF(size_type i, gbwt::range_type& run, size_type& run_id) const;

        // As above, but also sets 'run_end' to the last offset of the current logical run.
        gbwt::edge_type runLF(size_type i, size_type& run_end) const;

        // Returns invalid_offset() if there is no edge to the destination.
        size_type LF(size_type i, gbwt::comp_type to) const;

        // Returns Range::empty_range() if the range is empty or the destination is invalid.
        gbwt::range_type LF(gbwt::range_type range, gbwt::comp_type to) const;

        // As above, but also sets 'starts_with_to' if the range starts with node 'to', and
        // sets 'first_run' to the run identifier of the first run of to in overlapping
        // with the range.
        //gbwt::range_type LF(gbwt::range_type range, gbwt::node_type to, bool& starts_with_to, size_type& first_run) const;

        // As above, but also returns the number of characters x with
        // Node::reverse(x) < Node::reverse(to) in the range.
        //gbwt::range_type bdLF(gbwt::range_type range, gbwt::node_type to, size_type& reverse_offset) const;

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
        size_type offset(gbwt::rank_type outrank) const { return this->outgoing[outrank]; }

        //returns the run id by logical runs of an offset
        size_type logicalRunId(size_type i) const {
            if (i >= this->size())
                return gbwt::invalid_offset();
            size_type concRunId = this->first.predecessor(i)->first;
            if (!this->hasEdge(gbwt::ENDMARKER))
                return concRunId;
            size_type num0RunsBefore = this->firstByAlphabet.successor(i)->first;
            size_type num0Before = this->firstByAlphComp.select_iter(num0RunsBefore+1)->second;
            return concRunId - num0RunsBefore + num0Before;
        }
    };

    class GBWT
    {
        public:
            typedef CompressedRecord::size_type size_type;
//
//            //------------------------------------------------------------------------------
//
            GBWT() = default;
            explicit GBWT(const gbwt::GBWT& source);
//            GBWT(const DynamicGBWT& source);
//            GBWT(GBWT&& source);
//            ~GBWT();
//
//            // Merge the sources, assuming that node ids do not overlap.
//            // Also merges the metadata if all indexes contain it.
//            explicit GBWT(const std::vector<GBWT>& sources);
//
//            void swap(GBWT& another);
//            GBWT& operator=(const GBWT& source);
//            GBWT& operator=(const DynamicGBWT& source);
//            GBWT& operator=(GBWT&& source);
//
//            void resample(size_type sample_interval);
//
            size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
            void load(std::istream& in);
//
//            void simple_sds_serialize(std::ostream& out) const;
//            void simple_sds_load(std::istream& in);
//            size_t simple_sds_size() const;
//
            const static std::string EXTENSION; // .lfgbwt
//
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
            //size_type samples() const { return this->da_samples.size(); }

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

//            SearchState find(node_type node) const { return gbwt::find(*this, node); }
//
//            template<class Iterator>
//                SearchState find(Iterator begin, Iterator end) const { return gbwt::find(*this, begin, end); }
//
//            SearchState prefix(node_type node) const { return gbwt::prefix(*this, node); }
//
//            template<class Iterator>
//                SearchState prefix(Iterator begin, Iterator end) const { return gbwt::prefix(*this, begin, end); }
//
//            SearchState extend(SearchState state, node_type node) const { return gbwt::extend(*this, state, node); }
//
//            template<class Iterator>
//                SearchState extend(SearchState state, Iterator begin, Iterator end) const { return gbwt::extend(*this, state, begin, end); }
//
//            size_type locate(node_type node, size_type i) const { return gbwt::locate(*this, edge_type(node, i)); }
//            size_type locate(edge_type position) const { return gbwt::locate(*this, position); }
//
//            std::vector<size_type> locate(node_type node, range_type range) const { return this->locate(SearchState(node, range)); }
//            std::vector<size_type> locate(SearchState state) const;
//
            gbwt::vector_type extract(size_type sequence) const;
//            vector_type extract(edge_type position) const { return gbwt::extract(*this, position); }
//            vector_type extract(edge_type position, size_type max_length) const { return gbwt::extract(*this, position, max_length); }
//
//            //------------------------------------------------------------------------------
//
//            /*
//               Bidirectional search interface. The queries check that the parameters are valid.
//               On error or failed search, the return value is an empty bidirectional search state.
//            */
//
//            BidirectionalState bdFind(node_type node) const { return gbwt::bdFind(*this, node); }
//
//            BidirectionalState bdExtendForward(BidirectionalState state, node_type node) const { return gbwt::bdExtendForward(*this, state, node); }
//
//            BidirectionalState bdExtendBackward(BidirectionalState state, node_type node) const { return gbwt::bdExtendBackward(*this, state, node); }
//
//            //------------------------------------------------------------------------------
//
//            /*
//               Low-level interface: Nodes. The interface assumes that node identifiers are valid,
//               except in contains() / hasEdge(). This can be checked with contains().
//               */
//
            bool contains(gbwt::node_type node) const
            {
                return ((node < this->sigma() && node > this->header.offset) || node == gbwt::ENDMARKER);
            }
//
//            bool contains(edge_type position) const
//            {
//                return (this->contains(position.first) && position.second < this->nodeSize(position.first));
//            }
//
//            bool contains(SearchState state) const
//            {
//                return (this->contains(state.node) && !(state.empty()) && state.range.second < this->nodeSize(state.node));
//            }
//
//            bool hasEdge(node_type from, node_type to) const
//            {
//                return (this->contains(from) && this->record(from).hasEdge(to));
//            }
//
//            std::vector<edge_type> edges(node_type from) const
//            {
//                return this->record(from).outgoing;
//            }
//
            gbwt::node_type firstNode() const { return this->header.offset + 1; }
            gbwt::comp_type toComp(gbwt::node_type node) const { return (node == 0 ? node : node - this->header.offset); }
            gbwt::node_type toNode(gbwt::comp_type comp) const { return (comp == 0 ? comp : comp + this->header.offset); }

            size_type nodeSize(gbwt::node_type node) const { return this->record(node).size(); }
            bool empty(gbwt::node_type node) const { return this->nodeSize(node) == 0; }
//
//            //------------------------------------------------------------------------------
//
//            /*
//               Low-level interface: Navigation and searching. The interface assumes that node
//               identifiers are valid. This can be checked with contains().
//               */
//
//            // On error: invalid_edge().
//            edge_type LF(node_type from, size_type i) const
//            {
//                if(from == ENDMARKER) { return this->endmarker().LF(i); }
//                return this->record(from).LF(i);
//            }
//
            // On error: invalid_edge().
            gbwt::edge_type LF(gbwt::edge_type position) const
            {
                //if(position.first == ENDMARKER) { return this->endmarker().LF(position.second); }
                gbwt::edge_type ans = this->record(position.first).LF(position.second);
                ans.first = this->toNode(ans.first);
                return ans;
            }
            
            gbwt::edge_type LF(gbwt::edge_type position, gbwt::range_type& run, gbwt::size_type& run_id) const {
                gbwt::edge_type next = this->record(position.first).LF(position.second, run, run_id);
                next.first = this->toNode(next.first);
                return next;
            }
//
//            // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
//            // On error: invalid_edge().
            gbwt::edge_type inverseLF(gbwt::node_type from, size_type i) const;

            // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
            // On error: invalid_edge().
            gbwt::edge_type inverseLF(gbwt::edge_type position) const
            {
                return this->inverseLF(position.first, position.second);
            }
            
            gbwt::node_type predecessorAt(gbwt::node_type, size_type i) const;
//
//            // On error: invalid_offset().
//            size_type LF(node_type from, size_type i, node_type to) const
//            {
//                return this->record(from).LF(i, to);
//            }
//
//            // On error: invalid_offset().
//            size_type LF(edge_type position, node_type to) const
//            {
//                return this->record(position.first).LF(position.second, to);
//            }
//
//            // On error: Range::empty_range().
//            range_type LF(node_type from, range_type range, node_type to) const
//            {
//                return this->record(from).LF(range, to);
//            }
//
//            // On error: Range::empty_range().
//            range_type LF(SearchState state, node_type to) const
//            {
//                return this->record(state.node).LF(state.range, to);
//            }
//
//            // On error: Range::empty_range().
//            range_type bdLF(SearchState state, node_type to, size_type& reverse_offset) const
//            {
//                return this->record(state.node).bdLF(state.range, to, reverse_offset);
//            }
//
//            //------------------------------------------------------------------------------
//
//            /*
//               Low-level interface: Sequences. The interface assumes that node identifiers are
//               valid. This can be checked with contains().
//               */
//
//            // Starting position of the sequence or invalid_edge() if something fails.
//            edge_type start(size_type sequence) const { return this->LF(ENDMARKER, sequence); }
//
//            // Returns the sampled document identifier or invalid_sequence() if there is no sample.
//            size_type tryLocate(node_type node, size_type i) const
//            {
//                return this->da_samples.tryLocate(this->toComp(node), i);
//            }
//
//            // Returns the sampled document identifier or invalid_sequence() if there is no sample.
//            size_type tryLocate(edge_type position) const
//            {
//                return this->da_samples.tryLocate(this->toComp(position.first), position.second);
//            }
//
            //------------------------------------------------------------------------------

            gbwt::GBWTHeader  header;
            gbwt::Tags        tags;
            std::vector<CompressedRecord> bwt;
            //gbwt::RecordArray bwt;
            //gbwt::DASamples   da_samples;
            gbwt::Metadata    metadata;

            // Decompress and cache the endmarker, because decompressing it is expensive.
            //CompressedRecord *endmarker_record;

            //------------------------------------------------------------------------------
//
//            /*
//               Internal interface. Do not use.
//               */
//
//        private:
//            void copy(const GBWT& source);
//            void resetTags();
//            void addSource();
//            void cacheEndmarker();
//
//        public:
            const CompressedRecord& record(gbwt::node_type node) const { 
                if (!this->contains(node))
                    throw std::invalid_argument("node given to lf_gbwt.record invalid!");
                return this->bwt[this->toComp(node)]; 
            }
//            const DecompressedRecord& endmarker() const { return this->endmarker_record; }

            bool verify(const gbwt::GBWT& g) {
                //std::cout << "In LF GBWT verify()" << std::endl;
                bool equal = (this->size() == g.size()) && (this->sequences() == g.sequences()) && (this->sigma() == g.sigma()) && (this->effective() == g.effective()) && (this->runs() == g.runs()) && (this->bidirectional() == g.bidirectional());
                //std::cout << "equal " << equal << std::endl;
                equal = equal && this->verifyBWT(g);
                //std::cout << "equal " << equal << std::endl;
                equal = equal && this->verifyLF(g);
                //std::cout << "equal " << equal << std::endl;
                equal = equal && (!this->bidirectional() || this->verifyInverseLF(g));
                //std::cout << "equal " << equal << std::endl;
                return equal;
            }

            //assumes number of nodes and node ids are equivalent
            bool verifyBWT(const gbwt::GBWT& g) {
                //std::cout << "in verifyBWT" << std::endl;
                bool equal = true;
                #pragma omp parallel for schedule(dynamic, 1)
                for (gbwt::comp_type i = 0; i < this->effective(); ++i) {
                    gbwt::node_type node = this->toNode(i);
                    bool nodeEqual = node == g.toNode(i);
                    gbwt::CompressedRecord rec = g.record(node);
                    const CompressedRecord& lfrec = this->record(node);
                    nodeEqual = nodeEqual && rec.size() == lfrec.size();
                    //std::cout << "before node " << node << " bwt verification nodeEqual "<< nodeEqual  << std::endl;
                    for (gbwt::size_type j = 0; nodeEqual && j < lfrec.size(); ++j)
                        nodeEqual = nodeEqual && rec[j] == this->toNode(lfrec[j]);
                    #pragma omp critical 
                    {
                       equal = equal && nodeEqual;
                    }
                }
                //std::cout << "after bwt verification, equal " << equal << std::endl;
                return equal;
            }

            bool verifyLF(const gbwt::GBWT& g) {
                //std::cout << "in verifyLF" << std::endl;
                bool equal = true;
                #pragma omp parallel for schedule(dynamic, 1)
                for (gbwt::comp_type i = 0; i < this->effective(); ++i) {
                    gbwt::node_type node = this->toNode(i);
                    bool nodeEqual = node == g.toNode(i);
                    gbwt::size_type nodeSize = g.nodeSize(node);
                    nodeEqual = nodeEqual && nodeSize == this->nodeSize(node);
                    //std::cout << "before node " << node << " bwt verification nodeEqual "<< nodeEqual  << std::endl;
                    for (gbwt::size_type j = 0; nodeEqual && j < nodeSize; ++j)
                        nodeEqual = nodeEqual && this->LF({node, j}) == g.LF({node, j});
                    #pragma omp critical 
                    {
                       equal = equal && nodeEqual;
                    }
                }
                //std::cout << "after LF verification, equal " << equal << std::endl;
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

    CompressedRecord::CompressedRecord(const gbwt::CompressedRecord & rec, const GBWT *source){
        //std::cout << "Compressed Record initialized, effectiveAlphSize " << source->effective();
        size_type n = rec.size(), sigma = rec.outgoing.size(), runs = rec.runs().first;
        //std::cout << " n " << n << " sigma " << sigma << " runs " << runs << std::endl;
        if (n == 0)
            return;

        //outgoing
        this->outgoing.resize(sigma);
        for (unsigned i = 0; i < sigma; ++i)
            this->outgoing[i] = rec.outgoing[i].second;
        sdsl::util::bit_compress(this->outgoing);

        //std::cout << "outgoing construction done" << std::endl;

        //alphabet
        sdsl::sd_vector_builder alphabetBuilder(source->effective(), sigma);
        for (unsigned i = 0; i < sigma; ++i)
            alphabetBuilder.set(source->toComp(rec.outgoing[i].first));
        this->alphabet = sdsl::sd_vector<>(alphabetBuilder);

        //std::cout << "alphabet construction done" << std::endl;

        //first, firstByAlphabet, firstByAlphComp, and alphabetByRun;
        sdsl::sd_vector_builder 
            firstBuilder(n, runs),
            firstByAlphabetBuilder(n*sigma, runs),
            firstByAlphCompBuilder(n, runs);
        this->alphabetByRun.resize(runs);

        //std::vector<std::pair<size_type,size_type>> firstByAlphabetAssist;
        //start location of run in firstByAlphabetAssist and length of run
        std::vector<std::vector<std::pair<size_type,size_type>>> firstByAlphabetAssist;
        firstByAlphabetAssist.resize(sigma);

        gbwt::CompressedRecordFullIterator iter(rec);
        size_type start, runInd = 0;
        while (!iter.end()){
            start = iter.offset() - iter.run.second;
            firstBuilder.set(start);
            firstByAlphabetAssist[iter.run.first].emplace_back(iter.run.first*n + start, iter.run.second);
            //std::cout << "runInd " << runInd << " runId " << iter.runId() << std::endl;
            this->alphabetByRun[runInd] = iter.run.first;
            ++iter;
            ++runInd;
        }
        this->first = sdsl::sd_vector<>(firstBuilder);
        sdsl::util::bit_compress(this->alphabetByRun);

        //std::cout << "first and alphabetByRun construction done" << std::endl;

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
        if (!this->hasEdge(gbwt::ENDMARKER)) // (outgoing.size() == 0 || this->successor(0) != 0) 
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
        //size_type outrank = this->compAlphabetAt(i);
        //auto run = this->firstByAlphabet.predecessor(outrank*this->size() + i);
        //auto firstOutrankRun = this->firstByAlphabet.successor(outrank*this->size());
        //size_type numOutrank = firstByAlphCompSelect(run->first+1) - firstByAlphCompSelect(firstOutrankRun->first+1)
        //    + (i - (run->second-(outrank*this->size())));
        //return {this->alphabetSelect(outrank+1), outgoing[outrank] + numOutrank};
        gbwt::comp_type next = (*this)[i];
        //std::cout << "{IN LF, next = " << next << ", LF = " << this->LF(i, next) << "}" << std::flush;
        return {next, this->LF(i, next)};
    }

    gbwt::node_type GBWT::predecessorAt(gbwt::node_type revFrom, size_type i) const {
        auto revNode = [this](gbwt::comp_type x) { 
            if (x == gbwt::ENDMARKER) 
                return x;
            return gbwt::Node::reverse(this->toNode(x));
        };
        const CompressedRecord & rev = this->record(revFrom);
        if (i >= rev.size())
            return gbwt::invalid_node();
        size_type predoutrank = rev.firstByAlphabet.select_iter(
                rev.firstByAlphComp.predecessor(i)->first+1
                )->second/rev.size();
        //std::cout << "(predoutrank " << predoutrank << ")" << std::flush;

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

        //std::cout << "(predoutrank " << predoutrank << ")" << std::flush;
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

        //std::cout << "(predoutrank " << predoutrank << ")" << std::flush;
        return revNode(alphIter->second);
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

    // As above, but also sets 'run_end' to the last offset of the current logical run.
    //gbwt::edge_type CompressedRecord::runLF(size_type i, size_type& run_end) const {
        //if (i >= this->size())
        //    return gbwt::invalid_edge();
        //auto run_iter = this->first.predecessor(i);
        //size_type outrank = this->alphabetByRun[run_iter->first];
        //auto nextRun = this->firstByAlphabet.successor(outrank*this->size() + i);
        //auto firstOutrankRun = this->firstByAlphabet.successor(outrank*this->size());
        //size_type numOutrank = this->firstByAlphComp.select_iter(nextRun->first+1)->second 
        //    - this->firstByAlphComp.select_iter(firstOutrankRun->first+1)->second
        //    + (i - (nextRun->second-(outrank*this->size())));

        //run_end = (++run_iter)->second-1;

        //return {this->alphabet.select_iter(outrank+1)->second, this->outgoing[outrank] + numOutrank};
    //}

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
        //std::cout << "{IN LF, numOutrank: " << numOutrank << " ";
        if ((*this)[i] == to){
            size_type afteriInRun = this->first.successor(i)->second - i;
            numOutrank -= afteriInRun;
        }
            //+ (i - (nextRun->second-(outrank*this->size())));
        //std::cout << numOutrank << "}" << std::flush;
        
        return this->outgoing[outrank] + numOutrank;
    }

    // Returns Range::empty_range() if the range is empty or the destination is invalid.
    /*gbwt::range_type CompressedRecord::LF(gbwt::range_type range, gbwt::comp_type to) const {
        if (gbwt::Range::empty(range) || !this->hasEdge(to))
            return gbwt::Range::empty_range();
        //return {this->LF(range.first, to), this->LF(range.second-1, to)}; //wrong
    }*/


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
        this->bwt.resize(source.effective());
        source.bwt.forEach([this] (size_type a, const gbwt::CompressedRecord& rec) {this->bwt[a] = CompressedRecord(rec, this);});
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            double seconds = gbwt::readTimer() - bwtStart;
            std::cerr << "lf_GBWT::GBWT::GBWT(): Constructed " << source.effective() << " nodes of total length " << this->size() << " in " << seconds << " seconds" << std::endl;
        }

        //this->bwt.shrink_to_fit();
        //this->endmarker_record = &(this->bwt[0]);
        if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
        {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "lf_GBWT::GBWT::GBWT(): Processed " << source.effective() << " nodes of total length " << this->size() << " in " << seconds << " seconds" << std::endl;
        }
    }

    GBWT::size_type GBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += sdsl::serialize(this->header, out, child, "header");
        written_bytes += sdsl::serialize(this->tags, out, child, "tags");
        written_bytes += sdsl::serialize(this->bwt, out, child, "bwt");
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
        sdsl::load(this->bwt, in);
        if(this->hasMetadata())
        {
            sdsl::load(this->metadata, in);
        }
    }

    std::pair<gbwt::size_type, gbwt::size_type> GBWT::runs() const {
        std::pair<size_type, size_type> result(0,0);
        for (size_type i = 0; i < this->effective(); ++i){
            std::pair<size_type, size_type> temp = this->record(this->toNode(i)).runs();
            result.first  += temp.first; 
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
            //std::cout << "(" << position.first << "," << position.second << ") " << std::flush;
            position = this->LF(position);
        }
        ans.shrink_to_fit();
        return ans;
    }

    gbwt::edge_type GBWT::inverseLF(gbwt::node_type from, size_type i) const {
        if (!this->bidirectional() || from == gbwt::ENDMARKER) { return gbwt::invalid_edge(); }

        //std::cout << "{inverseLF finding pred: " << std::flush;
        //find the predecessor node id
        gbwt::node_type revFrom = gbwt::Node::reverse(from);
        gbwt::node_type pred = this->predecessorAt(revFrom, i);
        if (pred == gbwt::invalid_node()) { return gbwt::invalid_edge(); }

        //std::cout << pred << " finding offset: " << std::flush;
        //determine the offset 
        size_type offset = this->record(pred).offsetTo(this->toComp(from), i);
        if (offset == gbwt::invalid_offset()) { return {pred, offset}; return gbwt::invalid_edge(); }

        //std::cout << offset << "}" << std::flush;
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
//
//    //------------------------------------------------------------------------------
//
//    void printStatistics(const GBWT& gbwt, const std::string& name, std::ostream& out = std::cout);
//    std::string indexType(const GBWT&);
//
//    //------------------------------------------------------------------------------
    
} // namespace lf_gbwt
#endif //GBWT_QUERY_LF_GBWT_H

