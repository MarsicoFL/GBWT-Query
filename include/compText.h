#ifndef GBWT_QUERY_COMPTEXT_H
#define GBWT_QUERY_COMPTEXT_H

#include"fast_lcp.h"
#include<vector>
#include<iostream>

//implementation of the compressed data structure by Gagie et. al. that
//provides random access to l consecutive characters of the text in O(log(n/r) + l) time 
//in O(r log(n/r)) space. We modify it to only answer point queries instead of range queries.
//This changes the number of half blocks per sample per level from 7 to 4 (4 to 2 that are explicitly stored).
//Furthermore, in implementation we use GBWT logical runs, r_GBWT, instead of BWT runs, r, although 
//r_GBWT <= r+\Sigma <= 2r, therefore this implementation provides random access to character i
//of the text in O(log(n/r)) time in O(r log(n/r)) space. We also attempt to prune the structure
//during construction to reduce space.

class CompText{
    public:
    typedef gbwt::size_type size_type;
    private:
    //stores text in the GBWT compressed node format
    struct Level{
        typedef CompText::size_type size_type;
        //if the i-th and i+1-th blocks in this level corresponds to position k in the text that was sampled, and s_l is the length of a half block in this level,
        //then the i-th block corresponds to T[k-s_l, k) and the i+1-th block corresponds to T[k, k+s_l)
        //then the half blocks corresponding to the   i-th block are the     2i-th and the     2i+1-th half blocks, with ranges T[k-s_l, k-s_l+s_{l+1} = k-s_{l+1}) and T[k-s_l+s{l+1} = k-s_{l+1}, k    ) respectively.
        // and the half blocks corresponding to the i+1-th block are the 2(i+1)-th and the 2(i+1)+1-th half blocks, with ranges T[    k, k+s_l-s_{l+1} = k+s_{l+1}) and T[k+s_l-s{l+1} = k+s_{l+1}, k+s_l) respectively. 
        //where s_{l+1} = s{l}/2.
        //
        //
        //the i-th half-block in this level maps to a primary occurrence that contains the mapsTo[i]-th sample in the next level 
        //the starting position, i', of this primary occurrence is i' = j - s_{l+1} + offset
        //where j is the position of the mapsTo[i]-th sample in the text
        sdsl::int_vector<0> mapsTo;
        sdsl::int_vector<0> offset;

        void resize(size_type a) { mapsTo.resize(a); offset.resize(a); }
        void compress() { sdsl::util::bit_compress(mapsTo); sdsl::util::bit_compress(offset); }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type bytes = 0;
            bytes += this->mapsTo.serialize(out, child, "mapsTo");
            bytes += this->offset.serialize(out, child, "offset");
            sdsl::structure_tree::add_size(child, bytes);
            return bytes;
        }

        void load(std::istream& in) {
            sdsl::load(this->mapsTo, in);
            sdsl::load(this->offset, in);
        }
    };


    size_type s_0;
    public:
    const FastLCP* source;
    private:
    sdsl::sd_vector<> pathStarts;
    Level firstLevel;
    std::vector<Level> levels;
    //block size of firstLevel and levels[0]
    sdsl::int_vector<0> textSamples;

    public:
    CompText() = default;
    explicit CompText(const FastLCP &);
    void buildFullMem(const FastLCP &);
    gbwt::node_type at(size_type) const;
    gbwt::node_type atFLsuff(size_type suff) const { return this->at(this->FLsuffToTrueSuff(suff)); }
    size_type textLength() const { return this->pathStarts.size(); }

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;
        bytes += sdsl::serialize(this->s_0, out, child, "s_0");
        bytes += sdsl::serialize(this->pathStarts, out, child, "pathStarts");
        bytes += sdsl::serialize(this->firstLevel, out, child, "firstLevel");
        bytes += sdsl::serialize(this->levels, out, child, "levels");
        bytes += sdsl::serialize(this->textSamples, out, child, "textSamples");
        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load (std::istream& in) {
        sdsl::load(this->s_0, in);
        sdsl::load(this->pathStarts, in);
        sdsl::load(this->firstLevel, in);
        sdsl::load(this->levels, in);
        sdsl::load(this->textSamples, in);
    }

    const static std::string EXTENSION; //.compText
    //for debugging. remove later.
    //std::vector<size_type> temp_sl;
    /*
    gbwt::vector_type TrueSl(const gbwt::vector_type & fullText) const {
        gbwt::vector_type t;
        size_type s_l = this->s_0 >> levels.size();
        for (size_type i = 0; i < temp_sl.size(); ++i){
            //std::cout << temp_sl[i] << ':';
            for (size_type j = temp_sl[i] - s_l; j - temp_sl[i] + s_l < 2*s_l; ++j)
                if (j < fullText.size())
                    t.push_back(fullText[j]);
                else 
                    t.push_back(0);
        }
        return t;
    }
    */
    template<class T> bool compareTextSamples(const T & toCompare) const {
        if (toCompare.size() != this->textSamples.size())
            return false;
        for (size_type i = 0; i < toCompare.size(); ++i)
            if (toCompare[i] != this->source->rindex->index->toNode(this->textSamples[i]))
                return false;
        return true;
    }

    bool verifyText() const {
        bool good = true;
        #pragma omp parallel for schedule(dynamic, 1)
        for (gbwt::size_type i = 0; i < this->source->rindex->index->sequences(); ++i){
            //gbwt::size_type start = gbwt::readTimer();
            bool seqGood = true;
            gbwt::edge_type position = {gbwt::ENDMARKER, i};
            gbwt::size_type location = FLsuffToTrueSuff(i+1, 0);
            do {
                seqGood = seqGood && (position.first == this->at(--location));
                position = this->source->rindex->index->LF(position);
            } while (location != size_type(-1) && position.first != gbwt::ENDMARKER);
            #pragma omp critical
            {
                good = good && seqGood;
                //std::cout << "seq i done in " << gbwt::readTimer() - start << " seconds" << std::endl;
            }
        }
        //do {
        //    if (position.first != this->at(--ind))
        //        return false;
        //    this->fullLF(position, sequence);
        //} while (ind != 0 && sequence != size_type(-1));
        //return ind == 0 && sequence == size_type(-1);
        return good;
    }




    size_type FLsuffToTrueSuff(const size_type suff) const { return this->FLsuffToTrueSuff(this->source->rindex->seqId(suff), this->source->rindex->seqOffset(suff)); }
    private:
    size_type FLsuffToTrueSuff(const std::pair<size_type, size_type> suff) const { return this->FLsuffToTrueSuff(suff.first, suff.second); }
    size_type FLsuffToTrueSuff(const size_type seqId, const size_type seqOffset) const { return seqOffset + pathStarts.select_iter(seqId+1)->second; }

    void fullLF (gbwt::edge_type & position, size_type & sequence) const { 
        position = this->source->rindex->index->LF(position);
        if (position.first == gbwt::ENDMARKER) {
            --sequence; //currently, if LF on suffix 0 is called, then sequence becomes gbwt::invalid_offset(). This is purposeful, although this can be easily checked for and made to become sequences-1, becoming a circular bwt
            position = {gbwt::ENDMARKER, sequence};
        }
    };

};

const std::string CompText::EXTENSION = ".compText"; //.compText

CompText::CompText(const FastLCP & l): source(&l) {
    double start = gbwt::readTimer();
    //unpruned construction
    if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    {
        std::cerr << "CompText::CompText(): Unpruned construction of compressed text data structure" << std::endl;
    }
    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Computing path starts" << std::endl;
        }
        //compute path starts 
        std::vector<size_type> pathLengthPrefixSums;
        size_type prefixSum = 0, first = l.rindex->locateFirst(gbwt::ENDMARKER);
        for (size_type i = 0; i < l.rindex->index->sequences()-1; ++i){
            pathLengthPrefixSums.push_back(prefixSum);
            prefixSum += l.rindex->seqOffset(first) + 1; // + 1 for endmarker character
            first = l.rindex->locateNext(first);
        }
        pathLengthPrefixSums.push_back(prefixSum);
        prefixSum += l.rindex->seqOffset(first) + 1; // + 1 for endmarker character
        sdsl::sd_vector_builder builder(prefixSum, pathLengthPrefixSums.size());
        for (const auto & a : pathLengthPrefixSums)
            builder.set(a);
        this->pathStarts = sdsl::sd_vector<>(builder);
    }
    std::vector<size_type> sampleLocations;
    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Computing sample locations" << std::endl;
        }
        //compute sorted sampleLocations by true suffix values
        std::vector<size_type> temp = {0, this->pathStarts.size()-1};

        //add samples for beginning and end of concrete runs
        for (const auto & a : this->source->rindex->samples)
            temp.push_back(this->FLsuffToTrueSuff(a));
        for (const auto & a : this->source->samples_bot)
            temp.push_back(this->FLsuffToTrueSuff(a));

        //possible bug if samples and samples_bot don't include samples for logical runs,only concrete runs
        //i.e. if bwt[i] is endmarker, samples need to be stored
        //if missing, add samples for endmarkers in GBWT below

        gbwt::parallelQuickSort(temp.begin(), temp.end());
        for (const auto & a : temp)
            if (sampleLocations.size() == 0 || a != sampleLocations.back())
                sampleLocations.push_back(a);
//        std::cout << "sampleLocations ";
//           for (auto a:  sampleLocations) 
//              std::cout << a << ' ';
//        std::cout << std::endl;
    }
    //this->temp_sl = sampleLocations;


    //sizeEstimate[i] is an estimate for the size of the whole data structure if it had i intermediate levels
    //(if levels.size() was i)
    std::vector<size_type> sizeEstimate;
    
    auto ceilDiv = [] (size_type a, size_type b) { return (a/b) + (a%b != 0); };
    size_type r = l.samples_bot.size(), n = this->pathStarts.size();
    if (n == r) { throw std::invalid_argument("Number of characters and runs in text are the same, compressing won't work"); }
    this->s_0 = 1 << sdsl::bits::hi(ceilDiv(n,r));
    if (this->s_0 < ceilDiv(n,r)) { this->s_0 *= 2; }
    if (gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    { 
        std::cerr << "CompText::CompText(): The text of the GBWT has length " << n << " and " << r << " logical runs. Average logical run length is " << double(n)/r 
            << ". The smallest power of 2 larger than or equal to the average run length is " << this->s_0 << ". This is the width of the blocks in the first level of the compressed text data structure." << std::endl;
    }

    //if the i-th sample is at position k (sampleLocations[i] = k), 
    //blockPositionStarts[j][4i]   stores the position of the suffix T[k-s_{j+1} ... n]
    //blockPositionStarts[j][4i+1] stores the position of the suffix T[k         ... n]
    //blockPositionStarts[j][4i+2] stores the position of the suffix T[k+s_{j+1} ... n]
    //blockPositionStarts[j][4i+3] stores the position of the suffix T[k+s_{j}   ... n]
    //for j >= 0 s.t. s{j+1} >= 1
    std::vector<std::vector<gbwt::edge_type>> blockPositionStarts(sdsl::bits::hi(s_0));
    for (size_type i = 0; i < blockPositionStarts.size(); ++i)
        blockPositionStarts[i].resize(4*sampleLocations.size());
    if (gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    {
        std::cerr << "CompText::CompText(): The maximum number of levels in the tree is " << blockPositionStarts.size() << ". This is the log_2 of the first level's block width." << std::endl;
    }

    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Building block tree" << std::endl;
        }
        //build block tree 

        //auto ithStart = [] (size_type i, size_type n, size_type r) {
        //    //assumes 0 <= i <= r <= n
        //    return i*(n/r) + ((i*(n%r))/r);
        //};

        double firstLevelStart = gbwt::readTimer();
        {
            if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Building the first level of the block tree" << std::endl;
            }
            //build first level
            size_type topLevelRuns = ceilDiv(n, this->s_0), end = n, start = (topLevelRuns-1)*this->s_0;
            this->firstLevel.resize(topLevelRuns);

            size_type sequence = this->pathStarts.ones()-1;
            gbwt::edge_type position = {gbwt::ENDMARKER, sequence};

            gbwt::vector_type T(this->s_0);
            std::vector<bool> sampleUsed(sampleLocations.size(), false);
            size_type usedCount = 0;
            size_type runStart = gbwt::readTimer();
            for (size_type i = topLevelRuns; i != 0; end = start, start -= this->s_0){
                --i;
                if (i%10000 == 0) {
                    std::cout << "On run " << i << " after " << gbwt::readTimer() - runStart << " seconds" << std::endl;
                    runStart = gbwt::readTimer();
                }
//                std::cout << "start " << start << " end " << end << std::endl << "sampleLocations.size() " << sampleLocations.size() << std::endl;
//                std::cout << "run " << i << std::endl;
                //T.clear();
                for (size_type u = end; u != start;){
                    --u;
                    //T is reversed to allow searching using gbwt::GBWT::find
                    T[end-1-u] = position.first;
                    for (size_type j = 0, s_j = this->s_0; s_j != 1; ++j, s_j /= 2){
                        //VERIFY CORRECTNESS OF THIS
                        //check if current position can fill any values of blockPositionStarts
                        auto safeMinus = [] (size_type a, size_type b) { return (a > b)? a-b : 0; };
                        size_type sjplus1 = s_j/2;
                        auto it = sampleLocations.begin();
                        it = std::lower_bound(it, sampleLocations.end(), safeMinus(u, s_j));
                        if (it != sampleLocations.end() && u == *it+s_j    )
                            blockPositionStarts[j][4*(it - sampleLocations.begin()) + 3] = position;
                        it = std::lower_bound(it, sampleLocations.end(), safeMinus(u, sjplus1));
                        if (it != sampleLocations.end() && u == *it+sjplus1)
                            blockPositionStarts[j][4*(it - sampleLocations.begin()) + 2] = position;
                        it = std::lower_bound(it, sampleLocations.end(), u        );
                        if (it != sampleLocations.end() && u == *it        )
                            blockPositionStarts[j][4*(it - sampleLocations.begin()) + 1] = position;
                        it = std::lower_bound(it, sampleLocations.end(), u+sjplus1);
                        if (it != sampleLocations.end() && u == *it-sjplus1)
                            blockPositionStarts[j][4*(it - sampleLocations.begin())    ] = position;
                    }

                    this->fullLF(position, sequence);
                }
//                std::cout << "LF done " << std::endl;
                //assign block i to a run, block i corresponds to [i*s_0, (i+1)*s_0) = [start, end)
                //where s_0 is the smallest power of 2 larger than or equal to \ceiling(n/r) (and is the block size of the first level of the tree
                //where n is the total length and r is the number of GBWT logical runs
                auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                if (it != sampleLocations.end() && *it < end){
//                    std::cout << "no LF needed, marker found" << std::endl;
                    firstLevel.mapsTo[i] = it - sampleLocations.begin();
                    firstLevel.offset[i] = start + this->s_0 - *it;
                    usedCount += !sampleUsed[firstLevel.mapsTo[i]];
                    sampleUsed[firstLevel.mapsTo[i]] = true;
//                    std::cout << " mapsto " << firstLevel.mapsTo[i] << " offset " << firstLevel.offset[i] << std::endl;
                    continue;
                }
                //all endmarkers  are sampled therefore [start, end) contains no endmarkers
                //if [start, end) is not a primary occurrence of T[start, end).
                //Here we search for a primary occurence

//                std::cout << "LF needed, marker not found" << std::endl;
                    
                //assert(this->s_0 == end-start);
                size_type first = 0;
//                std::cout << "looking for T["<< start << ',' << end << ") = ";
//                    for (auto a : T)
//                       std::cout << a << ' ';
//                std::cout  << " in GBWT" << std::endl;
                gbwt::SearchState range = l.rindex->find(T.begin(), T.end(), first);

                for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                    std::cout << "found in range, seqId " << l.rindex->seqId(a) << " offset " << l.rindex->seqOffset(a);
                    a = this->FLsuffToTrueSuff(a);
//                    std::cout << " true suff " << a << std::endl; 
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                    if (it != sampleLocations.end() && *it < a + (end-start)){
//                        std::cout << " is in range, *it " << *it;
                        firstLevel.mapsTo[i] = it - sampleLocations.begin();
                        firstLevel.offset[i] = a + this->s_0 - *it;
//                        std::cout << " mapsto " << firstLevel.mapsTo[i] << " offset " << firstLevel.offset[i] << std::endl;
                        usedCount += !sampleUsed[firstLevel.mapsTo[i]];
                        sampleUsed[firstLevel.mapsTo[i]] = true;
                        break;
                    }
                }
            }
            //need to compress mapsto to rank using sampleUnused here. Do once compressing block tree
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Out of " << sampleLocations.size() << " sampled locations, " << usedCount << " unique blocks were mapped to by the " << topLevelRuns << " blocks in the first level of the block tree." << std::endl;
            }
        }
//        std::cout << "firstLevel ";
//        for (size_type i = 0; i < firstLevel.mapsTo.size(); ++i)
//            std::cout << '(' << firstLevel.mapsTo[i] << ',' << firstLevel.offset[i] << ") ";
//        std::cout << std::endl;
//        std::cout << "before compression firstLevel size " << std::flush << sdsl::size_in_bytes(firstLevel) << std::endl;
        firstLevel.compress();

        //build the rest of the levels
        
        size_type treeAndFirstLevelSize = sdsl::size_in_bytes(this->pathStarts) + 
                sdsl::size_in_bytes(firstLevel);
        size_type textSamplesEstimateBits = 2*s_0*firstLevel.mapsTo.size()*sdsl::bits::length(l.rindex->index->effective()); //+64?
        sizeEstimate.push_back(treeAndFirstLevelSize + ceilDiv(textSamplesEstimateBits, 8));
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): First level takes " << treeAndFirstLevelSize << " bytes. If the tree had no more levels, the text samples would take roughly " << ceilDiv(textSamplesEstimateBits, 8) << " bytes." 
                << " This is a total of " << sizeEstimate[0] << " bytes if the tree was truncated at this level." << std::endl;
            std::cerr << "CompText::CompText(): Construction of the first level took " << gbwt::readTimer() - firstLevelStart << " seconds" << std::endl;
        }
        std::vector<bool> sampleUsed(sampleLocations.size());
        size_type usedCount = 0;
        for(size_type l = 0, s_l = s_0; s_l != 1; ++l, s_l /= 2){
            double levelStartTime = gbwt::readTimer();
            sampleUsed.assign(sampleLocations.size(), false);
            usedCount = 0;
            //compute intermediate level l (levels[l])
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Computing level " << l << " of the tree. Blocks at this level have a width of " << s_l << std::endl;
            }
            size_type slplus1 = s_l/2;

            //since this is the unpruned construction, half blocks adjacent to samples are
            //not stored explicitly, they are implicitly mapped to the same sample in the below level
            //therefore each level has 2*samples half blocks
            this->levels.push_back(Level());
            this->levels.back().resize(2*sampleLocations.size());
            //compute mapsTo and offset for each explicitly stored half block
            for (size_type i = 0; i < sampleLocations.size(); ++i){
                size_type k = sampleLocations[i];
                size_type firstBlock = 2*i, fourthBlock = firstBlock+1;
                //first half block
                size_type start = (k < s_l)? 0 : k-s_l;
                size_type end = (k < slplus1)? 0 : k-slplus1;
//                std::cout << "first half block k " << k << " s_l " << s_l << " slplus1 " << slplus1 << " start "
//                    << start << " end " << end << std::endl;
                if (end > start) {
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                    if (it != sampleLocations.end() && *it < end){
//                        std::cout << "no LF needed, marker found" << std::endl;
                        levels.back().mapsTo[firstBlock] = it - sampleLocations.begin();
                        levels.back().offset[firstBlock] = start + slplus1 - *it;
                        //offset setting above is a possible bug/problem area here and other places it is
                        //set/used (construction/query respectively)
                        usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                        sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                    }
                    else{
                        //never happens at the beginnning or end, so the block is always of size s_l
                        //assert(end-start == s_l);
//                        std::cout << "LF needed, marker not found" << std::endl;
                        gbwt::vector_type T(end-start);
                        gbwt::edge_type position = blockPositionStarts[l][4*i];
                        size_type sequence = this->pathStarts.predecessor(end)->first;
                        this->fullLF(position, sequence);
                        for (size_type u = end; u != start; ){
                            --u;
                            //T is reversed to allow searching using gbwt::GBWT::find
                            T[end-1-u] = position.first;
                            this->fullLF(position, sequence);
                        }
//                        std::cout << "looking for T["<< start << ',' << end << ") = ";
//                            for (auto a : T)
//                            std::cout << a << ' ';
//                        std::cout << "in GBWT" << std::endl;
                        size_type first = 0;
                        gbwt::SearchState range = this->source->rindex->find(T.begin(), T.end(), first);

                        for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                            std::cout << "found in range, seqId " << this->source->rindex->seqId(a) << " offset " << this->source->rindex->seqOffset(a);
                            a = this->FLsuffToTrueSuff(a);
//                            std::cout << " true suff " << a << std::endl; 
                            auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                            if (it != sampleLocations.end() && *it < a + (end-start)){
//                                std::cout << " is in range, *it " << *it << std::flush;
                                levels.back().mapsTo[firstBlock] = it - sampleLocations.begin();
                                levels.back().offset[firstBlock] = a + slplus1 - *it;
//                                std::cout << " mapsto " << levels.back().mapsTo[firstBlock] << " offset " << levels.back().offset[firstBlock] << std::endl;
                                usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                                sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                                break;
                            }
                        }
                    }
                }
                else { /*std::cout << "end <= start" << std::endl;*/ levels.back().mapsTo[firstBlock] = levels.back().offset[firstBlock] = 0; 
                    usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                    sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                }
                //fourth half block
                end = (k > n-s_l)? n : k+s_l;
                start = (k > n-slplus1)? n : k+slplus1;
//                std::cout << "fourth half block k " << k << " s_l " << s_l << " slplus1 " << slplus1 << " start "
//                    << start << " end " << end << std::endl;
                if (end > start) {
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                    if (it != sampleLocations.end() && *it < end){
//                        std::cout << "no LF needed, marker found" << std::endl;
                        levels.back().mapsTo[fourthBlock] = it - sampleLocations.begin();
                        levels.back().offset[fourthBlock] = start + slplus1 - *it;
                        //offset setting above is a possible bug/problem area here and other places it is
                        //set/used (construction/query respectively)
                        usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                        sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                    }
                    else{
                        //never happens at the beginnning or end, so the block is always of size s_l
//                        std::cout << "LF needed, marker not found" << std::endl;
                        gbwt::vector_type T(end-start);
                        gbwt::edge_type position = blockPositionStarts[l][4*i+3];
                        size_type sequence = this->pathStarts.predecessor(end)->first;
                        this->fullLF(position, sequence);
                        for (size_type u = end; u != start; ){
                            --u;
                            //T is reversed to allow searching using gbwt::GBWT::find
                            T[end-1-u] = position.first;
                            this->fullLF(position, sequence);
                        }
//                        std::cout << "looking for T["<< start << ',' << end << ") = ";
//                            for (auto a : T)
//                            std::cout << a << ' ';
//                        std::cout << " in GBWT" << std::endl;
                        size_type first = 0;
                        gbwt::SearchState range = this->source->rindex->find(T.begin(), T.end(), first);

                        for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                            std::cout << "found in range, seqId " << this->source->rindex->seqId(a) << " offset " << this->source->rindex->seqOffset(a);
                            a = this->FLsuffToTrueSuff(a);
//                            std::cout << " true suff " << a << std::endl; 
                            auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                            if (it != sampleLocations.end() && *it < a + (end-start)){
//                                std::cout << " is in range, *it " << *it << std::flush;
                                levels.back().mapsTo[fourthBlock] = it - sampleLocations.begin();
                                levels.back().offset[fourthBlock] = a + slplus1 - *it;
//                                std::cout << " mapsto " << levels.back().mapsTo[fourthBlock] << " offset " << levels.back().offset[fourthBlock] << std::endl;
                                usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                                sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                                break;
                            }
                        }
                    }
                }
                else { levels.back().mapsTo[fourthBlock] = levels.back().offset[fourthBlock] = 0; 
                    usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                    sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                }
            }

            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Out of " << sampleLocations.size() << " sampled locations, " << usedCount << " unique blocks were mapped to by the " << sampleLocations.size()*2 << " first and fourth half blocks from this level of the block tree." << std::endl;
            }
//            std::cout << "levels.back() ";
//            for (size_type i = 0; i < levels.back().mapsTo.size(); ++i)
//                std::cout << '(' << levels.back().mapsTo[i] << ',' << levels.back().offset[i] << ") ";
//            std::cout << std::endl;
//            std::cout << "before compression levels.back() size " << sdsl::size_in_bytes(levels.back()) << std::endl;
            levels.back().compress();
            size_type thisLevelSize = sdsl::size_in_bytes(levels.back());
            treeAndFirstLevelSize += thisLevelSize;
            textSamplesEstimateBits = s_l*sampleLocations.size()*sdsl::bits::length(this->source->rindex->index->effective()); //not keeping track of +64 
            sizeEstimate.push_back(treeAndFirstLevelSize + ceilDiv(textSamplesEstimateBits, 8 ));
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): This level takes " << thisLevelSize << " bytes. If the tree had no more levels, the text samples would take roughly " << ceilDiv(textSamplesEstimateBits, 8) << " bytes." 
                    << " This is a total of " << sizeEstimate.back() << " bytes if the tree was truncated at this level." << std::endl;
                std::cerr << "CompText::CompText(): Construction of this level took " << gbwt::readTimer() - levelStartTime << " seconds" << std::endl;
            }
        }
//        std::cout << "sizeEstimate ";
//        for (auto a: sizeEstimate)
//            std::cout << a << ' ';
//        std::cout << std::endl;
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Construction of all levels took " << gbwt::readTimer() - firstLevelStart << " seconds" << std::endl;
        }
    }
    size_type numLevels;
    {
        //trim tree to smallest size, always has at least one internal level.
        auto minIter = std::min_element(sizeEstimate.begin()+1, sizeEstimate.end());
        numLevels = minIter - sizeEstimate.begin();
        //numLevels = 3; //REMOVE LATER
        
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): The tree will have " << numLevels << " levels since this results in the smallest estimated size (" << sizeEstimate[numLevels] << " bytes)" << std::endl;
        }
        levels.resize(numLevels);
        levels.shrink_to_fit();
    }
    {
        double textSamplesStart = gbwt::readTimer();
        //make text samples
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Making and storing text samples" << std::endl;
        }
        //note, samples of gbwt::GBWT::toComp(node_type)
        size_type s_l = s_0 >> numLevels;
        this->textSamples.resize(sampleLocations.size()*2*s_l);
        for (size_type i = 0; i < sampleLocations.size(); ++i){
            size_type k = sampleLocations[i];
            size_type blockStart = 2*i*s_l;
//            std::cout << "\ntext Sample block " << i << " location " << k << " blockStart " << blockStart << std::flush;
            for (size_type j = 0; j < s_l*2; ++j)
                textSamples[blockStart + j] = 0;
            //first block
            size_type start = (k < s_l)? 0 : k-s_l;
            size_type end = k;
            gbwt::edge_type position = blockPositionStarts[numLevels-1][4*i+1];
            size_type sequence = this->pathStarts.predecessor(end)->first;
//            std::cout << "\nbefore lf position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
            if (end != n)
                this->fullLF(position, sequence);
            else {sequence = l.rindex->index->sequences() - 1; position = {gbwt::ENDMARKER, sequence};}
//            std::cout << "first block start " << start << " end " << end;
            for (size_type u = 0; u != end-start; ++u){
                textSamples[blockStart+s_l-1-u] = l.rindex->index->toComp(position.first);
//                std::cout << " position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
                this->fullLF(position, sequence);
            }
            //second block
            start = k;
            end = (k > n-s_l)? n : k+s_l;
            blockStart += s_l;
            position = blockPositionStarts[numLevels-1][4*i+2];
            sequence = this->pathStarts.predecessor(end)->first;
//            std::cout << "before lf position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
            if (end != n)
                this->fullLF(position, sequence);
            else {sequence = l.rindex->index->sequences() - 1; position = {gbwt::ENDMARKER, sequence};}
//            std::cout << "second block start " << start << " end " << end;
            for (size_type u = 0; u != end-start; ++u){
                textSamples[blockStart+(end-start)-1-u] = l.rindex->index->toComp(position.first);
//                std::cout << " position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
                this->fullLF(position, sequence);
            }
        }
//        std::cout << "textSamples ";
//        for (auto i : textSamples)
//            std::cout << l.rindex->index->toNode(i) << ' ';
//        std::cout << std::endl;
//        std::cout << "before compression textSamples size " << sdsl::size_in_bytes(textSamples) << std::endl;
        sdsl::util::bit_compress(textSamples);
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): The text samples take " << sdsl::size_in_bytes(textSamples) << " bytes of memory." << std::endl;
            std::cerr << "CompText::CompText(): Text samples computed in " << gbwt::readTimer() - textSamplesStart << " seconds." << std::endl;
        }
    }
    if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
    {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "CompText::CompText(): Processed " << this->source->rindex->index->sequences() << " sequences of total length " << this->textLength() << " in " << seconds << " seconds" << std::endl;
        std::cerr << "CompText::CompText(): Final size in bytes: " << sdsl::size_in_bytes(*this) << std::endl;
    }
}

void CompText::buildFullMem(const FastLCP & l) {
    *this = CompText();
    double start = gbwt::readTimer();
    this->source = &l;
    //unpruned construction
    if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    {
        std::cerr << "CompText::CompText(): Unpruned construction of compressed text data structure" << std::endl;
    }
    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Computing path starts" << std::endl;
        }
        //compute path starts 
        std::vector<size_type> pathLengthPrefixSums;
        size_type prefixSum = 0, first = l.rindex->locateFirst(gbwt::ENDMARKER);
        for (size_type i = 0; i < l.rindex->index->sequences()-1; ++i){
            pathLengthPrefixSums.push_back(prefixSum);
            prefixSum += l.rindex->seqOffset(first) + 1; // + 1 for endmarker character
            first = l.rindex->locateNext(first);
        }
        pathLengthPrefixSums.push_back(prefixSum);
        prefixSum += l.rindex->seqOffset(first) + 1; // + 1 for endmarker character
        sdsl::sd_vector_builder builder(prefixSum, pathLengthPrefixSums.size());
        for (const auto & a : pathLengthPrefixSums)
            builder.set(a);
        this->pathStarts = sdsl::sd_vector<>(builder);
    }
    std::vector<size_type> sampleLocations;
    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Computing sample locations" << std::endl;
        }
        //compute sorted sampleLocations by true suffix values
        std::vector<size_type> temp = {0, this->pathStarts.size()-1};

        //add samples for beginning and end of concrete runs
        for (const auto & a : this->source->rindex->samples)
            temp.push_back(this->FLsuffToTrueSuff(a));
        for (const auto & a : this->source->samples_bot)
            temp.push_back(this->FLsuffToTrueSuff(a));

        //possible bug if samples and samples_bot don't include samples for logical runs,only concrete runs
        //i.e. if bwt[i] is endmarker, samples need to be stored
        //if missing, add samples for endmarkers in GBWT below

        gbwt::parallelQuickSort(temp.begin(), temp.end());
        for (const auto & a : temp)
            if (sampleLocations.size() == 0 || a != sampleLocations.back())
                sampleLocations.push_back(a);
//        std::cout << "sampleLocations ";
//           for (auto a:  sampleLocations) 
//              std::cout << a << ' ';
//        std::cout << std::endl;
    }
    //this->temp_sl = sampleLocations;


    //sizeEstimate[i] is an estimate for the size of the whole data structure if it had i intermediate levels
    //(if levels.size() was i)
    std::vector<size_type> sizeEstimate;
    
    auto ceilDiv = [] (size_type a, size_type b) { return (a/b) + (a%b != 0); };
    size_type r = l.samples_bot.size(), n = this->pathStarts.size();
    if (n == r) { throw std::invalid_argument("Number of characters and runs in text are the same, compressing won't work"); }
    this->s_0 = 1 << sdsl::bits::hi(ceilDiv(n,r));
    if (this->s_0 < ceilDiv(n,r)) { this->s_0 *= 2; }
    if (gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    { 
        std::cerr << "CompText::CompText(): The text of the GBWT has length " << n << " and " << r << " logical runs. Average logical run length is " << double(n)/r 
            << ". The smallest power of 2 larger than or equal to the average run length is " << this->s_0 << ". This is the width of the blocks in the first level of the compressed text data structure." << std::endl;
    }

    size_type maxlevels = sdsl::bits::hi(s_0);
    if (gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    {
        std::cerr << "CompText::CompText(): The maximum number of levels in the tree is " << maxlevels << ". This is the log_2 of the first level's block width." << std::endl;
        std::cerr << "CompText::CompText(): Extracting full text" << std::endl;
    }

    double textExtractStart = gbwt::readTimer();
    gbwt::vector_type fullText(n, gbwt::ENDMARKER);
    {
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_type i = 0; i < this->pathStarts.ones(); ++i){
            //extract sequence i
            double seqStart = gbwt::readTimer();
            size_type end = this->pathStarts.select_iter(i+2)->second - 1;
            gbwt::edge_type curr = {gbwt::ENDMARKER, i};
            curr = this->source->rindex->index->LF(curr);
            while (curr.first != gbwt::ENDMARKER) 
            {
                fullText[--end] = curr.first;
                curr = this->source->rindex->index->LF(curr);
            }
            #pragma omp critical
            {
                std::cout << "Extracted seq " << i << " in " << gbwt::readTimer() - seqStart << " seconds" << std::endl;
            }
        }
//        std::cout << "Extracted fullText: " << std::flush;
//        for (size_type i = 0; i<fullText.size(); ++i)
//            std::cout << fullText[i] << ',' << std::flush;
//        std::cout << std::endl;
    }
    if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
    {
        std::cerr << "CompText::CompText(): Extracted full text of length " << n << " with " << this->pathStarts.ones() << " sequences in " << gbwt::readTimer() - textExtractStart << " seconds"  << std::endl;
    }

    {
        if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Building block tree" << std::endl;
        }
        //build block tree 

        //auto ithStart = [] (size_type i, size_type n, size_type r) {
        //    //assumes 0 <= i <= r <= n
        //    return i*(n/r) + ((i*(n%r))/r);
        //};

        double firstLevelStart = gbwt::readTimer();
        {
            if(gbwt::Verbosity::level >= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Building the first level of the block tree" << std::endl;
            }
            //build first level
            size_type topLevelRuns = ceilDiv(n, this->s_0);
            this->firstLevel.resize(topLevelRuns);

            std::vector<bool> sampleUsed(sampleLocations.size(), false);
            size_type usedCount = 0;
            size_type runStart = gbwt::readTimer();
            auto ithStart = [n, this] (size_type i) { return std::min(this->s_0*i, n); };
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_type i = 1; i <= topLevelRuns; ++i){
                size_type blockInd = i-1;
                size_type start = ithStart(blockInd);
                size_type end = ithStart(blockInd+1);
                
                if (blockInd%10000 == 0) {
                    std::cout << "On run " << blockInd << " after " << gbwt::readTimer() - runStart << " seconds" << std::endl;
                    runStart = gbwt::readTimer();
                }
//                std::cout << "start " << start << " end " << end << std::endl << "sampleLocations.size() " << sampleLocations.size() << std::endl;
//                std::cout << "run " << blockInd << std::endl;
//                std::cout << "LF done " << std::endl;
                //assign block blockInd to a run, block blockInd corresponds to [blockInd*s_0, (blockInd+1)*s_0) = [start, end)
                //where s_0 is the smallest power of 2 larger than or equal to \ceiling(n/r) (and is the block size of the first level of the tree
                //where n is the total length and r is the number of GBWT logical runs
                auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                if (it != sampleLocations.end() && *it < end){
//                    std::cout << "no LF needed, marker found" << std::endl;
//                    std::cout << "it - sampleLocations.begin() " << it - sampleLocations.begin() << " *it " << *it << " start " << start << std::endl;
//                    std::cout << " blockInd " << blockInd << std::endl;
                    firstLevel.mapsTo[blockInd] = it - sampleLocations.begin();
                    firstLevel.offset[blockInd] = start + this->s_0 - *it;
//                    std::cout << " mapsto " << firstLevel.mapsTo[blockInd] << " offset " << firstLevel.offset[blockInd] << std::endl;
                    usedCount += !sampleUsed[firstLevel.mapsTo[blockInd]];
                    sampleUsed[firstLevel.mapsTo[blockInd]] = true;
//                    std::cout << " mapsto " << firstLevel.mapsTo[blockInd] << " offset " << firstLevel.offset[blockInd] << std::endl;
                    continue;
                }
                //all endmarkers  are sampled therefore [start, end) contains no endmarkers
                //if [start, end) is not a primary occurrence of T[start, end).
                //Here we search for a primary occurence

//                std::cout << "LF needed, marker not found" << std::endl;
                    
                //assert(this->s_0 == end-start);
                size_type first = 0;
//                std::cout << "looking for T["<< start << ',' << end << ") = ";
//                    for (auto a : T)
//                       std::cout << a << ' ';
//                std::cout  << " in GBWT" << std::endl;
                gbwt::SearchState range = l.rindex->find(fullText.rend() - end, fullText.rend() - start, first);

                for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                    std::cout << "found in range, seqId " << l.rindex->seqId(a) << " offset " << l.rindex->seqOffset(a);
                    a = this->FLsuffToTrueSuff(a);
//                    std::cout << " true suff " << a << std::endl; 
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                    if (it != sampleLocations.end() && *it < a + (end-start)){
//                        std::cout << " is in range, *it " << *it;
                        firstLevel.mapsTo[blockInd] = it - sampleLocations.begin();
                        firstLevel.offset[blockInd] = a + this->s_0 - *it;
//                        std::cout << " mapsto " << firstLevel.mapsTo[i] << " offset " << firstLevel.offset[i] << std::endl;
                        usedCount += !sampleUsed[firstLevel.mapsTo[blockInd]];
                        sampleUsed[firstLevel.mapsTo[blockInd]] = true;
                        break;
                    }
                }
            }
            //need to compress mapsto to rank using sampleUnused here. Do once compressing block tree
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Out of " << sampleLocations.size() << " sampled locations, " << usedCount << " unique blocks were mapped to by the " << topLevelRuns << " blocks in the first level of the block tree." << std::endl;
            }
        }
//        std::cout << "firstLevel ";
//        for (size_type i = 0; i < firstLevel.mapsTo.size(); ++i)
//            std::cout << '(' << firstLevel.mapsTo[i] << ',' << firstLevel.offset[i] << ") ";
//        std::cout << std::endl;
//        std::cout << "before compression firstLevel size " << std::flush << sdsl::size_in_bytes(firstLevel) << std::endl;
        firstLevel.compress();

        //build the rest of the levels
        
        size_type treeAndFirstLevelSize = sdsl::size_in_bytes(this->pathStarts) + 
                sdsl::size_in_bytes(firstLevel);
        size_type textSamplesEstimateBits = 2*s_0*firstLevel.mapsTo.size()*sdsl::bits::length(this->source->rindex->index->effective()); //+64?
        sizeEstimate.push_back(treeAndFirstLevelSize + ceilDiv(textSamplesEstimateBits, 8));
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): First level takes " << treeAndFirstLevelSize << " bytes. If the tree had no more levels, the text samples would take roughly " << ceilDiv(textSamplesEstimateBits, 8) << " bytes." 
                << " This is a total of " << sizeEstimate[0] << " bytes if the tree was truncated at this level." << std::endl;
            std::cerr << "CompText::CompText(): Construction of the first level took " << gbwt::readTimer() - firstLevelStart << " seconds" << std::endl;
        }
        std::vector<bool> sampleUsed(sampleLocations.size());
        size_type usedCount = 0;
        for(size_type l = 0, s_l = s_0; s_l != 1; ++l, s_l /= 2){
            double levelStartTime = gbwt::readTimer();
            sampleUsed.assign(sampleLocations.size(), false);
            usedCount = 0;
            //compute intermediate level l (levels[l])
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Computing level " << l << " of the tree. Blocks at this level have a width of " << s_l << std::endl;
            }
            size_type slplus1 = s_l/2;

            //since this is the unpruned construction, half blocks adjacent to samples are
            //not stored explicitly, they are implicitly mapped to the same sample in the below level
            //therefore each level has 2*samples half blocks
            this->levels.push_back(Level());
            this->levels.back().resize(2*sampleLocations.size());
            //compute mapsTo and offset for each explicitly stored half block
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_type i = 0; i < sampleLocations.size(); ++i){
                size_type k = sampleLocations[i];
                size_type firstBlock = 2*i, fourthBlock = firstBlock+1;
                //first half block
                size_type start = (k < s_l)? 0 : k-s_l;
                size_type end = (k < slplus1)? 0 : k-slplus1;
//                std::cout << "first half block k " << k << " s_l " << s_l << " slplus1 " << slplus1 << " start "
//                    << start << " end " << end << std::endl;
                if (end > start) {
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                    if (it != sampleLocations.end() && *it < end){
//                        std::cout << "no LF needed, marker found" << std::endl;
                        levels.back().mapsTo[firstBlock] = it - sampleLocations.begin();
                        levels.back().offset[firstBlock] = start + slplus1 - *it;
                        //offset setting above is a possible bug/problem area here and other places it is
                        //set/used (construction/query respectively)
                        usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                        sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                    }
                    else{
                        //never happens at the beginnning or end, so the block is always of size s_l
                        //assert(end-start == s_l);
//                        std::cout << "LF needed, marker not found" << std::endl;
//                        std::cout << "looking for T["<< start << ',' << end << ") = ";
//                            for (auto a : T)
//                            std::cout << a << ' ';
//                        std::cout << "in GBWT" << std::endl;
                        size_type first = 0;
                        gbwt::SearchState range = this->source->rindex->find(fullText.rend() - end, fullText.rend() - start, first);

                        for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                            std::cout << "found in range, seqId " << this->source->rindex->seqId(a) << " offset " << this->source->rindex->seqOffset(a);
                            a = this->FLsuffToTrueSuff(a);
//                            std::cout << " true suff " << a << std::endl; 
                            auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                            if (it != sampleLocations.end() && *it < a + (end-start)){
//                                std::cout << " is in range, *it " << *it << std::flush;
                                levels.back().mapsTo[firstBlock] = it - sampleLocations.begin();
                                levels.back().offset[firstBlock] = a + slplus1 - *it;
//                                std::cout << " mapsto " << levels.back().mapsTo[firstBlock] << " offset " << levels.back().offset[firstBlock] << std::endl;
                                usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                                sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                                break;
                            }
                        }
                    }
                }
                else { /*std::cout << "end <= start" << std::endl;*/ levels.back().mapsTo[firstBlock] = levels.back().offset[firstBlock] = 0; 
                    usedCount += !sampleUsed[levels.back().mapsTo[firstBlock]];
                    sampleUsed[levels.back().mapsTo[firstBlock]] = true;
                }
                //fourth half block
                end = (k > n-s_l)? n : k+s_l;
                start = (k > n-slplus1)? n : k+slplus1;
//                std::cout << "fourth half block k " << k << " s_l " << s_l << " slplus1 " << slplus1 << " start "
//                    << start << " end " << end << std::endl;
                if (end > start) {
                    auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), start);
                    if (it != sampleLocations.end() && *it < end){
//                        std::cout << "no LF needed, marker found" << std::endl;
                        levels.back().mapsTo[fourthBlock] = it - sampleLocations.begin();
                        levels.back().offset[fourthBlock] = start + slplus1 - *it;
                        //offset setting above is a possible bug/problem area here and other places it is
                        //set/used (construction/query respectively)
                        usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                        sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                    }
                    else{
                        //never happens at the beginnning or end, so the block is always of size s_l
//                        std::cout << "LF needed, marker not found" << std::endl;
//                        std::cout << "looking for T["<< start << ',' << end << ") = ";
//                            for (auto a : T)
//                            std::cout << a << ' ';
//                        std::cout << " in GBWT" << std::endl;
                        size_type first = 0;
                        gbwt::SearchState range = this->source->rindex->find(fullText.rend() - end, fullText.rend() - start, first);

                        for (size_type j = 0, a = first; j < range.size(); ++j, first = this->source->rindex->locateNext(first)){
//                            std::cout << "found in range, seqId " << this->source->rindex->seqId(a) << " offset " << this->source->rindex->seqOffset(a);
                            a = this->FLsuffToTrueSuff(a);
//                            std::cout << " true suff " << a << std::endl; 
                            auto it = std::lower_bound(sampleLocations.begin(), sampleLocations.end(), a);
                            if (it != sampleLocations.end() && *it < a + (end-start)){
//                                std::cout << " is in range, *it " << *it << std::flush;
                                levels.back().mapsTo[fourthBlock] = it - sampleLocations.begin();
                                levels.back().offset[fourthBlock] = a + slplus1 - *it;
//                                std::cout << " mapsto " << levels.back().mapsTo[fourthBlock] << " offset " << levels.back().offset[fourthBlock] << std::endl;
                                usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                                sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                                break;
                            }
                        }
                    }
                }
                else { levels.back().mapsTo[fourthBlock] = levels.back().offset[fourthBlock] = 0; 
                    usedCount += !sampleUsed[levels.back().mapsTo[fourthBlock]];
                    sampleUsed[levels.back().mapsTo[fourthBlock]] = true;
                }
            }

            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): Out of " << sampleLocations.size() << " sampled locations, " << usedCount << " unique blocks were mapped to by the " << sampleLocations.size()*2 << " first and fourth half blocks from this level of the block tree." << std::endl;
            }
//            std::cout << "levels.back() ";
//            for (size_type i = 0; i < levels.back().mapsTo.size(); ++i)
//                std::cout << '(' << levels.back().mapsTo[i] << ',' << levels.back().offset[i] << ") ";
//            std::cout << std::endl;
//            std::cout << "before compression levels.back() size " << sdsl::size_in_bytes(levels.back()) << std::endl;
            levels.back().compress();
            size_type thisLevelSize = sdsl::size_in_bytes(levels.back());
            treeAndFirstLevelSize += thisLevelSize;
            textSamplesEstimateBits = s_l*sampleLocations.size()*sdsl::bits::length(this->source->rindex->index->effective()); //not keeping track of +64 
            sizeEstimate.push_back(treeAndFirstLevelSize + ceilDiv(textSamplesEstimateBits, 8 ));
            if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
            {
                std::cerr << "CompText::CompText(): This level takes " << thisLevelSize << " bytes. If the tree had no more levels, the text samples would take roughly " << ceilDiv(textSamplesEstimateBits, 8) << " bytes." 
                    << " This is a total of " << sizeEstimate.back() << " bytes if the tree was truncated at this level." << std::endl;
                std::cerr << "CompText::CompText(): Construction of this level took " << gbwt::readTimer() - levelStartTime << " seconds" << std::endl;
            }
        }
//        std::cout << "sizeEstimate ";
//        for (auto a: sizeEstimate)
//            std::cout << a << ' ';
//        std::cout << std::endl;
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Construction of all levels took " << gbwt::readTimer() - firstLevelStart << " seconds" << std::endl;
        }
    }
    size_type numLevels;
    {
        //trim tree to smallest size
        auto minIter = std::min_element(sizeEstimate.begin(), sizeEstimate.end());
        numLevels = minIter - sizeEstimate.begin();
        //numLevels = 3; //REMOVE LATER
        
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): The tree will have " << numLevels << " levels since this results in the smallest estimated size (" << sizeEstimate[numLevels] << " bytes)" << std::endl;
        }
        levels.resize(numLevels);
        levels.shrink_to_fit();
    }
    {
        double textSamplesStart = gbwt::readTimer();
        //make text samples
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): Making and storing text samples" << std::endl;
        }
        //note, samples of gbwt::GBWT::toComp(node_type)
        size_type s_l = s_0 >> numLevels;
        this->textSamples.resize(sampleLocations.size()*2*s_l);
        for (size_type i = 0; i < sampleLocations.size(); ++i){
            size_type k = sampleLocations[i];
            size_type blockStart = 2*i*s_l;
//            std::cout << "\ntext Sample block " << i << " location " << k << " blockStart " << blockStart;
            for (size_type j = 0; j < s_l*2; ++j)
                textSamples[blockStart + j] = 0;
            //first block
            size_type start = (k < s_l)? 0 : k-s_l;
            size_type end = k;
            auto it = fullText.begin() + start;
//            std::cout << "\nbefore lf position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
//            std::cout << "first block start " << start << " end " << end;
            for (size_type u = start; u != end; ++it, ++u){
                textSamples[blockStart+s_l-(end-u)] = l.rindex->index->toComp(*it);
//                std::cout << " position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
            }
            //second block
            start = k;
            end = (k > n-s_l)? n : k+s_l;
            blockStart += s_l;
            it = fullText.begin() + k;
            for (size_type u = start; u != end; ++it, ++u){
                textSamples[blockStart+(u-start)] = l.rindex->index->toComp(*it);
//                std::cout << " position " << '(' << position.first << ',' << position.second  << ") sequence " << sequence << std::endl;
            }
        }
//        std::cout << "textSamples ";
//        for (auto i : textSamples)
//            std::cout << l.rindex->index->toNode(i) << ' ';
//        std::cout << std::endl;
//        std::cout << "before compression textSamples size " << sdsl:size_in_bytes(textSamples) << std::endl;
        sdsl::util::bit_compress(textSamples);
        if (gbwt::Verbosity::level>= gbwt::Verbosity::FULL)
        {
            std::cerr << "CompText::CompText(): The text samples take " << sdsl::size_in_bytes(textSamples) << " bytes of memory." << std::endl;
            std::cerr << "CompText::CompText(): Text samples computed in " << gbwt::readTimer() - textSamplesStart << " seconds." << std::endl;
        }
    }
    if(gbwt::Verbosity::level >= gbwt::Verbosity::BASIC)
    {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "CompText::CompText(): Processed " << this->source->rindex->index->sequences() << " sequences of total length " << this->textLength() << " in " << seconds << " seconds" << std::endl;
        std::cerr << "CompText::CompText(): Final size in bytes: " << sdsl::size_in_bytes(*this) << std::endl;
    }
}

gbwt::node_type CompText::at(size_type i) const{
//    std::cout << "in CompText::at(" << i << ")" << std::endl;
    if (i >= this->textLength())
        return gbwt::invalid_offset();
    

    size_type blockIndex, offset;
    {
        //firstLevel
        blockIndex = i/this->s_0;
        offset = i%this->s_0;
        //levels[0]
        offset += firstLevel.offset[blockIndex];
        //assert(offset < 2*this->s_0);
        blockIndex = 2*firstLevel.mapsTo[blockIndex] + (offset/this->s_0);
        offset = offset%this->s_0;
    }


    //within blockIndex
    size_type level = 0;
    size_type s_l = this->s_0;
    size_type slplus1 = s_l/2;
//    std::cout << "blockIndex " << blockIndex << " offset " << offset << " s_0 " << this->s_0 << " s_l " << s_l << " slplus1 " << slplus1 << std::endl;
//    std::cout << "firstLevel:\n";

    for (;level < levels.size(); ++level, s_l = slplus1, slplus1/=2){
        //currently at the blockIndex-th block of levels[level]
        //at an offset of offset within the block
//        std::cout << "level " << level << " s_l " << s_l << " slplus1 " << slplus1 << " blockIndex " << blockIndex << " offset " << offset << std::endl;
        if ((blockIndex%2 == 0 && offset >= slplus1) || (blockIndex%2 == 1 && offset < slplus1)){
            if (blockIndex%2 == 0)
                offset -= slplus1;
            continue; //the second and third half blocks implicitly move to the next level
        }

        //this line needs to change if explicitly storing all four half blocks
        if (blockIndex%2 == 1)
            offset -= slplus1;
        offset += levels[level].offset[blockIndex];
//        std::cout << "offset before mod " << offset << std::endl;
        blockIndex = 2*levels[level].mapsTo[blockIndex]+(offset/slplus1);
        offset = offset % slplus1;
    }
//    std::cout << "level " << level << " s_l " << s_l << " slplus1 " << slplus1 << " blockIndex " << blockIndex << " offset " << offset << std::endl;
//    std::cout << "return " << this->source->rindex->index->toNode(this->textSamples[blockIndex*s_l + offset]) << std::endl;
    return this->source->rindex->index->toNode(this->textSamples[blockIndex*s_l + offset]);
}


#endif //GBWT_QUERY_COMPTEXT_H
