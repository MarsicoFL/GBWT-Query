#ifndef GBWT_QUERY_FAST_LCP_H
#define GBWT_QUERY_FAST_LCP_H

//heavily derived/modified from FastLocate

#include<gbwt/gbwt.h>
#include<gbwt/fast_locate.h>
#include<lf_gbwt.h>

struct FastLCP{
    typedef gbwt::size_type size_type;

    const gbwt::FastLocate* rindex;

    sdsl::int_vector<0> samples_bot;
    sdsl::sd_vector<> first;
    sdsl::int_vector<0> first_to_run;

    sdsl::int_vector<0> samples_lcp;

    FastLCP() = default;
    explicit FastLCP(const gbwt::FastLocate& source);
    explicit FastLCP(const lf_gbwt::GBWT&, const gbwt::FastLocate&);

    size_type serialize(std::ostream &out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
    void load(std::istream& in);

    const static std::string EXTENSION; //.flcp

    size_type locatePrev(size_type next) const;
    size_type LCP(size_type suff) const;

    size_type globalRunId(gbwt::node_type node, size_type run_id) const
    {
        return this->rindex->comp_to_run[this->rindex->index->toComp(node)] + run_id;
    }

    size_type getSample(gbwt::node_type node, size_type run_id) const
    {
        return this->rindex->samples[this->globalRunId(node, run_id)];
    }

    size_type getSampleBot(gbwt::node_type node, size_type run_id) const
    {
        return this->samples_bot[this->globalRunId(node, run_id)];
    }
    
    size_type locateLast(gbwt::node_type node) const
    {
        size_type runId = (this->rindex->index->toComp(node)+1 == this->rindex->index->effective())?
            samples_bot.size()-1:
            this->rindex->comp_to_run[this->rindex->index->toComp(node)+1] - 1;
        //std::cout << "samples_bot.size " << this->samples_bot.size()
            //<< " samples.size " << this->rindex->samples.size()
            //<< " runId " << runId
            //<< " node " << node
            //<< " this->globalRunId(node+1, 0) - 1 " << this->globalRunId(node+1, 0) - 1 
            //<< " samples_bot.size()-1: " << samples_bot.size()-1
            //<< std::endl;
        return this->samples_bot[runId];
    }

    //returns path length not including endmarker
    size_type pathLength(gbwt::size_type i) {
        if (i >= this->rindex->seqId(this->first.size()))
            return gbwt::invalid_offset();
        auto iter = this->first.predecessor(this->rindex->pack(i+1,0)-1);
        return iter->second + this->samples_lcp[iter->first];
    }

    std::vector<gbwt::size_type> decompressSA(const gbwt::node_type node) {
        if (!this->rindex->index->contains(node)) { return {}; }
        gbwt::size_type n = this->rindex->index->nodeSize(node);
        if (n == 0) { return {}; }
        std::vector<gbwt::size_type> rec(n);
        rec.shrink_to_fit();
        rec.back() = this->locateLast(node);
        for (gbwt::size_type i = n - 1; i != 0; --i) { rec[i-1] = this->locatePrev(rec[i]); }
        return rec;
    }

    bool verifySuff() {
        bool good = true;
        #pragma omp parallel for schedule(dynamic, 1)
        for (gbwt::comp_type i = 0; i < this->rindex->index->effective(); ++i) {
            gbwt::size_type node = this->rindex->index->effective();
            bool thisNodeGood = (this->decompressSA(node) == this->rindex->decompressSA(node));
            #pragma omp critical 
            {
                good = good && thisNodeGood;
            }
        }
        return good;
    }
};

FastLCP::size_type FastLCP::locatePrev(size_type next) const{
    auto iter = this->first.predecessor(next);
    return this->samples_bot[this->first_to_run[iter->first]-1] + (next - iter->second); 

}

FastLCP::size_type FastLCP::LCP(size_type suff) const{
    auto iter = this->first.predecessor(suff);
    //return this->samples_lcp[this->first_to_run[iter->first]] - (suff - iter->second);
    return this->samples_lcp[iter->first] - (suff - iter->second);
}


FastLCP::FastLCP(const gbwt::FastLocate& source) : rindex(&source) {
    using namespace gbwt;
    double start = readTimer();
    const char msgPrefix[] = "FastLCP::FastLCP(FastLocate): ";

    if(this->rindex->empty())
    {
        if(Verbosity::level >= Verbosity::FULL)
        {
            std::cerr << msgPrefix << "The input FastLocate data structure is empty" << std::endl;
        }
        return;
    }

    // Determine the number of logical runs before each record.
    size_type total_runs = this->rindex->samples.size();
    //this->rindex->index->bwt.forEach([&](size_type comp, const CompressedRecord& record)
            //{
            //total_runs += record.runs().second;
            //});
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << total_runs << " logical runs in the GBWT" << std::endl;
    }

    // Global sample buffers.
    struct sample_record
    {
        size_type seq_id, seq_offset, run_id;

        // Sort by text position.
        bool operator<(const sample_record& another) const
        {
            return (this->seq_id < another.seq_id || (this->seq_id == another.seq_id && this->seq_offset < another.seq_offset));
        }
    };
    std::vector<sample_record> head_samples, tail_samples;
    head_samples.reserve(total_runs);
    tail_samples.reserve(total_runs);

    // Run identifier for each offset in the endmarker. We cannot get this
    // information efficiently with random access.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Processing the endmarker record" << std::endl;
    }
    std::vector<size_type> endmarker_runs(this->rindex->index->sequences(), 0);
    {
        size_type run_id = 0;
        edge_type prev = this->rindex->index->start(0);
        for(size_type i = 1; i < this->rindex->index->sequences(); i++)
        {
            edge_type curr = this->rindex->index->start(i);
            if(curr.first == ENDMARKER || curr.first != prev.first) { run_id++; prev = curr; }
            endmarker_runs[i] = run_id;
        }
    }

    // Extract the samples from each sequence.
    double extract_start = readTimer();
   if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Extracting head/tail samples" << std::endl;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for(size_type i = 0; i < this->rindex->index->sequences(); i++)
    {
        size_type start = gbwt::readTimer();
        std::vector<sample_record> head_buffer, tail_buffer;
        size_type seq_offset = 0, run_id = endmarker_runs[i];
        if(i == 0 || run_id != endmarker_runs[i - 1])
        {
            head_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
        }
        if(i + 1 >= this->rindex->index->sequences() || run_id != endmarker_runs[i + 1])
        {
            tail_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
        }
        edge_type curr = this->rindex->index->start(i); seq_offset++;
        range_type run(0, 0);
        while(curr.first != ENDMARKER)
        {
            edge_type next = this->rindex->index->record(curr.first).LF(curr.second, run, run_id);
            if(curr.second == run.first)
            {
                head_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
            }
            if(curr.second == run.second)
            {
                tail_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
            }
            curr = next; seq_offset++;
        }
        // GBWT is an FM-index of the reverse paths. The sequence offset r-index needs
        // is the distance to the BWT position with the endmarker (to the end of the
        // path, to the start of the string).
        for(sample_record& record : head_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
        for(sample_record& record : tail_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
        #pragma omp critical
        {
            head_samples.insert(head_samples.end(), head_buffer.begin(), head_buffer.end());
            tail_samples.insert(tail_samples.end(), tail_buffer.begin(), tail_buffer.end());
            std::cout << "Extracted seq " << i << " in " << gbwt::readTimer() - start << " seconds" << std::endl;
        }
    }
    sdsl::util::clear(endmarker_runs);
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - extract_start;
        std::cerr << msgPrefix << "Extracted " << head_samples.size() << " / " << tail_samples.size() << " head/tail samples in " << seconds << " seconds" << std::endl;
    }

    // Store the tail samples.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the tail samples" << std::endl;
    }
    parallelQuickSort(tail_samples.begin(), tail_samples.end(), [](const sample_record& a, const sample_record& b)
            {
            return (a.run_id < b.run_id);
            });
    this->samples_bot.width(sdsl::bits::length(this->rindex->pack(this->rindex->index->sequences() - 1, this->rindex->header.max_length - 1)));
    this->samples_bot.resize(total_runs);
    for(size_type i = 0; i < total_runs; i++)
    {
        this->samples_bot[i] = this->rindex->pack(tail_samples[i].seq_id, tail_samples[i].seq_offset);
    }
    sdsl::util::clear(tail_samples);

    // Store the head samples.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the head samples" << std::endl;
    }
    parallelQuickSort(head_samples.begin(), head_samples.end());
    sdsl::sd_vector_builder builder(this->rindex->index->sequences() * this->rindex->header.max_length, total_runs);
    this->first_to_run.width(sdsl::bits::length(total_runs - 1));
    this->first_to_run.resize(total_runs);
    for(size_type i = 0; i < total_runs; i++)
    {
        builder.set_unsafe(this->rindex->pack(head_samples[i].seq_id, head_samples[i].seq_offset));
        this->first_to_run[i] = head_samples[i].run_id;
    }
    this->first = sdsl::sd_vector<>(builder);


    //lcp_sample[i] is the sample for head_sample[i] (after head_samples is sorted by the sample_record < operator)
    std::vector<size_type> lcp_samples;
    lcp_samples.resize(total_runs);

    // Compute LCP samples.
    double compute_start = readTimer();
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Computing the LCP samples" << std::endl;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for(size_type i = 0; i < this->rindex->index->sequences(); i++)
    {
        size_type start = gbwt::readTimer();
        sample_record srec = {i+1, 0, 0};
        size_type sampIndex = std::lower_bound(head_samples.begin(), head_samples.end(), srec) - head_samples.begin();
        size_type run_id;
        edge_type curr = {ENDMARKER, i};
        range_type run(0,0);
        do {
            edge_type next = this->rindex->index->record(curr.first).LF(curr.second, run, run_id);
            if (curr.second == run.first){
                //compute lcp
                size_type lcp = 0;
                edge_type above = {curr.first, curr.second-1};
                while (curr.second != 0 && curr.first != ENDMARKER && above.first == curr.first){
                    ++lcp;
                    above = this->rindex->index->inverseLF(above);
                    curr = this->rindex->index->inverseLF(curr);
                }
                lcp_samples[--sampIndex] = lcp;
            }
            curr = next;
        } while (curr.first != ENDMARKER);
        #pragma omp critical
        {
        std::cout << "Extracted seq " << i << " in " << gbwt::readTimer() - start << " seconds" << std::endl;
        }
    }
    
    // Store the LCP samples
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the lcp samples" << std::endl;
    }

    this->samples_lcp.resize(total_runs);
    for (size_type i = 0; i<total_runs; ++i)
        this->samples_lcp[i] = lcp_samples[i];
    sdsl::util::bit_compress(this->samples_lcp);

    sdsl::util::clear(head_samples);
    sdsl::util::clear(lcp_samples);
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - compute_start;
        std::cerr << msgPrefix << "Extracted " << this->samples_lcp.size() << " lcp samples in " << seconds << " seconds" << std::endl;
    }

    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        std::cerr << msgPrefix << "Processed " << this->rindex->index->sequences() << " sequences of total length " << this->rindex->index->size() << " in " << seconds << " seconds" << std::endl;
    }
}

FastLCP::FastLCP(const lf_gbwt::GBWT& lfg, const gbwt::FastLocate& r) : rindex(&r) {
    using namespace gbwt;
    double start = readTimer();
    const char msgPrefix[] = "FastLCP::FastLCP(lf gbwt, FastLocate): ";
    
    if(this->rindex->empty() || lfg.empty())
    {
       if (Verbosity::level >= Verbosity::FULL)
       {
           if (this->rindex->empty())
               std::cerr << msgPrefix << "The input FastLocate data structure is empty" << std::endl;
           if (lfg.empty())
               std::cerr << msgPrefix << "The input lf gbwt data structure is empty" << std::endl;
       }
       return;
    }

    // Determine the number of logical runs before each record.
    size_type total_runs = this->rindex->samples.size();
    //this->rindex->index->bwt.forEach([&](size_type comp, const CompressedRecord& record)
            //{
            //total_runs += record.runs().second;
            //});
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << total_runs << " logical runs in the GBWT" << std::endl;
    }

    // Global sample buffers.
    struct sample_record
    {
        size_type seq_id, seq_offset, run_id;

        // Sort by text position.
        bool operator<(const sample_record& another) const
        {
            return (this->seq_id < another.seq_id || (this->seq_id == another.seq_id && this->seq_offset < another.seq_offset));
        }
    };
    std::vector<sample_record> head_samples, tail_samples;
    head_samples.reserve(total_runs);
    tail_samples.reserve(total_runs);

    // Run identifier for each offset in the endmarker. We cannot get this
    // information efficiently with random access.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Processing the endmarker record" << std::endl;
    }
    std::vector<size_type> endmarker_runs(lfg.sequences(), 0);
    {
        size_type run_id = 0;
        edge_type prev = lfg.LF({ENDMARKER, 0});
        for(size_type i = 1; i < lfg.sequences(); i++)
        {
            edge_type curr = lfg.LF({ENDMARKER, i});
            if(curr.first == ENDMARKER || curr.first != prev.first) { run_id++; prev = curr; }
            endmarker_runs[i] = run_id;
        }
    }

    // Extract the samples from each sequence.
    double extract_start = readTimer();
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Extracting head/tail samples" << std::endl;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for(size_type i = 0; i < lfg.sequences(); i++)
    {
        size_type start = gbwt::readTimer();
        std::vector<sample_record> head_buffer, tail_buffer;
        size_type seq_offset = 0, run_id = endmarker_runs[i];
        if(i == 0 || run_id != endmarker_runs[i - 1])
        {
            head_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
        }
        if(i + 1 >= lfg.sequences() || run_id != endmarker_runs[i + 1])
        {
            tail_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
        }
        edge_type curr = lfg.LF({ENDMARKER, i}); seq_offset++;
        range_type run(0, 0);
        while(curr.first != ENDMARKER)
        {
            edge_type next = lfg.LF(curr, run, run_id);
            if(curr.second == run.first)
            {
                head_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
            }
            if(curr.second == run.second)
            {
                tail_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
            }
            curr = next; seq_offset++;
        }
        // GBWT is an FM-index of the reverse paths. The sequence offset r-index needs
        // is the distance to the BWT position with the endmarker (to the end of the
        // path, to the start of the string).
        for(sample_record& record : head_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
        for(sample_record& record : tail_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
        #pragma omp critical
        {
            head_samples.insert(head_samples.end(), head_buffer.begin(), head_buffer.end());
            tail_samples.insert(tail_samples.end(), tail_buffer.begin(), tail_buffer.end());
            std::cout << "Extracted seq " << i << " in " << gbwt::readTimer() - start << " seconds" << std::endl;
        }
    }
    sdsl::util::clear(endmarker_runs);
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - extract_start;
        std::cerr << msgPrefix << "Extracted " << head_samples.size() << " / " << tail_samples.size() << " head/tail samples in " << seconds << " seconds" << std::endl;
    }

    // Store the tail samples.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the tail samples" << std::endl;
    }
    parallelQuickSort(tail_samples.begin(), tail_samples.end(), [](const sample_record& a, const sample_record& b)
            {
            return (a.run_id < b.run_id);
            });
    this->samples_bot.width(sdsl::bits::length(this->rindex->pack(this->rindex->index->sequences() - 1, this->rindex->header.max_length - 1)));
    this->samples_bot.resize(total_runs);
    for(size_type i = 0; i < total_runs; i++)
    {
        this->samples_bot[i] = this->rindex->pack(tail_samples[i].seq_id, tail_samples[i].seq_offset);
    }
    sdsl::util::clear(tail_samples);

    // Store the head samples.
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the head samples" << std::endl;
    }
    parallelQuickSort(head_samples.begin(), head_samples.end());
    sdsl::sd_vector_builder builder(this->rindex->index->sequences() * this->rindex->header.max_length, total_runs);
    this->first_to_run.width(sdsl::bits::length(total_runs - 1));
    this->first_to_run.resize(total_runs);
    for(size_type i = 0; i < total_runs; i++)
    {
        builder.set_unsafe(this->rindex->pack(head_samples[i].seq_id, head_samples[i].seq_offset));
        this->first_to_run[i] = head_samples[i].run_id;
    }
    this->first = sdsl::sd_vector<>(builder);


    //lcp_sample[i] is the sample for head_sample[i] (after head_samples is sorted by the sample_record < operator)
    std::vector<size_type> lcp_samples;
    lcp_samples.resize(total_runs);

    // Compute LCP samples.
    double compute_start = readTimer();
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Computing the LCP samples" << std::endl;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for(size_type i = 0; i < this->rindex->index->sequences(); i++)
    {
        size_type start = gbwt::readTimer();
        sample_record srec = {i+1, 0, 0};
        size_type sampIndex = std::lower_bound(head_samples.begin(), head_samples.end(), srec) - head_samples.begin();
        size_type run_id;
        edge_type curr = {ENDMARKER, i};
        range_type run(0,0);
        size_type pathLen = this->pathLength(i)+1; //+1 to count ENDMARKER
        gbwt::vector_type path(pathLen);
        do {
            path[--pathLen] = curr.first;
            edge_type next = lfg.LF(curr, run, run_id);
            if (curr.second == run.first){
                //compute lcp
                size_type lcp = 0;
                edge_type above = {curr.first, curr.second-1};
                while (curr.second != 0 && path[pathLen+lcp] != ENDMARKER && above.first == path[pathLen+lcp]){
                    ++lcp;
                    above = lfg.inverseLF(above);
                }
                lcp_samples[--sampIndex] = lcp;
            }
            curr = next;
        } while (curr.first != ENDMARKER);
        #pragma omp critical
        {
        std::cout << "Extracted seq " << i << " in " << gbwt::readTimer() - start << " seconds" << std::endl;
        }
    }
    
    // Store the LCP samples
    if(Verbosity::level >= Verbosity::FULL)
    {
        std::cerr << msgPrefix << "Storing the lcp samples" << std::endl;
    }

    this->samples_lcp.resize(total_runs);
    for (size_type i = 0; i<total_runs; ++i)
        this->samples_lcp[i] = lcp_samples[i];
    sdsl::util::bit_compress(this->samples_lcp);

    sdsl::util::clear(head_samples);
    sdsl::util::clear(lcp_samples);
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - compute_start;
        std::cerr << msgPrefix << "Extracted " << this->samples_lcp.size() << " lcp samples in " << seconds << " seconds" << std::endl;
    }

    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        std::cerr << msgPrefix << "Processed " << lfg.sequences() << " sequences of total length " << lfg.size() << " in " << seconds << " seconds" << std::endl;
    }
}

FastLCP::size_type FastLCP::serialize(std::ostream &out, sdsl::structure_tree_node* v, std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;
    
    written_bytes += sdsl::serialize(this->samples_bot, out, child, "samples_bot");
    written_bytes += sdsl::serialize(this->first, out, child, "first");
    written_bytes += sdsl::serialize(this->first_to_run, out, child, "first_to_run");
    written_bytes += sdsl::serialize(this->samples_lcp, out, child, "samples_lcp");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void FastLCP::load(std::istream& in) {
    sdsl::load(this->samples_bot, in);
    sdsl::load(this->first, in);
    sdsl::load(this->first_to_run, in);
    sdsl::load(this->samples_lcp, in);
}

const std::string FastLCP::EXTENSION = ".flcp";
//------------------------------------------------------------------------------

#endif //GBWT_QUERY_FAST_LCP_H
