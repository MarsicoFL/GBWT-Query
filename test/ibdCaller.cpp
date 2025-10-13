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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <tuple>
#include <utility>
#include <vector>

#include <gbwt/fast_locate.h>
#include <gbwt/gbwt.h>
#include <gbwt/metadata.h>

#include "compText.h"
#include "fast_lcp.h"
#include "lf_gbwt.h"
#include "longMatchQuery.h"

namespace
{

struct Options
{
    std::string base_name;
    gbwt::size_type min_length = 0;
    gbwt::size_type max_haplotypes = std::numeric_limits<gbwt::size_type>::max();
    bool verbose = false;
    int version = 4; // Supported: 2, 3, 4, 24
    std::string output_path;
};

void usage(const char* program)
{
    std::cerr << "Usage: " << program << " <index_base> <min_length> [options]\n\n"
              << "Options:\n"
              << "  --max-haps <N>   Limit the number of forward haplotypes processed (default: all)\n"
              << "  --version <V>    Select long match implementation (2, 3, 4, or 24; default: 4)\n"
              << "  --output <path>  Write results to file instead of stdout\n"
              << "  --verbose        Report progress to stderr\n";
}

bool parse_arguments(int argc, char** argv, Options& options)
{
    if (argc < 3)
    {
        usage(argv[0]);
        return false;
    }

    options.base_name = argv[1];
    try
    {
        options.min_length = std::stoull(argv[2]);
    }
    catch (const std::exception&)
    {
        std::cerr << "Invalid minimum length: '" << argv[2] << "'\n";
        return false;
    }

    for (int i = 3; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--max-haps")
        {
            if (++i == argc)
            {
                std::cerr << "--max-haps requires an argument\n";
                return false;
            }
            try
            {
                options.max_haplotypes = std::stoull(argv[i]);
            }
            catch (const std::exception&)
            {
                std::cerr << "Invalid --max-haps argument: '" << argv[i] << "'\n";
                return false;
            }
        }
        else if (arg == "--version")
        {
            if (++i == argc)
            {
                std::cerr << "--version requires an argument\n";
                return false;
            }
            try
            {
                options.version = std::stoi(argv[i]);
            }
            catch (const std::exception&)
            {
                std::cerr << "Invalid --version argument: '" << argv[i] << "'\n";
                return false;
            }
            if (options.version != 2 && options.version != 3 && options.version != 4 && options.version != 24)
            {
                std::cerr << "Unsupported version " << options.version << ". Expected 2, 3, 4, or 24.\n";
                return false;
            }
        }
        else if (arg == "--output")
        {
            if (++i == argc)
            {
                std::cerr << "--output requires a path argument\n";
                return false;
            }
            options.output_path = argv[i];
        }
        else if (arg == "--verbose")
        {
            options.verbose = true;
        }
        else
        {
            std::cerr << "Unrecognized argument: '" << arg << "'\n";
            usage(argv[0]);
            return false;
        }
    }

    if (options.min_length == 0)
    {
        std::cerr << "Minimum length must be > 0\n";
        return false;
    }

    return true;
}

struct Segment
{
    gbwt::size_type hap1_start = 0;
    gbwt::size_type hap2_start = 0;
    gbwt::size_type length = 0;
};

bool operator<(const Segment& a, const Segment& b)
{
    if (a.hap1_start != b.hap1_start) return a.hap1_start < b.hap1_start;
    if (a.hap2_start != b.hap2_start) return a.hap2_start < b.hap2_start;
    return a.length < b.length;
}

bool operator==(const Segment& a, const Segment& b)
{
    return a.hap1_start == b.hap1_start && a.hap2_start == b.hap2_start && a.length == b.length;
}

template<class MatchVector>
void record_matches(const MatchVector& matches,
                    gbwt::size_type query_hap,
                    gbwt::size_type orientation_factor,
                    gbwt::size_type hap_count,
                    std::map<std::pair<gbwt::size_type, gbwt::size_type>, std::vector<Segment>>& output)
{
    for (const auto& match : matches)
    {
        const gbwt::size_type query_start = std::get<0>(match);
        const gbwt::size_type length = std::get<1>(match);
        const gbwt::size_type target_sequence = std::get<2>(match);
        const gbwt::size_type target_start = std::get<3>(match);

        if (length == 0) { continue; }

        gbwt::size_type target_orientation = (orientation_factor == 1 ? 0 : (target_sequence % orientation_factor));
        gbwt::size_type target_hap = (orientation_factor == 1 ? target_sequence : target_sequence / orientation_factor);

        if (target_hap >= hap_count) { continue; }
        if (target_hap == query_hap) { continue; }
        if (target_hap < query_hap) { continue; }
        if (orientation_factor > 1 && target_orientation != 0) { continue; }

        Segment segment{query_start, target_start, length};
        output[{query_hap, target_hap}].push_back(segment);
    }
}

void deduplicate(std::map<std::pair<gbwt::size_type, gbwt::size_type>, std::vector<Segment>>& output)
{
    for (auto& entry : output)
    {
        auto& segments = entry.second;
        std::sort(segments.begin(), segments.end());
        segments.erase(std::unique(segments.begin(), segments.end()), segments.end());
    }
}

bool file_exists(const std::string& path)
{
    struct stat buffer{};
    return (::stat(path.c_str(), &buffer) == 0);
}

} // namespace

int main(int argc, char** argv)
{
    Options options;
    if (!parse_arguments(argc, argv, options))
    {
        return 1;
    }

    if (options.verbose)
    {
        std::cerr << "Loading indexes with base '" << options.base_name << "'\n";
    }

    gbwt::GBWT gbwt_index;
    {
        const std::string gbwt_path = options.base_name + gbwt::GBWT::EXTENSION;
        std::ifstream input(gbwt_path, std::ios::binary);
        if (!input)
        {
            std::cerr << "Failed to open required GBWT file '" << gbwt_path << "'\n";
            return 1;
        }
        if (options.verbose)
        {
            std::cerr << "Loading GBWT index from '" << gbwt_path << "'\n";
        }
        gbwt_index.load(input);
    }

    gbwt::FastLocate locate_index;
    {
        const std::string locate_path = options.base_name + gbwt::FastLocate::EXTENSION;
        if (file_exists(locate_path))
        {
            std::ifstream input(locate_path, std::ios::binary);
            locate_index.load(input);
            if (options.verbose)
            {
                std::cerr << "Loaded FastLocate from '" << locate_path << "'\n";
            }
        }
        else
        {
            if (options.verbose)
            {
                std::cerr << "Building FastLocate in memory" << std::endl;
            }
            locate_index = gbwt::FastLocate(gbwt_index);
        }
        locate_index.index = &gbwt_index;
    }

    FastLCP lcp_index;
    {
        const std::string lcp_path = options.base_name + FastLCP::EXTENSION;
        if (file_exists(lcp_path))
        {
            std::ifstream input(lcp_path, std::ios::binary);
            lcp_index.load(input);
            if (options.verbose)
            {
                std::cerr << "Loaded FastLCP from '" << lcp_path << "'\n";
            }
        }
        else
        {
            if (options.verbose)
            {
                std::cerr << "Building FastLCP in memory" << std::endl;
            }
            lcp_index = FastLCP(locate_index);
        }
        lcp_index.rindex = &locate_index;
    }

    lf_gbwt::GBWT lf_index;
    const bool needs_lf = (options.version == 3 || options.version == 4);
    if (needs_lf)
    {
        const std::string lf_path = options.base_name + lf_gbwt::GBWT::EXTENSION;
        if (file_exists(lf_path))
        {
            std::ifstream input(lf_path, std::ios::binary);
            lf_index.load(input);
            if (options.verbose)
            {
                std::cerr << "Loaded lf_gbwt::GBWT from '" << lf_path << "'\n";
            }
        }
        else
        {
            if (options.verbose)
            {
                std::cerr << "Building lf_gbwt::GBWT in memory" << std::endl;
            }
            lf_index = lf_gbwt::GBWT(gbwt_index);
        }
    }

    CompText comp_text;
    const bool needs_comp_text = (options.version == 4 || options.version == 24);
    if (needs_comp_text)
    {
        const std::string comp_path = options.base_name + CompText::EXTENSION;
        if (file_exists(comp_path))
        {
            std::ifstream input(comp_path, std::ios::binary);
            comp_text.load(input);
            if (options.verbose)
            {
                std::cerr << "Loaded CompText from '" << comp_path << "'\n";
            }
        }
        else
        {
            if (options.verbose)
            {
                std::cerr << "Building CompText in memory" << std::endl;
            }
            comp_text.buildFullMemPruned(lcp_index);
        }
        comp_text.source = &lcp_index;
    }

    if (!gbwt_index.bidirectional())
    {
        if (options.version == 3 || options.version == 4)
        {
            std::cerr << "Selected version requires a bidirectional GBWT\n";
            return 1;
        }
    }

    const gbwt::size_type orientation_factor = (gbwt_index.bidirectional() ? 2 : 1);
    const gbwt::size_type hap_count = gbwt_index.sequences() / orientation_factor;

    if (hap_count == 0)
    {
        std::cerr << "GBWT contains no haplotypes\n";
        return 1;
    }

    const gbwt::size_type hap_limit = std::min(options.max_haplotypes, hap_count);

    std::ofstream output_file;
    std::ostream* out = &std::cout;
    if (!options.output_path.empty())
    {
        output_file.open(options.output_path);
        if (!output_file)
        {
            std::cerr << "Failed to open output path '" << options.output_path << "'\n";
            return 1;
        }
        out = &output_file;
    }

    std::map<std::pair<gbwt::size_type, gbwt::size_type>, std::vector<Segment>> pair_segments;

    for (gbwt::size_type hap = 0; hap < hap_limit; ++hap)
    {
        if (options.verbose)
        {
            std::cerr << "Processing haplotype " << hap + 1 << " / " << hap_limit << "\n";
        }

        gbwt::size_type sequence = (orientation_factor == 2 ? gbwt::Path::encode(hap, false) : hap);
        gbwt::vector_type query = gbwt_index.extract(sequence);

        if (query.empty())
        {
            if (options.verbose)
            {
                std::cerr << "Skipping haplotype " << hap << " (empty path)\n";
            }
            continue;
        }

        std::vector<std::tuple<gbwt::size_type, gbwt::size_type, gbwt::size_type, gbwt::size_type>> matches;
        switch (options.version)
        {
            case 2:
                matches = longMatchQuery2(gbwt_index, locate_index, lcp_index, options.min_length, query);
                break;
            case 3:
                matches = longMatchQuery3(lf_index, locate_index, lcp_index, options.min_length, query);
                break;
            case 4:
                matches = longMatchQuery4(lf_index, locate_index, lcp_index, comp_text, options.min_length, query);
                break;
            case 24:
                matches = longMatchQuery2_4(gbwt_index, locate_index, lcp_index, comp_text, options.min_length, query);
                break;
            default:
                throw std::logic_error("Unsupported version encountered at runtime");
        }

        record_matches(matches, hap, orientation_factor, hap_count, pair_segments);
    }

    deduplicate(pair_segments);

    *out << "#hap1_id\thap2_id\thap1_start\thap2_start\tlength_nodes\n";

    gbwt::size_type reported_pairs = 0;
    gbwt::size_type reported_segments = 0;
    for (const auto& entry : pair_segments)
    {
        const auto& hap_pair = entry.first;
        const auto& segments = entry.second;
        if (segments.empty()) { continue; }

        ++reported_pairs;
        for (const Segment& segment : segments)
        {
            ++reported_segments;
            *out << hap_pair.first << '\t'
                 << hap_pair.second << '\t'
                 << segment.hap1_start << '\t'
                 << segment.hap2_start << '\t'
                 << segment.length << '\n';
        }
    }

    if (options.verbose)
    {
        std::cerr << "Reported " << reported_segments << " IBD candidate segments across "
                  << reported_pairs << " haplotype pairs\n";
    }

    return 0;
}
