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

#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<string>
#include<exception>
#include<random>
#include<algorithm>


#include<gbwt/gbwt.h>
#include<gbwt/dynamic_gbwt.h>
#include<gbwt/fast_locate.h>
#include<fast_lcp.h>
#include<lf_gbwt.h>
#include<compText.h>

template <class T>
T safeOpenFile(const std::string & name) {
    T file(name);
    if (!file.is_open()) {
        std::cout << "File '" << name << "' failed to open!\n";
        exit(1);
    }
    return file;
}

std::vector<gbwt::size_type> generateNUniqueRandomLessThanM(gbwt::size_type N, gbwt::size_type M){
    if (N > M)
        throw std::invalid_argument("N must be less than or equal to M!");
    std::vector<gbwt::size_type> ans(M);
    for (gbwt::size_type i = 0; i < M; ++i)
        ans[i] = i;

    std::random_device rd;
    std::mt19937_64 g(rd());
    std::shuffle(ans.begin(), ans.end(), g);
    ans.resize(N);
    return ans;
}

int main(int argc, char* argv[]) {
    gbwt::size_type haplotypes;
    std::string suffix = "Trimmed";
    if (argc != 3) {
        std::cout << "There must be exactly two arguments passed, first the base name for all the data structures, then the number of random haplotypes to remove from the GBWT and write new structures. " << argc - 1 << " arguments passed." << std::endl;
        for (int i = 0; i < argc; ++i)
            std::cout << "Argument " << i << ": '" << argv[i] << "'" << std::endl;
        exit(1);
    }
    haplotypes = std::stoull(argv[2]);
    if (haplotypes == 0){
        std::cout << "Number of haplotypes must be > 0, haplotypes: " << haplotypes << ". Arg passed: " << argv[2] << std::endl;
        exit(1);
    }

    double start = gbwt::readTimer();
    gbwt::size_type bytes;

    //open and load gbwt
    gbwt::GBWT g;
    std::ifstream gin = safeOpenFile<std::ifstream>(argv[1] + gbwt::GBWT::EXTENSION); 
    g.load(gin);
    gin.close(); 
    std::cout << "Read GBWT in " << gbwt::readTimer() - start << " seconds" << std::endl;
    gbwt::printStatistics(g, argv[1] + gbwt::GBWT::EXTENSION, std::cout);
    
    gbwt::size_type hapsInGBWT = g.sequences() / (1+g.bidirectional());
    if (haplotypes >= hapsInGBWT){
        std::cout << "Haplotypes to remove must be less than the number of haplotypes in the GBWT. Haplotypes to remove: " << haplotypes << ". Haplotypes in the GBWT: " << hapsInGBWT << std::endl;
        exit(1);
    }

    //building dynamic gbwwt
    start = gbwt::readTimer();
    gbwt::DynamicGBWT dg(g);
    std::cout << "Finished constructing dynamic GBWT in " << gbwt::readTimer() - start << " seconds." << std::endl;

    //getting haplotypes to remove
    start = gbwt::readTimer();
    std::vector<gbwt::size_type> hapsToRemove = generateNUniqueRandomLessThanM(haplotypes, hapsInGBWT);
    std::cout << "Finished selecting haplotypes to remove in " << gbwt::readTimer() - start << " seconds." << std::endl;

    //store raw haplotype paths to query
    start = gbwt::readTimer();
    sdsl::int_vector<> queryPaths;
    gbwt::size_type numPaths = 0;
    for (auto seqId : hapsToRemove){
        gbwt::vector_type path;
        gbwt::size_type size;
        if (g.bidirectional()){
            //forward orientation
            size = queryPaths.size();
            path = g.extract(gbwt::Path::encode(seqId, false));
            queryPaths.resize(size+path.size()+1);
            for (gbwt::size_type i = 0; i != path.size(); ++i)
                queryPaths[i+size] = path[i];
            queryPaths[path.size()+size] = gbwt::ENDMARKER;
            //reverse orientation
            size = queryPaths.size();
            path = g.extract(gbwt::Path::encode(seqId, true));
            queryPaths.resize(size+path.size()+1);
            for (gbwt::size_type i = 0; i != path.size(); ++i)
                queryPaths[i+size] = path[i];
            queryPaths[path.size()+size] = gbwt::ENDMARKER;

            numPaths += 2;
        }
        else{
            size = queryPaths.size();
            path = g.extract(seqId);
            queryPaths.resize(size+path.size()+1);
            for (gbwt::size_type i = 0; i != path.size(); ++i)
                queryPaths[i+size] = path[i];
            queryPaths[path.size()+size] = gbwt::ENDMARKER;

            numPaths += 1;
        }
    }
    std::ofstream queryOut = safeOpenFile<std::ofstream>(argv[1] + suffix + "QueryPaths");
    bytes = sdsl::serialize(hapsToRemove, queryOut);
    bytes += queryPaths.serialize(queryOut);
    std::cout << "Finished storing " << numPaths << " sequences of total length " << queryPaths.size() << " to remove in " << gbwt::readTimer() - start << " seconds and " << bytes << " space to " << argv[1] + suffix + "QueryPaths"  << std::endl;

    //remove haplotypes from dynamic gbwt
    start = gbwt::readTimer();
    dg.metadata.clear();
    dg.remove(hapsToRemove);
    std::cout << "Finished removing haplotypes to in " << gbwt::readTimer() - start << " seconds." << std::endl;
    
    //build compressed gbwt
    start = gbwt::readTimer();
    gbwt::GBWT trimG(dg);
    std::cout << "Finished building new GBWT in " << gbwt::readTimer() - start << " seconds." << std::endl;
    //writing compressed gbwt
    start = gbwt::readTimer();
    std::ofstream trimGout = safeOpenFile<std::ofstream>(argv[1] + suffix + gbwt::GBWT::EXTENSION);
    bytes = trimG.serialize(trimGout);
    trimGout.close();
    std::cout << "Finished writing new GBWT of size " << bytes << " bytes in " << gbwt::readTimer() - start << " seconds to " << argv[1] + suffix + gbwt::GBWT::EXTENSION << std::endl;
    
    //build FastLocate
    start = gbwt::readTimer();
    gbwt::FastLocate trimfl(trimG);
    std::cout << "Finished building new FastLocate in " << gbwt::readTimer() - start << " seconds." << std::endl;
    //writing FastLocate
    start = gbwt::readTimer();
    std::ofstream trimflout = safeOpenFile<std::ofstream>(argv[1] + suffix + gbwt::FastLocate::EXTENSION);
    bytes = trimfl.serialize(trimflout);
    trimflout.close();
    std::cout << "Finished writing new FastLocate of size " << bytes << " bytes in " << gbwt::readTimer() - start << " seconds to " << argv[1] + suffix + gbwt::FastLocate::EXTENSION << std::endl;

    //build FastLCP
    start = gbwt::readTimer();
    FastLCP trimflcp(trimfl);
    std::cout << "Finished building new FastLCP in " << gbwt::readTimer() - start << " seconds." << std::endl;
    //writing FastLCP
    start = gbwt::readTimer();
    std::ofstream trimflcpout = safeOpenFile<std::ofstream>(argv[1] + suffix + FastLCP::EXTENSION);
    bytes = trimflcp.serialize(trimflcpout);
    trimflcpout.close();
    std::cout << "Finished writing new FastLCP of size " << bytes << " bytes in " << gbwt::readTimer() - start << " seconds to " << argv[1] + suffix + FastLCP::EXTENSION << std::endl;

    //build lf gbwt
    start = gbwt::readTimer();
    lf_gbwt::GBWT trimlfg(trimG);
    std::cout << "Finished building new lf gbwt in " << gbwt::readTimer() - start << " seconds." << std::endl;
    //writing lf gbwt
    start = gbwt::readTimer();
    std::ofstream trimlfgout = safeOpenFile<std::ofstream>(argv[1] + suffix + lf_gbwt::GBWT::EXTENSION);
    bytes = trimlfg.serialize(trimlfgout);
    trimlfgout.close();
    std::cout << "Finished writing new lf gbwt of size " << bytes << " bytes in " << gbwt::readTimer() - start << " seconds to " << argv[1] + suffix + lf_gbwt::GBWT::EXTENSION << std::endl;

    //build comp text
    start = gbwt::readTimer();
    CompText trimct;
    trimct.buildFullMem(trimflcp);
    std::cout << "Finished building new comp text in " << gbwt::readTimer() - start << " seconds." << std::endl;
    //writing comp text
    start = gbwt::readTimer();
    std::ofstream trimctout = safeOpenFile<std::ofstream>(argv[1] + suffix + CompText::EXTENSION);
    bytes = trimct.serialize(trimctout);
    std::cout << "Finished writing new comp text of size " << bytes <<  " bytes in " << gbwt::readTimer() - start << " seconds to " << argv[1] + suffix + CompText::EXTENSION << std::endl;
    

    return 0;
}
