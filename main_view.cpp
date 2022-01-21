//
// Created by huangzhibo on 2021/12/16.
//

#include "main_view.h"
#include "cxxopts.h"

int view(int argc, char *argv[]) {

    cxxopts::Options options("geftools view",
                             "About:  Show the contents of cell bin GEF\n");
    options
        .set_width(120)
        .add_options()
            ("i,input-file", "Input cell bin GEF file [request]", cxxopts::value<std::string>(), "FILE")
            ("o,output-gem", "Output cell bin gem file",
                    cxxopts::value<std::string>()->default_value("stdout"), "FILE")
            ("m,output-mask", "Output border of polygons to mask format file",
                    cxxopts::value<std::string>()->default_value(""), "FILE")
            ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list "
                         "of two vertex coordinates (minX,maxX,minY,maxY)",
                         cxxopts::value<std::string>()->default_value(""), "STR")
            ("g,genes", "Comma separated list of genes to include (or exclude with \"^\" prefix)",
                    cxxopts::value<std::string>(), "[^]STR")
            ("G,genes-file", "File of genes to include (or exclude with \"^\" prefix))",
                    cxxopts::value<std::string>(), "[^]FILE")
            ("force-genes", "Only warn about unknown subset genes",
                    cxxopts::value<bool>()->default_value("false"))
//            ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
            ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (argc <= 1 || result.count("help"))
    {
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    if (result.count("input-file") != 1){
        std::cerr << "[ERROR] The -i,--input-file parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    ViewOptions viewopts;

    viewopts.input_file = result["input-file"].as<string>();
    viewopts.output_gem = result["output-gem"].as<string>();
    viewopts.output_mask = result["output-mask"].as<string>();
    viewopts.verbose = result["verbose"].as<bool>();
    viewopts.force_genes = result["force-genes"].as<bool>();

    if (result.count("region") == 1){
        string region_tmp = result["region"].as<string>();
        vector<string> regions = split(region_tmp, ',');
        viewopts.restrict_region = true;
        viewopts.region[0] = static_cast<unsigned int>(strtol(regions[0].c_str(), nullptr, 10));
        viewopts.region[1] = static_cast<unsigned int>(strtol(regions[1].c_str(), nullptr, 10));
        viewopts.region[2] = static_cast<unsigned int>(strtol(regions[2].c_str(), nullptr, 10));
        viewopts.region[3] = static_cast<unsigned int>(strtol(regions[3].c_str(), nullptr, 10));
    }

    unsigned long n_genes = result.count("genes");
    if (n_genes == 1){
        string genes_tmp = result["genes"].as<string>();
        string::iterator it = genes_tmp.begin();
        if(*it == '^') {
            genes_tmp.erase(it);
            viewopts.exclude = true;
        }else {
            viewopts.exclude = false;
        }
        viewopts.genes = split(genes_tmp, ',');
    }else if(n_genes > 1){
        std::cerr << "[ERROR] The -g,--genes parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    unsigned long n_genes_file = result.count("genes-file");
    if(n_genes_file == 1){
        if(n_genes == 1){
            std::cerr << "[ERROR] You cannot specify both -g,--genes "
                         "and -G,--genes-file parameters at the same time.\n" << std::endl;
            exit(1);
        }
        string genes_file = result["genes-file"].as<string>();
        viewopts.genes = readLines(genes_file);
    }else if(n_genes_file > 1){
        std::cerr << "[ERROR] The -G,--genes-file parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    time_t prev;
    time(&prev);
    unsigned long cprev=clock();

    CgefReader cgef_reader = CgefReader(viewopts.input_file, viewopts.verbose);

    if(viewopts.restrict_region){
        cgef_reader.restrictRegion(viewopts.region[0],
                                   viewopts.region[1],
                                   viewopts.region[2],
                                   viewopts.region[3]);
    }

    cgef_reader.toGem(viewopts.output_gem, viewopts.genes, viewopts.force_genes, viewopts.exclude);

    if(viewopts.verbose){
        prev = printTime(prev, "CgefReader init");
        cprev = printCpuTime(cprev, "CgefReader init");
    }

    return 0;
}

void toGem(){

}
