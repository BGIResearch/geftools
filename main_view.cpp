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
            ("o,output-gem", "Output cell bin gem file", cxxopts::value<std::string>(), "FILE")
            ("m,output-mask", "Output border of polygons to mask format file", cxxopts::value<std::string>(), "FILE")
            ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list "
                         "of two vertex coordinates (minX,minY,maxX,maxY)", cxxopts::value<std::string>(), "STR")
            ("g,genes", "Comma separated list of genes to include (or exclude with \"^\" prefix)",
                    cxxopts::value<std::string>(), "[^]STR")
            ("G,genes-file", "File of genes to include (or exclude with \"^\" prefix))",
                    cxxopts::value<std::string>(), "[^]FILE")
//            ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
//            ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
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

    if (result.count("output-gem") != 1){
        std::cerr << "[ERROR] The -g,--output-gem parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    if (result.count("output-mask") != 1){
        std::cerr << "[ERROR] The -m,--output-mask parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }


    ViewOptions view_options = {
        result["input-file"].as<string>(),
        result["output-gem"].as<string>(),
        result["output-mask"].as<string>(),
//        result["threads"].as<int>(),
//        result["verbose"].as<bool>()
    };

    if (result.count("region") != 1){
        string region_tmp = result["region"].as<string>();
        view_options.region = split(region_tmp, ',');
    }

    unsigned long n_genes = result.count("genes");
    if (n_genes == 1){
        string genes_tmp = result["genes"].as<string>();
        view_options.genes = split(genes_tmp, ',');
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
        view_options.genes = readLines(genes_file);
    }else if(n_genes_file > 1){
        std::cerr << "[ERROR] The -G,--genes-file parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    return 0;
}

void toGem(){

}
