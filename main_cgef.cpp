//
// Created by huangzhibo on 2021/12/14.
//

#include <unistd.h>
#include "main_cgef.h"

int cgef(int argc, char *argv[]) {

    cxxopts::Options options("geftools cgef",
                       "About:  Generate cell bin GEF (.cgef) according to"
                       " common bin GEF (.bgef) file and mask file\n");
    options
    .set_width(120)
    .add_options()
    ("i,input-file", "input bin GEF file [request]", cxxopts::value<std::string>(), "FILE")
    ("m,mask-file", "input mask file [request]", cxxopts::value<std::string>(), "FILE")
    ("o,output-file", "output cell bin GEF file (.cgef) [request]", cxxopts::value<std::string>(), "FILE")
//    ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (argc <= 1 || result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    if (result.count("mask-file") != 1){
        std::cout << "[ERROR] The -m,--mask-file parameter must be given correctly.\n" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }

    if (result.count("input-file") != 1){
        std::cout << "[ERROR] The -i,--input-file parameter must be given correctly.\n" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }

    if (result.count("output-file") != 1){
        std::cout << "[ERROR] The -o,--output-file parameter must be given correctly.\n" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }

    cgefOptions opts = {
        result["input-file"].as<string>(),
        result["mask-file"].as<string>(),
        result["output-file"].as<string>(),
//        result["threads"].as<int>(),
    };
    opts.verbose = result["verbose"].as<bool>();

    time_t prev;
    time(&prev);
    unsigned long cprev=clock();
    cellBinWriter(opts.input_file, opts.mask_file, opts.output_file);
    if(opts.verbose){
        prev = printTime(prev, "cellBinWriter");
        cprev = printCpuTime(cprev, "cellBinWriter");
    }

    return 0;
}
