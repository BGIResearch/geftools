//
// Created by huangzhibo on 2021/12/14.
//

#include "main_gen.h"

int gen(int argc, char *argv[]) {

    cxxopts::Options options("geftools gen",
                       "About:  Generate cell bin GEF according to common bin GEF file and mask file\n");
    options
    .set_width(70)
    .add_options()
    ("i,input-file", "input bin GEF file [request]", cxxopts::value<std::string>(), "FILE")
    ("m,mask-file", "input mask file [request]", cxxopts::value<std::string>(), "FILE")
    ("o,output-file", "output cell bin GEF file [request]", cxxopts::value<std::string>(), "FILE")
//    ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
//    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
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

    genOptions gen_options = {
        result["input-file"].as<string>(),
        result["mask-file"].as<string>(),
        result["output-file"].as<string>(),
//        result["threads"].as<int>(),
//        result["verbose"].as<bool>()
    };

    cellBinWriter(gen_options.input_file, gen_options.mask_file, gen_options.output_file);

    return 0;
}
