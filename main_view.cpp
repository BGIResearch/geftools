//
// Created by huangzhibo on 2021/12/16.
//

#include "main_view.h"
#include "cxxopts.h"

int view(int argc, char *argv[]) {

    cxxopts::Options options("geftools view",
                             "About:  Show the contents of cell bin GEF\n");
    options
        .set_width(70)
        .add_options()
            ("i,input-file", "input cell bin GEF file [request]", cxxopts::value<std::string>(), "FILE")
            ("g,output-gem", "output cell bin gem file", cxxopts::value<std::string>(), "FILE")
            ("p,output-polygon-mask", "output border of polygons to mask format file", cxxopts::value<std::string>(), "FILE")
            ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
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

    ViewOptions view_options = {
        result["input-file"].as<string>(),
        result["mask-file"].as<string>(),
        result["output-file"].as<string>(),
        result["threads"].as<int>(),
        result["verbose"].as<bool>()
    };

    return 0;
}
