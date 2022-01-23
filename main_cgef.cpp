//
// Created by huangzhibo on 2021/12/14.
//

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
    ("b,block", "Pre block size", cxxopts::value<std::string>()->default_value("256,256"), "FILE")
    ("r,rand-celltype", "number of random cell type", cxxopts::value<int>()->default_value("0"), "INT")
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
        result["rand-celltype"].as<int>(),
//        result["threads"].as<int>(),
    };
    opts.verbose = result["verbose"].as<bool>();
    vector<string> block_size_tmp = split(result["block"].as<string>(), ',');

    if(block_size_tmp.size() != 2){
        std::cout << "[ERROR] The -b,--block parameter must be given correctly.\n" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }
    opts.block_size[0] = static_cast<int>(strtol(block_size_tmp[0].c_str(), nullptr, 10));
    opts.block_size[1] = static_cast<int>(strtol(block_size_tmp[1].c_str(), nullptr, 10));
    generateCgef(opts.output_file, opts.input_file, opts.mask_file, opts.block_size, opts.verbose);
    return 0;
}

int generateCgef(const string &cgef_file,
                 const string &bgef_file,
                 const string &mask_file,
                 const int* block_size,
                 bool verbose) {

    unsigned long cprev=clock();
    Mask mask = Mask(mask_file, block_size);
    if(verbose) cprev = printCpuTime(cprev, "Mask init");
    BgefReader common_bin_gef = BgefReader(bgef_file, 1, true);
    CgefWriter cgef_writer = CgefWriter(cgef_file, true);
    cgef_writer.setRandomCellTypeNum(20);
    cgef_writer.write(common_bin_gef, mask);
    if(verbose) printCpuTime(cprev, "generateCgef");
    return 0;
}
