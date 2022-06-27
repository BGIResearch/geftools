//
// Created by huangzhibo on 2021/12/16.
//

#include "main_view.h"
#include "cxxopts.h"
#include "geftogem.h"

int view(int argc, char *argv[]) {
    cxxopts::Options options("geftools view",
                             "About:  Show the contents of cell bin GEF\n");
    //TODO support restrict gene_list and region for bGEF 
    options
        .set_width(120)
        .add_options()
            ("i,input-file", "Input bGEF/cGEF file [request]", cxxopts::value<std::string>(), "FILE")
            ("o,output-gem", "Output gem file ",
                    cxxopts::value<std::string>()->default_value("stdout"), "FILE")
           ("d,exp_data", "Input bgef for cgem",
                   cxxopts::value<std::string>()->default_value(""), "FILE")
            // ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list "
            //              "of two vertex coordinates (minX,maxX,minY,maxY). just support cGEF.",
            //              cxxopts::value<std::string>()->default_value(""), "STR")
            // ("g,genes", "Comma separated list of genes to include (or exclude with \"^\" prefix). just support cGEF.",
            //         cxxopts::value<std::string>(), "[^]STR")
            // ("G,genes-file", "File of genes to include (or exclude with \"^\" prefix)). just support cGEF.",
            //         cxxopts::value<std::string>(), "[^]FILE")
            // ("force-genes", "Only warn about unknown subset genes, just support cGEF.",
            //         cxxopts::value<bool>()->default_value("false"))
            ("b,bin-size", "Set bin size for bgef file, just support bGEF.", cxxopts::value<int>()->default_value("1"), "INT")
            ("s,serial-number", "Serial number [request]", cxxopts::value<std::string>(), "STR")
//            ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
            //("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("e,exon", "whether or not output exon", cxxopts::value<int>()->default_value("1"), "INT")
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
    if (result.count("serial-number") != 1){
        std::cerr << "[ERROR] The -s,--serial-number parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    // if(result.count("output-gem") != 1)
    // {
    //     std::cerr << "[ERROR] The -o,--output-gem parameter must be given correctly.\n" << std::endl;
    //     std::cerr << options.help() << std::endl;
    //     exit(1);
    // }
    bool boutexon = result["exon"].as<int>();
    string strin = result["input-file"].as<string>();
    string snstr = result["serial-number"].as<string>();
    string strout = result["output-gem"].as<string>();
    

    geftogem gem(strout, snstr, boutexon);
    if(is_bgef(strin))
    {
        gem.bgeftogem(strin);
    }
    else
    {
        string strexp = result["exp_data"].as<string>();
        gem.cgeftogem(strin, strexp);
    }
    
    return 0;

    ViewOptions viewopts;

    viewopts.input_file = result["input-file"].as<string>();
    viewopts.output_gem = result["output-gem"].as<string>();
//    viewopts.output_mask = result["output-mask"].as<string>();
    viewopts.verbose = result["verbose"].as<bool>();
    viewopts.force_genes = result["force-genes"].as<bool>();
    viewopts.bin_size = result["bin-size"].as<int>();
    string strdata = result["data"].as<string>();

    if (result.count("region") == 1){
        string region_tmp = result["region"].as<string>();
        vector<string> regions = split(region_tmp, ',');
        viewopts.restrict_region = true;
        for(int i =0; i < 4; i++){
            viewopts.region[i] = static_cast<unsigned int>(strtol(regions[i].c_str(), nullptr, 10));
        }
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

    if(is_bgef(viewopts.input_file)){
        BgefReader bgef_reader = BgefReader(viewopts.input_file, viewopts.bin_size, 1, viewopts.verbose);
        bgef_reader.toGem(viewopts.output_gem, snstr);
    }else{

        CgefReader cgef_reader = CgefReader(viewopts.input_file, viewopts.verbose);

        if(viewopts.restrict_region){
            cgef_reader.restrictRegion(viewopts.region[0],
                                       viewopts.region[1],
                                       viewopts.region[2],
                                       viewopts.region[3]);
        }

        cgef_reader.toGem(viewopts.output_gem, viewopts.genes, viewopts.force_genes, viewopts.exclude);
    }

    if(viewopts.verbose){
        prev = printTime(prev, "CgefReader init");
        cprev = printCpuTime(cprev, "CgefReader init");
    }

    return 0;
}
