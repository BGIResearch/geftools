//
// Created by huangzhibo on 2021/12/16.
//

#include "main_view.h"
#include "cxxopts.h"
int test1()
{
    CgefReader cgef_reader = CgefReader("/jdfssz2/ST_BIGDATA/Stomics/auto_analysis/tmppath/users/st_stomics_uat/P20Z10200N0039/S2022042010007/FP200000364TL_D1/FP200000364TL_D1_result/FP200000364TL_D1.cellbin.gef", true);
    // cgef_reader.restrictRegion(7680, 9728, 8704, 10752);
    // cgef_reader.getCellBorders(true, 0);
    
    CellData *arryptr = cgef_reader.getCell();
    int cellnum = cgef_reader.getCellNum();
    for(int i=0;i<cellnum;i++)
    {
        if(arryptr[i].gene_count == 0)
        {
            printf("%d\n", i);
        }
    }

    printf("end\n");
    return 0;
}
int test2()
{
    BgefReader bgef_reader("/ldfssz1/ST_BI/USER/zhaozijian/geftool/build/FP200000364TL_D1.raw.gef",1);
    Expression *arrptr = bgef_reader.getExpression();
    int num = bgef_reader.getExpressionNum();
    int kk = 0;
    for(int i=0;i<num;i++)
    {
        if(arrptr[i].x >= 8857 && arrptr[i].x < 8866  &&
        arrptr[i].y >= 9071 && arrptr[i].y < 9080)
        {
            kk++;
        }
    }
    printf("%d\n",kk);
    return 0;
}

int view(int argc, char *argv[]) {
    //return test2();
    cxxopts::Options options("geftools view",
                             "About:  Show the contents of cell bin GEF\n");
    //TODO support restrict gene_list and region for bGEF 
    options
        .set_width(120)
        .add_options()
            ("i,input-file", "Input bGEF/cGEF file [request]", cxxopts::value<std::string>(), "FILE")
            ("o,output-gem", "Output gem file",
                    cxxopts::value<std::string>()->default_value("stdout"), "FILE")
//            ("m,output-mask", "Output border of polygons to mask format file",
//                    cxxopts::value<std::string>()->default_value(""), "FILE")
            ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list "
                         "of two vertex coordinates (minX,maxX,minY,maxY). just support cGEF.",
                         cxxopts::value<std::string>()->default_value(""), "STR")
            ("g,genes", "Comma separated list of genes to include (or exclude with \"^\" prefix). just support cGEF.",
                    cxxopts::value<std::string>(), "[^]STR")
            ("G,genes-file", "File of genes to include (or exclude with \"^\" prefix)). just support cGEF.",
                    cxxopts::value<std::string>(), "[^]FILE")
            ("force-genes", "Only warn about unknown subset genes, just support cGEF.",
                    cxxopts::value<bool>()->default_value("false"))
            ("b,bin-size", "Set bin size for bgef file, just support bGEF.", cxxopts::value<int>()->default_value("1"), "INT")
            ("s,serial-number", "Serial number", cxxopts::value<std::string>(), "STR")
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
    if (result.count("serial-number") != 1){
        std::cerr << "[ERROR] The -s,--serial-number parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    string snstr = result["serial-number"].as<string>();
    ViewOptions viewopts;

    viewopts.input_file = result["input-file"].as<string>();
    viewopts.output_gem = result["output-gem"].as<string>();
//    viewopts.output_mask = result["output-mask"].as<string>();
    viewopts.verbose = result["verbose"].as<bool>();
    viewopts.force_genes = result["force-genes"].as<bool>();
    viewopts.bin_size = result["bin-size"].as<int>();

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

void toGem(){

}
