// Created by huangzhibo on 2021/12/14.
//

#include "main_cgef.h"
#include "opencv2/opencv.hpp"
#include "cgefParam.h"
#include "cgefCellgem.h"

int cgef(int argc, char *argv[]) {

    cxxopts::Options options("geftools cgef",
                       "About:  Generate cell bin GEF (.cgef) according to"
                       " common bin GEF (.bgef) file and mask file\n");
    options
    .set_width(120)
    .add_options()
    ("i,input-file", "input GEF file [request]", cxxopts::value<std::string>(), "FILE")
    ("m,mask-file", "input mask file [request]", cxxopts::value<std::string>(), "FILE")
    ("o,output-file", "output cell bin GEF file (.cgef) [request]", cxxopts::value<std::string>(), "FILE")
    ("b,block", "Pre block size", cxxopts::value<std::string>()->default_value("256,256"), "FILE")
    ("r,rand-celltype", "number of random cell type", cxxopts::value<int>()->default_value("0"), "INT")
    ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("p,patch", "add the path to cgef", cxxopts::value<int>()->default_value("0"))
    ("n,cnum", "top level cell num", cxxopts::value<int>()->default_value("5000"), "INT")
    ("R,ratio", "other level cell num ratio", cxxopts::value<int>()->default_value("20"), "FLOAT")
    ("a,allocat", "allocation strategy ", cxxopts::value<int>()->default_value("2"), "INT")
    ("g,gem", "raw gem file", cxxopts::value<std::string>(), "FILE")
    ("c,canvas", "set canvas size", cxxopts::value<std::string>()->default_value("90000,90000"), "FILE")
    ("l,limit", "set blk limit", cxxopts::value<std::string>()->default_value("16,16"), "FILE")
    ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (argc <= 1 || result.count("help"))
    {
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    CgefOptions opts;
    if (result.count("mask-file") != 1){
        // std::cerr << "[ERROR] The -m,--mask-file parameter must be given correctly.\n" << std::endl;
        // std::cerr << options.help() << std::endl;
        // exit(1);
        opts.mask_file = "";
    }
    else
    {
        opts.mask_file = result["mask-file"].as<string>();
    }

    if (result.count("input-file") != 1){
        std::cerr << "[ERROR] The -i,--input-file parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    if (result.count("output-file") != 1){
        // std::cerr << "[ERROR] The -o,--output-file parameter must be given correctly.\n" << std::endl;
        // std::cerr << options.help() << std::endl;
        // exit(1);
        opts.output_file = "";
    }
    else
    {
        opts.output_file = result["output-file"].as<string>();
    }

//     CgefOptions opts = {
//         result["input-file"].as<string>(),
//         result["mask-file"].as<string>(),
//         result["output-file"].as<string>(),
//         result["rand-celltype"].as<int>(),
// //        result["threads"].as<int>(),
//     };
    int patch = result["patch"].as<int>();
    opts.input_file = result["input-file"].as<string>();
    opts.rand_celltype_num = result["rand-celltype"].as<int>();
    opts.verbose = result["verbose"].as<bool>();
    opts.cellnum = result["cnum"].as<int>();
    int tmpr = result["ratio"].as<int>();
    opts.ratio = tmpr*1.0/100;
    int allocat = result["allocat"].as<int>();
    vector<string> block_size_tmp = split(result["block"].as<string>(), ',');
    vector<string> canvas_size_tmp = split(result["canvas"].as<string>(), ',');
    vector<string> limit_tmp = split(result["limit"].as<string>(), ',');

    if(block_size_tmp.size() != 2){
        std::cerr << "[ERROR] The -b,--block parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    opts.block_size[0] = static_cast<int>(strtol(block_size_tmp[0].c_str(), nullptr, 10));
    opts.block_size[1] = static_cast<int>(strtol(block_size_tmp[1].c_str(), nullptr, 10));
    int canvas_size[2]={0,0};
    canvas_size[0] = static_cast<int>(strtol(canvas_size_tmp[0].c_str(), nullptr, 10));
    canvas_size[1] = static_cast<int>(strtol(canvas_size_tmp[1].c_str(), nullptr, 10));
    int limit_blk[2] = {0,0};
    limit_blk[0] = static_cast<int>(strtol(limit_tmp[0].c_str(), nullptr, 10));
    limit_blk[1] = static_cast<int>(strtol(limit_tmp[1].c_str(), nullptr, 10));
    
    cgefParam::GetInstance()->m_cellgemstr = opts.input_file;
    cgefParam::GetInstance()->m_maskstr = opts.mask_file;
    cgefParam::GetInstance()->m_block_size[0] = opts.block_size[0];
    cgefParam::GetInstance()->m_block_size[1] = opts.block_size[1];
    cgefParam::GetInstance()->m_intype = (InputType)patch;
    switch (patch)
    {
    case 0:
        generateCgef(opts.output_file, opts.input_file, opts.mask_file, opts.block_size,canvas_size,limit_blk,
                    opts.rand_celltype_num, allocat, opts.cellnum, opts.ratio, opts.verbose);
        break;
    case 1:
        break;
    case 2:
        cgefParam::GetInstance()->m_rawgemstr = result["gem"].as<string>();
    case 3:
    {
        cgefParam::GetInstance()->m_threadcnt = result["threads"].as<int>();
        CgefWriter *pcgef_writer = new CgefWriter(true);
        pcgef_writer->setOutput(opts.output_file);
        pcgef_writer->setRandomCellTypeNum(opts.rand_celltype_num);

        cgefCellgem cgem;
        cgem.writeFile(pcgef_writer);

        //pcgef_writer->addLevel(allocat, opts.cellnum, opts.ratio, canvas_size, limit_blk);
        delete pcgef_writer;
    }
        break;
    default:
        break;
    }

    return 0;
}

int generateCgef(const string &cgef_file,
                 const string &bgef_file,
                 const string &mask_file,
                 const int* block_size,
                 int* canvas_size,
                 int* limit_blk,
                 int rand_cell_type_num,
                 int allocat,
                 int cellnum,
                 float ratio,
                 bool verbose) {
    unsigned long cprev=clock();
    if(!cgef_file.empty()) //从bgef生成cgef
    {
        BgefReader common_bin_gef = BgefReader(bgef_file, 1, true);
        ExpressionAttr expression_attr = common_bin_gef.getExpressionAttr();

        unsigned int mask_size[2]; // rows, cols
        mask_size[0] = expression_attr.max_y - expression_attr.min_y + 1;
        mask_size[1] = expression_attr.max_x - expression_attr.min_x + 1;

        Mask mask = Mask(mask_file, block_size, mask_size);
        if(verbose) cprev = printCpuTime(cprev, "Mask init");
        cout << "The number of cells (from mask file): " << mask.getCellNum() << endl;
        CgefWriter cgef_writer = CgefWriter(true);
        cgef_writer.setOutput(cgef_file);
        cgef_writer.setRandomCellTypeNum(rand_cell_type_num);
        cgef_writer.write(common_bin_gef, mask);
        cgef_writer.addLevel(allocat, cellnum, ratio, canvas_size, limit_blk);
    }
    else //为cgef 添加level层次
    {
        CgefWriter cgef_writer = CgefWriter(true);
        cgef_writer.setInput(bgef_file);
        cgef_writer.addLevel(allocat, cellnum, ratio, canvas_size, limit_blk);
    }

    if(verbose) printCpuTime(cprev, "generateCgef");
    return 0;
}