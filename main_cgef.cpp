// Created by huangzhibo on 2021/12/14.
//

#include "main_cgef.h"
#include "opencv2/opencv.hpp"
#include "cgefParam.h"
#include "cgefCellgem.h"
#include "cellAdjust.h"
#include "cgef_reader.h"

int ctest(const char *path)
{
    cellAdjust ca;
    ca.readBgef("/ldfssz1/ST_BI/USER/stereopy/test/tanliwei/test/test_data/FP200000443TL_E2.bgef");
    ca.readCgef("/ldfssz1/ST_BI/USER/stereopy/test/tanliwei/test/test_data/FP200000443TL_E2.cgef");
    vector<string> genename;
    vector<cellgem_label> vecCellgem;
    ca.getCellLabelgem(genename, vecCellgem);

    unordered_map<uint32_t, vector<DnbExpression>> map_dnb;
    DnbExpression tdnb;
    for(cellgem_label &cl : vecCellgem)
    {
        if(map_dnb.find(cl.cellid) == map_dnb.end())
        {
            vector<DnbExpression> tvec;
            map_dnb.emplace(cl.cellid, tvec);
        }
        tdnb.gene_id = cl.geneid;
        tdnb.count = cl.midcnt;
        tdnb.x = cl.x;
        tdnb.y = cl.y;
        map_dnb[cl.cellid].push_back(tdnb);
    }

    vector<Cell> veccell;
    vector<DnbExpression> vecdnb;
    auto itor = map_dnb.begin();
    uint32_t offset = 0;
    for(;itor != map_dnb.end();itor++)
    {
        Cell ce;
        ce.cellid = itor->first+1;
        ce.count = itor->second.size();
        ce.offset = offset;
        offset += ce.count;
        veccell.emplace_back(std::move(ce));
        vecdnb.insert(vecdnb.end(), itor->second.begin(), itor->second.end());
    }
    map_dnb.clear();

    ca.writeCellAdjust("jdkjsfl", (Cell*)veccell.data(), veccell.size(),
        (DnbExpression*)vecdnb.data(), vecdnb.size());
    return 0;
}

int ctest2(const char *path)
{
    // CgefReader cr(path);
    // vector<int> region={0,10000,0,10000};
    // vector<string> genelist={"AC240274.1","AC145212.1","KDM5D","TTTY14","LINC00278"};
    // vector<string> vec_gene;
    // vector<unsigned long long> uniq_cells;
    // vector<unsigned int> cell_ind;
    // vector<unsigned int> gene_ind;
    // vector<unsigned int> count;

    // cr.getfiltereddata(region, genelist, vec_gene, uniq_cells, cell_ind, gene_ind, count);



    return 0;
}

int split2int(char *str, const char *d, int *arry)
{
    int i = 0;
    char *p = strtok(str, d);
    while (p)
    {
        if(memcmp(p, "not exists", 10) == 0)
        {
            arry[i++]=USHRT_MAX;
        }
        else
        {
            arry[i++]=atoi(p);
        }
        p = strtok(NULL, d);
    }
    return i;
}

int change(const char *gef, const char *gem)
{
    std::map<uint32_t, uint16_t> map_cid_cluster;
    int cid_i = 0;
    int clu_i = 0;
    FILE *fgem = fopen(gem, "r");
    if(fgem)
    {
        char buf[128]={0};
        fgets(buf, 128, fgem);
        int count = 0;
        char *p = strtok(buf, "\t");
        while (p)
        {
            if(memcmp(p,"CellID", 6) == 0)
            {
                cid_i = count;
            }
            else if(memcmp(p,"Cellcluster", 11) == 0)
            {
                clu_i = count;
            }
            p = strtok(NULL, "\t");
            count++;
        }
        printf("%s %d %d\n", buf, cid_i, clu_i);

        int arry[10]; 
        while (fgets(buf, 128, fgem)!=NULL)
        {
            split2int(buf, "\t", arry);
            if(map_cid_cluster.find(arry[cid_i]) == map_cid_cluster.end())
            {
                map_cid_cluster.emplace(arry[cid_i], arry[clu_i]);
            }
        }
        fclose(fgem);
        printf("read gem end\n");
        hid_t fid = H5Fopen(gef, H5F_ACC_RDWR, H5P_DEFAULT);
        hid_t cell_dataset_id = H5Dopen(fid, "/cellBin/cell", H5P_DEFAULT);

        hsize_t dims[1];
        hid_t sid = H5Dget_space(cell_dataset_id);
        H5Sget_simple_extent_dims(sid, dims, nullptr);

        hid_t memtype = getMemtypeOfCellData();
        CellData *cell_array = (CellData *) malloc(dims[0] * sizeof(CellData));
        H5Dread(cell_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array);

        for(int i=0;i<dims[0];i++)
        {
            cell_array[i].cluster_id = map_cid_cluster[cell_array[i].id];
        }

        H5Dwrite(cell_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array);

        H5Sclose(sid);
        H5Tclose(memtype);
        H5Dclose(cell_dataset_id);
        H5Fclose(fid);
        free(cell_array);
    }
}

int cgef(int argc, char *argv[]) {
    //return change(argv[1], argv[2]);
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
    //("p,patch", "add the patch to cgef", cxxopts::value<int>()->default_value("0"))
    ("n,cnum", "top level cell num", cxxopts::value<int>()->default_value("5000"), "INT")
    ("R,ratio", "other level cell num ratio", cxxopts::value<int>()->default_value("20"), "FLOAT")
    ("a,allocat", "allocation strategy ", cxxopts::value<int>()->default_value("2"), "INT")
    ("g,raw-gem", "raw gem file", cxxopts::value<std::string>(), "FILE")
    ("c,canvas", "set canvas size, minx,miny,maxx,maxy", cxxopts::value<std::string>()->default_value("0,0,90000,90000"), "FILE")
    ("l,limit", "set blk limit", cxxopts::value<std::string>()->default_value("16,16"), "FILE")
    ("S,split", "split cellid to layers and blks", cxxopts::value<int>()->default_value("0"))
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

    if (result.count("mask-file") != 1){
        cgefParam::GetInstance()->m_maskstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_maskstr = result["mask-file"].as<string>();
    }

    if (result.count("output-file") != 1){
        cgefParam::GetInstance()->m_outputstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_outputstr = result["output-file"].as<string>();
    }

    if(result.count("raw-gem") != 1)
    {
        cgefParam::GetInstance()->m_rawgemstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_rawgemstr = result["raw-gem"].as<string>();
    }

    int rand_celltype_num = result["rand-celltype"].as<int>();
    cgefParam::GetInstance()->m_inputstr = result["input-file"].as<string>();
    cgefParam::GetInstance()->m_threadcnt = result["threads"].as<int>();
    
    vector<string> block_size_tmp = split(result["block"].as<string>(), ',');
    if(block_size_tmp.size() != 2){
        std::cerr << "[ERROR] The -b,--block parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        exit(1);
    }
    cgefParam::GetInstance()->m_block_size[0] = static_cast<int>(strtol(block_size_tmp[0].c_str(), nullptr, 10));
    cgefParam::GetInstance()->m_block_size[1] = static_cast<int>(strtol(block_size_tmp[1].c_str(), nullptr, 10));

    //分层分块参数-----
    int allocat = 2;
    int topcellnum = 5000;
    float ratio = 0.2;
    int canvas_size[4]={0,0,0,0}; //全局画布大小
    int limit_blk[2] = {0,0};//分块限制
    int isplit = result["split"].as<int>();
    if(isplit > 0)
    {
        int tmpr = result["ratio"].as<int>();
        ratio = tmpr*1.0/100;
        topcellnum = result["cnum"].as<int>();
        allocat = result["allocat"].as<int>();
        vector<string> canvas_size_tmp = split(result["canvas"].as<string>(), ',');
        vector<string> limit_tmp = split(result["limit"].as<string>(), ',');

        canvas_size[0] = static_cast<int>(strtol(canvas_size_tmp[0].c_str(), nullptr, 10));
        canvas_size[1] = static_cast<int>(strtol(canvas_size_tmp[1].c_str(), nullptr, 10));
        canvas_size[2] = static_cast<int>(strtol(canvas_size_tmp[2].c_str(), nullptr, 10));
        canvas_size[3] = static_cast<int>(strtol(canvas_size_tmp[3].c_str(), nullptr, 10));

        limit_blk[0] = static_cast<int>(strtol(limit_tmp[0].c_str(), nullptr, 10));
        limit_blk[1] = static_cast<int>(strtol(limit_tmp[1].c_str(), nullptr, 10));

        CgefWriter cgef_writer(true);
        cgef_writer.setInput(cgefParam::GetInstance()->m_inputstr);
        if(isplit == 1)
        {
            cgef_writer.addLevel_1();
        }
        else
        {
            cgef_writer.addLevel(allocat, topcellnum, ratio, canvas_size, limit_blk);
        }
        
        return 0;
    }
    //-------------------
    
    generateCgef(cgefParam::GetInstance()->m_outputstr,
                cgefParam::GetInstance()->m_inputstr, 
                cgefParam::GetInstance()->m_maskstr,
                cgefParam::GetInstance()->m_rawgemstr,
                cgefParam::GetInstance()->m_block_size,
                rand_celltype_num);
    return 0;
}

int generateCgef(const string &cgef_file,
                 const string &bgef_file,
                 const string &mask_file,
                 const string& raw_gem,
                 const int* block_size,
                 int rand_celltype_num,
                 bool verbose) {
    unsigned long cprev=clock();
    CgefWriter cgef_writer(verbose);
    cgef_writer.setOutput(cgef_file);
    cgef_writer.setRandomCellTypeNum(rand_celltype_num);

    cgefCellgem cgem;
    cgem.writeFile(&cgef_writer, mask_file, bgef_file, raw_gem);

    if(verbose) printCpuTime(cprev, "generateCgef");
    return 0;
}