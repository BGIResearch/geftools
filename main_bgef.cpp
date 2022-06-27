//
// Created by huangzhibo on 2021/12/14.
//

#include "main_bgef.h"

#include <utility>
#include "thread_pool.h"
#include "read_task.h"
#include "dnb_merge_task.h"
#include "bin_task.h"
#include "special_bin.h"
#include "utils.h"
#include "timer.h"

int test(const char *path)
{
    timer st1("");
    BgefReader bgef_reader(path, 200, 1);
    bgef_reader.getReduceExpression();

    printf("end\n");

    return 0;
}

int bgef(int argc, char *argv[]) {
    //return test(argv[1]);
    cxxopts::Options options("geftools bgef",
                       "About:  Generate common bin GEF(.bgef) according to gem file or bin1 GEF\n");
    options
    .set_width(120)
    .add_options()
    ("i,input-file", "input gene expression matrix file(.gem/.gem.gz) or bin1 bGEF file [request]",
            cxxopts::value<std::string>(), "FILE")
    ("o,output-file", "output bin GEF file (.bgef) [request]", cxxopts::value<std::string>(), "FILE")
    ("b,bin-size", "Set bin size by the comma-separated list [request]",
            cxxopts::value<std::string>()->default_value("1,10,20,50,100,200,500"), "STR")
    ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list"
                 " of two vertex coordinates (minX,maxX,minY,maxY)",
                 cxxopts::value<std::string>()->default_value(""), "STR")
    ("t,threads", "number of threads", cxxopts::value<int>()->default_value("8"), "INT")
    ("s,stat", "create stat group", cxxopts::value<bool>()->default_value("true"))
    ("O,omics", "input omics [request]", cxxopts::value<std::string>()->default_value("Transcriptomics"), "STR")
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (argc <= 1 || result.count("help"))
    {
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

    if (result.count("omics") != 1){
        std::cout << "[ERROR] The -O,--omics parameter must be given correctly.\n" << std::endl;
        std::cout << options.help() << std::endl;
        exit(1);
    }

    BgefOptions *opts = BgefOptions::GetInstance();

    opts->input_file_ = result["input-file"].as<string>();
    opts->output_file_ = result["output-file"].as<string>();
    bool bstat= result["stat"].as<bool>();

    vector<string> bs_tmp = split(result["bin-size"].as<string>(), ',');

    for(auto & binsize : bs_tmp){
        opts->bin_sizes_.emplace_back(static_cast<unsigned int>(strtol(binsize.c_str(), nullptr, 10)));
    }

    bool b100 = false;
    for(int bin : opts->bin_sizes_)
    {
        if(bin == 100)
        {
            b100 = true; //用户指定了bin100
            opts->m_stattype = 2;
            break;
        }
    }
    if(!b100 && bstat) //用户没有指定bin100，但是需要生成stat
    {
        opts->bin_sizes_.emplace_back(100);
        opts->m_stattype = 1;
    }

    if (result.count("region") == 1){
        string region_tmp = result["region"].as<string>();
        vector<string> regions = split(region_tmp, ',');
        for(auto & r : regions){
            opts->region_.emplace_back(static_cast<unsigned int>(strtol(r.c_str(), nullptr, 10)));
        }
    }

    opts->thread_ = result["threads"].as<int>();
    opts->verbose_ = result["verbose"].as<bool>();
    opts->m_stromics = result["omics"].as<string>();

    gem2gef(opts);
//    generateBgef(opts.output_file, opts.input_file, opts.verbose);
    return 0;
}

int generateBgef(const string &input_file,
                 const string &bgef_file,
                 const string &stromics,
                 int n_thread,
                 vector<unsigned int> bin_sizes,
                 vector<int> region,
                 bool verbose,
                 bool bstat) {
    unsigned long cprev=clock();
    BgefOptions *opts = BgefOptions::GetInstance();
    opts->input_file_ = input_file;
    opts->output_file_ = bgef_file;
    opts->bin_sizes_ = std::move(bin_sizes);
    opts->region_ = std::move(region);
    opts->thread_ = n_thread;
    opts->verbose_ = verbose;
    opts->m_stromics = stromics;

    bool b100 = false;
    for(int bin : opts->bin_sizes_)
    {
        if(bin == 100)
        {
            b100 = true; //用户指定了bin100
            opts->m_stattype = 2;
            break;
        }
    }
    if(!b100 && bstat) //用户没有指定bin100，但是需要生成stat
    {
        opts->bin_sizes_.emplace_back(100);
        opts->m_stattype = 1;
    }
    
    gem2gef(opts);
    if(verbose) printCpuTime(cprev, "generateBgef");
    return 0;
}

void gem2gef(BgefOptions *opts)
{
    unsigned long cprev0 =clock(), cprev;

    unsigned int resolution;
    if(!H5Fis_hdf5(opts->input_file_.c_str())){
        mRead(opts);
        resolution = parseResolutin(opts->input_file_);
    }else{
        BgefReader bgef_reader(opts->input_file_, 1, opts->verbose_);
        ExpressionAttr expression_attr = bgef_reader.getExpressionAttr();

        if(opts->region_.empty()){
            bgef_reader.getGeneExpression(opts->map_gene_exp_);
            opts->m_bexon = bgef_reader.isExonExist();
            opts->range_ = {expression_attr.min_x, expression_attr.max_x, expression_attr.min_y, expression_attr.max_y};
            opts->offset_x_ = expression_attr.min_x;
            opts->offset_y_ = expression_attr.min_y;
        }else{
            bgef_reader.getGeneExpression(opts->map_gene_exp_, opts->region_);
            opts->m_bexon = bgef_reader.isExonExist();
            int min_x = opts->region_[0];
            int max_x = opts->region_[1];
            int min_y = opts->region_[2];
            int max_y = opts->region_[3];

            opts->range_ = {expression_attr.min_x + min_x, min(expression_attr.max_x, max_x+expression_attr.min_x),
                            expression_attr.min_y + min_y, min(expression_attr.max_y, max_y+expression_attr.min_y)};
            opts->offset_x_ = expression_attr.min_x + min_x;
            opts->offset_y_ = expression_attr.min_y + min_y;
        }

        resolution = expression_attr.resolution;
    }

    if(opts->verbose_) cprev = printCpuTime(cprev0, "read gene expression file");
    if(opts->map_gene_exp_.empty())
    {
        printf("the exp is empty\n");
        return;
    }

    opts->m_genes_queue.init(opts->map_gene_exp_.size());
    ThreadPool thpool(opts->thread_ * 2);

    BgefWriter bgef_writer(opts->output_file_, opts->verbose_, opts->m_bexon, opts->m_stromics);
    bgef_writer.setResolution(resolution);
    
    int genecnt = 0;
    for(unsigned int bin : opts->bin_sizes_)
    {
        cprev=clock();
        auto& dnb_matrix = opts->dnbmatrix_;
        auto& dnbAttr = opts->dnbmatrix_.dnb_attr;
        auto& range = opts->range_;

        dnbAttr.min_x = (opts->offset_x_ / bin) * bin;
        dnbAttr.len_x = (range[1] - range[0])/bin + 1;
        dnbAttr.min_y = (opts->offset_y_ / bin) * bin;
        dnbAttr.len_y = (range[3]-range[2])/bin + 1;
        dnbAttr.max_gene = 0;
        dnbAttr.max_mid = 0;
        dnbAttr.number = 0;
        unsigned long matrix_len = (unsigned long)(dnbAttr.len_x) * dnbAttr.len_y;
        printf("bin %d matrix: min_x=%d len_x=%d min_y=%d len_y=%d matrix_len=%lu\n",
               bin, dnbAttr.min_x, dnbAttr.len_x, dnbAttr.min_y, dnbAttr.len_y, matrix_len);
        if (bin == 1)
        {
            dnb_matrix.pmatrix_us = (BinStatUS*)calloc(matrix_len, sizeof(BinStatUS));
            assert(dnb_matrix.pmatrix_us);
            if(opts->m_bexon)
            {
                dnb_matrix.pexon16 = (unsigned short*)calloc(matrix_len, 2);
                assert(dnb_matrix.pexon16);
            }
        }
        else
        {
            dnb_matrix.pmatrix = (BinStat*)calloc(matrix_len, sizeof(BinStat));
            assert(dnb_matrix.pmatrix);
            if(opts->m_bexon)
            {
                dnb_matrix.pexon32 = (unsigned int*)calloc(matrix_len, 4);
                assert(dnb_matrix.pexon32);
            }
        }


        for(int i=0; i < opts->thread_; i++)
        {
            auto *task = new DnbMergeTask(opts->map_gene_exp_.size(), i, bin);
            thpool.addTask(task);
        }

        auto itor = opts->map_gene_exp_.begin();
        for(;itor != opts->map_gene_exp_.end();itor++)
        {
            auto *task = new BinTask(bin, itor->first.c_str());
            thpool.addTask(task);
        }

        if (bin == 1)
        {
            unsigned long totalSize = 0;
            for (auto& p : opts->map_gene_exp_)
                totalSize += p.second.size();
            opts->expressions_.reserve(totalSize);
            opts->genes_.reserve(opts->map_gene_exp_.size());
        }

        unsigned int offset = 0;
        unsigned int maxexp = 0;
        unsigned int maxexon = 0;
        genecnt = 0;
        while (true) //write gene
        {
            GeneInfo *pgeneinfo = opts->m_geneinfo_queue.getPtr();
            if (bin == 1){
                opts->expressions_.insert(opts->expressions_.end(), pgeneinfo->vecptr->begin(), pgeneinfo->vecptr->end());
            }
            else
            {
                if(bin != 100 || opts->m_stattype == 2)
                {
                    for (auto g : *pgeneinfo->vecptr)
                    {
                        g.x *= bin;
                        g.y *= bin;
                        opts->expressions_.push_back(std::move(g));
                    }
                }
            }

            if(bin != 100 || opts->m_stattype == 2)
            {
                opts->genes_.emplace_back(pgeneinfo->geneid, offset, static_cast<unsigned int>(pgeneinfo->vecptr->size()));
                offset += pgeneinfo->vecptr->size();
                maxexp = std::max(maxexp, pgeneinfo->maxexp);
                maxexon = std::max(maxexon, pgeneinfo->maxexon);
            }

            if(bin == 100)
            {
                GeneErank erank(pgeneinfo->geneid);
                erank.umicnt = pgeneinfo->umicnt;
                erank.e10 = pgeneinfo->e10;
                //erank.c50 = pgeneinfo->c50;
                opts->total_umicnt_ += pgeneinfo->umicnt;
                opts->vec_bin100_.emplace_back(erank);
            }
            delete pgeneinfo;
            genecnt++;
            if(genecnt == opts->map_gene_exp_.size())
            {
                break;
            }
        }

        if(bin != 100 || opts->m_stattype == 2)
        {
            bgef_writer.storeGene(opts->expressions_, opts->genes_, dnb_matrix.dnb_attr, maxexp, bin);
            bgef_writer.storeGeneExon(opts->expressions_, maxexon, bin);
            opts->expressions_.clear();
            opts->genes_.clear();
        }

        thpool.waitTaskDone();
        opts->m_genes_queue.clear(bin);
        //write dnb
        writednb(opts, bgef_writer, bin);


        if (bin == 1)
        {
            if (dnb_matrix.pmatrix_us != nullptr)
            {
                free(dnb_matrix.pmatrix_us);
                dnb_matrix.pmatrix_us = nullptr;
                if(opts->m_bexon)
                {
                    free(dnb_matrix.pexon16);
                    dnb_matrix.pexon16 = nullptr;
                }
            }
        }
        else
        {
            if (dnb_matrix.pmatrix != nullptr)
            {
                free(dnb_matrix.pmatrix);
                dnb_matrix.pmatrix = nullptr;
                if(opts->m_bexon)
                {
                    free(dnb_matrix.pexon32);
                    dnb_matrix.pexon32 = nullptr;
                }
            }
        }
        if(opts->verbose_) printCpuTime(cprev, "bin process");
    }

    if(opts->verbose_) printCpuTime(cprev0, "gem2gef");
}


int mRead(BgefOptions *opts) //多线程读
{
    opts->infile_ = gzopen(opts->input_file_.c_str(), "r");
    gzbuffer(opts->infile_, READLEN);

    // Process the header lines
    std::string line;
    while (readline(opts->infile_, line))
    {
        if (line[0] == '#')
        {
            // Skip the offset parameter
            if (line.substr(0, 9) == "#OffsetX=")
                opts->offset_x_ = stoi(line.substr(9));
            else if (line.substr(0, 9) == "#OffsetY=")
                opts->offset_y_ = stoi(line.substr(9));
            continue;
        }
        if (line.substr(0, 6) == "geneID") break;
    }

    ThreadPool thpool(opts->thread_);
    for(int i=0;i<opts->thread_;i++)
    {
        auto *rtask = new ReadTask();
        thpool.addTask(rtask);
    }

    while (true)
    {
        sleep(1);
        if(thpool.idlCount() == opts->thread_)
        {
            break;
        }
    }
    gzclose(opts->infile_);

    // Subtract min value of coordinates
    int minx = opts->range_[0];
    int miny = opts->range_[2];
    if (minx != 0 || miny != 0)
    {
        opts->offset_x_ += minx;
        opts->offset_y_ += miny;
        for (auto& p : opts->map_gene_exp_)
        {
            for (auto& g : p.second)
            {
                g.x -= minx;
                g.y -= miny;
            }
        }
    }

    return 0;
}

unsigned int parseResolutin(string& filename) {
    std::unordered_map<string, unsigned int> pitch({
        {"CL1", 900},{"N1", 900},{"V3", 715},{"K2", 715},{"S2", 715},
        {"S1", 900},{"F3", 715},{"F1", 800},{"V1", 800},{"DP84", 715},
        {"DP8", 850},{"FP2", 500},{"SS2", 500},{"FP1", 600},{"E1", 700},{"DP40", 700},
        {"G1", 700},{"A", 500},{"B", 500},{"C", 500},{"D", 500},
        {"U", 715},{"V", 715},{"W", 715},{"X", 715},{"Y", 500}});

    auto pos = filename.find_last_of('/');
    if (pos == std::string::npos)
        pos = -1;

    unsigned int result = 0;
    string chip_prefix = filename.substr(pos+1, 4);
    while (chip_prefix.size() >= 2){
        if (pitch.count(chip_prefix) != 0)
        {
            result = pitch[chip_prefix];
            break;
        }
        chip_prefix.pop_back();
    }

    return result;
}

// sort gene name by MIDCount
void sortGeneByCnt(std::unordered_map <std::string, std::vector<Expression>>& data,
                   vector<pair<string, unsigned int>>& geneCnts)
{
    auto itor = data.begin();
    unsigned int umicnt = 0;
    for(;itor != data.end();itor++)
    {
        umicnt = 0;
        for(Expression &exp : itor->second)
        {
            umicnt += exp.count;
        }
        geneCnts.emplace_back(std::make_pair(itor->first, umicnt));
    }

    typedef pair<string, unsigned int> MyPair;
    std::sort(geneCnts.begin(), geneCnts.end(), [](const MyPair& p1, const MyPair& p2){
        if (p1.second > p2.second)
            return true;
        else if (p1.second == p2.second)
            return p1.first < p2.first;
        else
            return false;
    });
}

unsigned int Boxplot(vector<unsigned int> &vecmid)
{
    sort(vecmid.begin(), vecmid.end(), 
        [](const unsigned int a, const unsigned int b){return a<b;});
    int sz = vecmid.size();
    int pos_q1 = ceil( (sz + 1)*1.0 /4 );
    float q1 = vecmid[pos_q1-2]*0.25+vecmid[pos_q1-1]*0.75; //下四分位数Q1

    int pos_q2 = ceil( (sz + 1)*2.0 /4 );
    float q2 = vecmid[pos_q2-2]*0.5+vecmid[pos_q2-1]*0.5; //中位数

    int pos_q3 = ceil( (sz + 1)*3.0 /4 );
    float q3 = vecmid[pos_q3-2]*0.75+vecmid[pos_q3-1]*0.25; //上四分位数Q3

    float IQR = q3 - q1; //四分位距
    unsigned int uplimit = ceil(q3+1.5*IQR);
    //float lowlimit = q1-1.5*IQR;

    sz--;
    return vecmid[sz] < uplimit ? vecmid[sz] : uplimit;
}

void writednb(BgefOptions *opts, BgefWriter &bgef_writer, int bin)
{
    unsigned long cprev =clock();
    if(bin == 100)
    {
        vector<pair<string, unsigned int>> geneCnts;
        sortGeneByCnt(opts->map_gene_exp_, geneCnts);
        
        std::vector<float> vec_e10_result;
        SpecialBin sbin;
        sbin.calcE10(geneCnts, vec_e10_result);

        vector<GeneStat> geneStat;
        size_t sz = geneCnts.size();
        geneStat.reserve(sz);
        for (int i = 0; i < sz; ++i){
            geneStat.emplace_back(geneCnts[i].first, geneCnts[i].second, vec_e10_result[i]);
        }
        bgef_writer.storeStat(geneStat);
        if(opts->m_stattype != 2)
        {
            return;
        }
    }

    vector<unsigned int> vec_mid;
    DnbMatrix &dnbM = opts->dnbmatrix_;
    unsigned long number = 0;
    unsigned long matrix_len = (unsigned long)(dnbM.dnb_attr.len_x) * dnbM.dnb_attr.len_y;
    if (bin == 1)
    {
        for(unsigned long i=0;i<matrix_len;i++)
        {
            if(dnbM.pmatrix_us[i].gene_count)
            {
                ++number;
                vec_mid.push_back(dnbM.pmatrix_us[i].mid_count);
            }
        }
    }
    else
    {
        for(unsigned long i=0;i<matrix_len;i++)
        {
            if(dnbM.pmatrix[i].gene_count)
            {
                ++number;
                vec_mid.push_back(dnbM.pmatrix[i].mid_count);
            }
        }
    }

    int sz = vec_mid.size();
    sort(vec_mid.begin(), vec_mid.end(), [](const unsigned int a, const unsigned int b){return a<b;});
    if(bin > 50)
    {
        //printf("%d\n", vec_mid[sz-1]);
        dnbM.dnb_attr.max_mid = vec_mid[sz-1];
    }
    else
    {
        int limit = sz*0.999;
        //printf("%d %d %d %d \n", sz, limit, vec_mid[sz-1], vec_mid[limit]);
        dnbM.dnb_attr.max_mid = vec_mid[limit];
    }
    
    dnbM.dnb_attr.number = number;
    bgef_writer.storeDnb(dnbM, bin);
    bgef_writer.storeWholeExon(dnbM, bin);

    if(opts->verbose_) printCpuTime(cprev, "writednb");
}



void StereoDataToGef(const string &output_file, int binsize, int sz, unsigned long *cellptr)
{
    for(int i=0;i<sz;i++)
    {
        printf("%ld\n", cellptr[i]);
    }
}
