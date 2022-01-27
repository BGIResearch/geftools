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

int bgef(int argc, char *argv[]) {

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

    BgefOptions *opts = BgefOptions::GetInstance();

//    BgefOptions opts = {
//            result["input-file"].as<string>(),
//            result["output-file"].as<string>(),
//    };

    opts->input_file_ = result["input-file"].as<string>();
    opts->output_file_ = result["output-file"].as<string>();

    vector<string> bs_tmp = split(result["bin-size"].as<string>(), ',');

    for(auto & binsize : bs_tmp){
        opts->bin_sizes_.emplace_back(static_cast<unsigned int>(strtol(binsize.c_str(), nullptr, 10)));
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

    gem2gef(opts);
//    generateBgef(opts.output_file, opts.input_file, opts.verbose);
    return 0;
}

int generateBgef(const string &input_file,
                 const string &bgef_file,
                 int n_thread,
                 vector<unsigned int> bin_sizes,
                 vector<unsigned int> region,
                 bool verbose) {
    unsigned long cprev=clock();
    BgefOptions *opts = BgefOptions::GetInstance();
    opts->input_file_ = input_file;
    opts->output_file_ = bgef_file;
    opts->bin_sizes_ = std::move(bin_sizes);
    opts->region_ = std::move(region);
    opts->thread_ = n_thread;
    opts->verbose_ = verbose;
    gem2gef(opts);
    if(verbose) printCpuTime(cprev, "generateBgef");
    return 0;
}

void gem2gef(BgefOptions *opts)
{
    unsigned long cprev0 =clock(), cprev;

    unsigned int resolution;
    if(decideSuffix(opts->input_file_, "gem") || decideSuffix(opts->input_file_, "gz")){
        mRead(opts);
        resolution = parseResolutin(opts->input_file_);
    }else{
        BgefReader bgef_reader(opts->input_file_, 1, opts->verbose_);
        ExpressionAttr expression_attr = bgef_reader.getExpressionAttr();

        if(opts->region_.empty()){
            bgef_reader.getGeneExpression(opts->map_gene_exp_);
            opts->range_ = {expression_attr.min_x, expression_attr.max_x, expression_attr.min_y, expression_attr.max_y};
            opts->offset_x_ = expression_attr.min_x;
            opts->offset_y_ = expression_attr.min_y;
        }else{
            bgef_reader.getGeneExpression(opts->map_gene_exp_, opts->region_);
            unsigned int min_x = opts->region_[0];
            unsigned int max_x = opts->region_[1];
            unsigned int min_y = opts->region_[2];
            unsigned int max_y = opts->region_[3];

            opts->range_ = {expression_attr.min_x + min_x, min(expression_attr.max_x, max_x+expression_attr.min_x),
                            expression_attr.min_y + min_y, min(expression_attr.max_y, max_y+expression_attr.min_y)};
            opts->offset_x_ = expression_attr.min_x + min_x;
            opts->offset_y_ = expression_attr.min_y + min_y;
        }

        resolution = expression_attr.resolution;
    }

    if(opts->verbose_) cprev = printCpuTime(cprev0, "read gene expression file");

    opts->gene_info_queue_.init(opts->map_gene_exp_.size());
    ThreadPool thpool(opts->thread_ * 2);

    BgefWriter bgef_writer(opts->output_file_, opts->verbose_);
    bgef_writer.setResolution(resolution);

    int genecnt = 0;
    for(unsigned int bin : opts->bin_sizes_)
    {
        cprev=clock();
        auto& dnb_matrix = opts->dnbmatrix_;
        auto& dnbAttr = opts->dnbmatrix_.dnb_attr;
        auto& range = opts->range_;

        dnbAttr.min_x = (opts->offset_x_ / bin) * bin;
        dnbAttr.len_x = int((float(range[1]) / bin) - (float(range[0]) / bin)) + 1;
        dnbAttr.min_y = (opts->offset_y_ / bin) * bin;
        dnbAttr.len_y = int((float(range[3]) / bin) - (float(range[2]) / bin)) + 1;
        unsigned long matrix_len = (unsigned long)(dnbAttr.len_x) * dnbAttr.len_y;
        if (bin == 1)
        {
            dnb_matrix.pmatrix_us = (BinStatUS*)calloc(matrix_len, sizeof(BinStatUS));
            assert(dnb_matrix.pmatrix_us);
        }
        else
        {
            dnb_matrix.pmatrix = (BinStat*)calloc(matrix_len, sizeof(BinStat));
            assert(dnb_matrix.pmatrix);
        }
        printf("bin %d matrix: min_x=%d len_x=%d min_y=%d len_y=%d matrix_len=%lu\n",
               bin,
               dnbAttr.min_x, dnbAttr.len_x,
               dnbAttr.min_y, dnbAttr.len_y,
               matrix_len);

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

        // timer time_write_gene("write gene");
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
        genecnt = 0;
        while (true) //write gene
        {
            GeneInfo2 *pgenedata = opts->gene_queue_.getGeneInfo2();
            if (bin == 1){
                opts->expressions_.insert(opts->expressions_.end(), pgenedata->vecdataptr->begin(), pgenedata->vecdataptr->end());
            }
            else
            {
                for (auto g : *pgenedata->vecdataptr)
                {
                    g.x *= bin;
                    g.y *= bin;
                    opts->expressions_.push_back(std::move(g));
                }
            }

            opts->genes_.emplace_back(pgenedata->geneid, offset, static_cast<unsigned int>(pgenedata->vecdataptr->size()));
            offset += pgenedata->vecdataptr->size();
            maxexp = std::max(maxexp, pgenedata->maxexp);
            // h5Writer.storeGene(pgenedata->geneid, *(pgenedata->vecdataptr));

            if(bin == 100)
            {
                GeneErank erank(pgenedata->geneid);
                erank.umicnt = pgenedata->umicnt;
                erank.e10 = pgenedata->e10;
                erank.c50 = pgenedata->c50;
                opts->total_umicnt_ += pgenedata->umicnt;
                opts->vec_bin100_.emplace_back(erank);
            }

            // delete pgenedata->vecdataptr;
            // delete pgenedata;

            genecnt++;
            if(genecnt == opts->map_gene_exp_.size())
            {
                break;
            }
        }

        bgef_writer.storeGene(opts->expressions_, opts->genes_, dnb_matrix.dnb_attr, maxexp, bin);
        opts->expressions_.clear();
        opts->genes_.clear();

//        cprev = printCpuTime(cprev, "wait");
        thpool.waitTaskDone();
//        printCpuTime(cprev, "waitTaskDone");
        // tm.stop();
        opts->gene_info_queue_.clear(bin);
        //write dnb
        writednb(opts, bgef_writer, bin);

        if (bin == 1)
        {
            if (dnb_matrix.pmatrix_us != nullptr)
            {
                free(dnb_matrix.pmatrix_us);
                dnb_matrix.pmatrix_us = nullptr;
            }
        }
        else
        {
            if (dnb_matrix.pmatrix != nullptr)
            {
                free(dnb_matrix.pmatrix);
                dnb_matrix.pmatrix = nullptr;
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
    unsigned int minx = opts->range_[0];
    unsigned int miny = opts->range_[2];
    if (opts->offset_x_ == 0 && opts->offset_y_ == 0)
    {
        opts->offset_x_ = minx;
        opts->offset_y_ = miny;
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
        {"DP8", 850},{"FP2", 500},{"FP1", 600},{"E1", 700},{"DP40", 700},
        {"G1", 700},{"A", 500},{"B", 500},{"C", 500},{"D", 500},
        {"U", 715},{"V", 715},{"W", 715},{"X", 715},{"Z", 500},
        {"Y", 900}});

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

void writednb(BgefOptions *opts, BgefWriter &bgef_writer, int bin)
{
    unsigned long cprev =clock();
    if(bin == 100)
    {
        // h5Writer.storeDnbRange(pcmd->m_range);
        vector<pair<string, unsigned int>> geneCnts;
        sortGeneByCnt(opts->map_gene_exp_, geneCnts);
        // h5Writer.storeGeneName(geneCnts);

        std::vector<float> vec_e10_result;
        SpecialBin sbin;
        sbin.calcE10(geneCnts, vec_e10_result);
        // h5Writer.storeE10(vec_e10_result);

        vector<GeneStat> geneStat;
        size_t sz = geneCnts.size();
        geneStat.reserve(sz);
        for (int i = 0; i < sz; ++i){
            geneStat.emplace_back(geneCnts[i].first, geneCnts[i].second, vec_e10_result[i]);
        }
        bgef_writer.storeStat(geneStat);
    }

    DnbMatrix &dnbM = opts->dnbmatrix_;
    unsigned long number = 0;
    unsigned long matrix_len = (unsigned long)(dnbM.dnb_attr.len_x) * dnbM.dnb_attr.len_y;
    if (bin == 1)
    {
        for(unsigned long i=0;i<matrix_len;i++)
            if(dnbM.pmatrix_us[i].gene_count)
                ++number;
    }
    else
    {
        for(unsigned long i=0;i<matrix_len;i++)
            if(dnbM.pmatrix[i].gene_count)
                ++number;
    }
    dnbM.dnb_attr.number = number;
    bgef_writer.storeDnb(dnbM, bin);

//    if(bin == 200)
//    {
//        special_bin sbin;
//        std::vector<int> vecdnb;
//        unsigned int x, y;
//        unsigned int y_len = dnbM.dnb_attr.len_y;
//        for(unsigned long i=0;i<matrix_len;i++)
//        {
//            if(dnbM.pmatrix[i].gene_count)
//            {
//                x = i/y_len;
//                y = i%y_len;
//
//                vecdnb.push_back(x*bin);
//                vecdnb.push_back(y*bin);
//                vecdnb.push_back(dnbM.pmatrix[i].mid_count);
//            }
//        }
//        sbin.createPNG_py(vecdnb);
//    }
    if(opts->verbose_) printCpuTime(cprev, "writednb");
}

