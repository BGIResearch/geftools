
#ifndef GENETOH5_COMMANDPARSE_H
#define GENETOH5_COMMANDPARSE_H

#include <unordered_map>
#include <zlib.h>
#include "gef.h"
#include "gene_info_queue.h"
#include "gene_queue.h"


class BgefOptions {
private:
    BgefOptions(){};
    ~BgefOptions(){};
public:
    static BgefOptions *GetInstance()
    {
        static BgefOptions instance;
        return &instance;
    }

    int thread_ = 8; //设置取bin线程数
//    int m_thread_dnb = 8; //设置dnbmerge线程数
    bool reverse_ = false; // true: gef to gem, false: gem to gef
    bool verbose_ = false;

    string input_file_;
    string output_file_;
    vector<int> bin_sizes_;

    std::unordered_map <std::string, std::vector<Expression>> map_gene_exp_;
    std::vector<GeneErank> vec_bin100_;
    unsigned long total_umicnt_ = 0;

    gzFile infile_;

    mutex dnbmtx_;
    DnbMatrix dnbmatrix_;
    std::vector<unsigned int> range_ = {UINT_MAX, 0, UINT_MAX, 0};
    GeneInfoQueue gene_info_queue_;
    GeneQueue gene_queue_;
    std::vector<Expression> expressions_;
    std::vector<Gene> genes_;
    unsigned int offset_x_ = 0; // offset of coordinate
    unsigned int offset_y_ = 0;
};

#endif //GENETOH5_COMMANDPARSE_H