#include "bin_task.h"


BinTask::BinTask(int bin, const char *ptr):
    m_bin(bin),
    m_geneid(ptr),
    m_maxexp(0)
{
    opts_ = BgefOptions::GetInstance();
}

BinTask::~BinTask()
{
}

void BinTask::bin1task()
{
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];

    auto *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = &exp_vec;
    opts_->gene_info_queue_.addqueue(pgeneinfo);

    auto *pgenedata = new GeneInfo2(m_geneid);
    pgenedata->vecdataptr = &exp_vec;
    for (auto& exp : exp_vec)
    {
        if (exp.count > m_maxexp)
            m_maxexp = exp.count;
    }
    pgenedata->maxexp = m_maxexp;
    // pgenedata->vecdataptr = new std::vector<int>;
    // pgenedata->vecdataptr->reserve(vecgene.size()*3);
    // for(auto gene : vecgene)
    // {
    //     pgenedata->vecdataptr->push_back(gene.x);
    //     pgenedata->vecdataptr->push_back(gene.y);
    //     pgenedata->vecdataptr->push_back(gene.cnt);
    // }
    opts_->gene_queue_.addqueue(pgenedata);
//    printf("bin1task\n");
}

void BinTask::bin100task()
{
    unsigned long umicnt = 0;
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];
    for(auto exp : exp_vec)
    {
        x = exp.x / m_bin;
        y = exp.y / m_bin;
        dnb = x<<32 | y;
        map_dnb[dnb]+=exp.count;
        umicnt += exp.count;
    }

    auto *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = new std::vector<Expression>;
    pgeneinfo->vecptr->reserve(map_dnb.size());

    auto *pgenedata = new GeneInfo2(m_geneid);
    // pgenedata->vecdataptr = new std::vector<int>;
    // pgenedata->vecdataptr->reserve(map_dnb.size()*3);
    pgenedata->umicnt = umicnt;
    auto itor_dnb = map_dnb.begin();
    Expression exp{0,0,0};
    for(;itor_dnb!=map_dnb.end();itor_dnb++)
    {
        exp.x = itor_dnb->first>>32;
        exp.y = itor_dnb->first & 0xFFFFFFFF;
        exp.count = itor_dnb->second;
        pgeneinfo->vecptr->emplace_back(exp);

        // pgenedata->vecdataptr->push_back(gene.x);
        // pgenedata->vecdataptr->push_back(gene.y);
        // pgenedata->vecdataptr->push_back(gene.cnt);
        if (exp.count > m_maxexp)
            m_maxexp = exp.count;
    }
    pgenedata->maxexp = m_maxexp;
    pgenedata->vecdataptr = pgeneinfo->vecptr;

    std::sort(pgeneinfo->vecptr->begin(), pgeneinfo->vecptr->end(), 
            [](const Expression &a, const Expression &b){return a.count > b.count;});

    int j = 0;
    int sz = pgeneinfo->vecptr->size() * 0.1;
    unsigned long midcnt = 0;
    auto itor = pgeneinfo->vecptr->begin();
    for(;itor != pgeneinfo->vecptr->end() && j< sz;itor++,j++)
    {
        midcnt += itor->count;
    }
    pgenedata->e10 = (midcnt*1.0/umicnt)*100;

    itor = pgeneinfo->vecptr->begin();
    midcnt = 0;
    j = 0;
    for(;itor != pgeneinfo->vecptr->end();itor++,j++)
    {
        midcnt += itor->count;
        if(midcnt*1.0/umicnt > 0.5) 
        {
            break;
        }
    }

    sz = pgeneinfo->vecptr->size();
    pgenedata->c50 = (j*1.0/sz)*100;

    opts_->gene_info_queue_.addqueue(pgeneinfo);
    opts_->gene_queue_.addqueue(pgenedata);
}

void BinTask::othertask()
{
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];
    for(auto exp : exp_vec)
    {
        x = exp.x / m_bin;
        y = exp.y / m_bin;
        dnb = x<<32 | y;
        map_dnb[dnb]+=exp.count;
    }

    GeneInfo *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = new std::vector<Expression>;
    pgeneinfo->vecptr->reserve(map_dnb.size());

    GeneInfo2 *pgenedata = new GeneInfo2(m_geneid);
    // pgenedata->vecdataptr = new std::vector<int>;
    // pgenedata->vecdataptr->reserve(map_dnb.size()*3);
    auto itor_dnb = map_dnb.begin();
    Expression exp{0,0,0};
    for(;itor_dnb!=map_dnb.end();itor_dnb++)
    {
        exp.x = itor_dnb->first>>32;
        exp.y = itor_dnb->first & 0xFFFFFFFF;
        exp.count = itor_dnb->second;
        pgeneinfo->vecptr->emplace_back(exp);

        // pgenedata->vecdataptr->push_back(gene.x);
        // pgenedata->vecdataptr->push_back(gene.y);
        // pgenedata->vecdataptr->push_back(gene.cnt);
        if (exp.count > m_maxexp)
            m_maxexp = exp.count;
    }
    pgenedata->maxexp = m_maxexp;
    pgenedata->vecdataptr = pgeneinfo->vecptr;

    opts_->gene_info_queue_.addqueue(pgeneinfo);
    opts_->gene_queue_.addqueue(pgenedata);
}

void BinTask::doTask()
{
    if(m_bin == 1)
    {
        bin1task();
    }
    else if(m_bin == 100)
    {
        bin100task();
    }
    else
    {
        othertask();
    }
}