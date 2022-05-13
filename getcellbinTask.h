/*
 * @Author: zhaozijian
 * @Date: 2022-05-05 17:30:29
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-10 15:32:54
 * @Description: file content
 */
#ifndef GEFTOOL_GETCELLBINTASK_H
#define GEFTOOL_GETCELLBINTASK_H

#include "thread_pool.h"
#include "gef.h"
#include "cgefUtil.h"
#include <map>
#include <vector>
#include <string>
#include "cgefParam.h"
#include "opencv2/opencv.hpp"
using namespace cv;

class getcellbinTask:public ITask
{
public:
    getcellbinTask(const std::string &gene, bgef_gene *ptr, 
                GeneData* gdptr, vector<GeneExpData> &vecgexp, 
                vector<bgef_cell *> &veccellexp,
                unordered_map<uint32_t, uint32_t> &cl2cid,
                Mat &outimg):
                m_gene(gene), m_gptr(ptr),m_gdptr(gdptr),
                m_vec_geneexp(vecgexp), m_vec_cellexp(veccellexp),
                m_map_clabel2cid(cl2cid), m_outimg(outimg)
    {
    }
    ~getcellbinTask()
    {

    }
    void doTask()
    {
        uint32_t cid = 0, cell_label = 0;
        std::map<uint32_t, std::vector<Expression*>> map_cell;
        for(Expression &exp : m_gptr->m_vecExp)
        {
            cell_label = m_outimg.at<uint32_t>(exp.x-cgefParam::GetInstance()->m_min_x, exp.y-cgefParam::GetInstance()->m_min_y);
            if(cell_label)
            {
                cid = m_map_clabel2cid[cell_label];
                map_cell[cid].emplace_back(&exp);

                m_expsum += exp.count;
                m_maxExp = std::max(m_maxExp, exp.count);
            }
        }

        addGeneArry(map_cell);
    }

    void addGeneArry(std::map<uint32_t, std::vector<Expression*>> &map_cell)
    {
        lock_guard<mutex> lock(m_mutex_gene);
        memcpy(m_gdptr[m_gid].gene_name, m_gene.c_str(), m_gene.length());
        m_gdptr[m_gid].cell_count = map_cell.size();
        m_gdptr[m_gid].exp_count = m_expsum;
        m_gdptr[m_gid].max_mid_count = m_maxExp;
        m_gdptr[m_gid].offset = m_offset;
        m_offset += map_cell.size();
        m_gid++;

        auto itor = map_cell.begin();
        for(;itor != map_cell.end();itor++)
        {
            m_vec_geneexp.emplace_back(itor->first, map_cell.size());

            m_vec_cellexp[itor->first]->add(m_gid, itor->second.expcnt, itor->second.dnbcnt);
        }
        
        cgefParam::GetInstance()->m_maxExp_gexp = std::max(cgefParam::GetInstance()->m_maxExp_gexp, m_maxExp);
        cgefParam::GetInstance()->m_minExp = std::min(cgefParam::GetInstance()->m_minExp, m_expsum);
        cgefParam::GetInstance()->m_maxExp = std::max(cgefParam::GetInstance()->m_maxExp, m_expsum);
        cgefParam::GetInstance()->m_minCell = std::min(cgefParam::GetInstance()->m_minCell, m_gdptr[m_gid].cell_count);
        cgefParam::GetInstance()->m_maxCell = std::max(cgefParam::GetInstance()->m_maxCell, m_gdptr[m_gid].cell_count);
    }

private:
    static uint16_t m_gid ;
    static uint32_t m_offset;
    static mutex m_mutex_gene;
    
    uint32_t m_maxExp = 0; //max MID count of current gene
    uint32_t m_expsum = 0; 
    bgef_gene *m_gptr = nullptr;
    GeneData *m_gdptr = nullptr;
    Mat &m_outimg;
    const std::string &m_gene;
    vector<GeneExpData> &m_vec_geneexp;
    vector<bgef_cell *> &m_vec_cellexp;
    unordered_map<uint32_t, uint32_t> &m_map_clabel2cid;
};

#endif