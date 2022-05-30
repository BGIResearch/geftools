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
#include <set>
#include <vector>
#include <string>
#include "cgefParam.h"
#include "gene_queue.h"
#include "opencv2/opencv.hpp"
using namespace cv;

class cellbin
{
public:
    cellbin(int x, int y, uint16_t area, uint32_t label, unsigned int *block_size):
    m_cx(x), m_cy(y), m_area(area), m_label(label)
    {
        m_blkid = cx/block_size[0] + (cy/block_size[1])*block_size[2];
    }
    ~cellbin(){};
    void add(vector<Dnbs_exon> & vecdnb)
    {
        for(Dnbs_exon &dnb : vecdnb)
        {
            if(m_setgid.find(dnb.geneid)==m_setgid.end())
            {
                m_setgid.insert(dnb.geneid);
                m_genecnt++;
            }
            m_dnbcnt++;
            m_expcnt+=dnb.midcnt;
            m_exoncnt+=dnb.exoncnt;
            m_vecCexp.push_back(dnb.geneid, dnb.midcnt);
            m_vecCExon.push_back(dnb.exoncnt);
        }
    }
public:
    uint16_t m_genecnt = 0;
    uint16_t m_expcnt = 0;
    uint16_t m_dnbcnt = 0;
    uint16_t m_area = 0;
    uint16_t m_exoncnt = 0;
    int m_cx = 0;
    int m_cy = 0;
    uint32_t m_label = 0;
    uint32_t m_blkid = 0;
    vector<CellExpData> m_vecCexp;
    vector<uint16_t> m_vecCExon;
    vector<Point> m_vecborder;
    set<uint16_t> m_setgid;
};

class cgefCellgem;
class getcellbinTask:public ITask
{
public:
    getcellbinTask(cgefCellgem *ptr, uint32_t label, Rect &rect, vector<Point> &vecpoint):
    m_cgefPtr(ptr),m_label(label),m_rect(rect),m_vecpoint(vecpoint)
    {
    }
    ~getcellbinTask()
    {

    }
    void doTask()
    {
        int cx = m_cgefPtr->m_centroids.at<double>(m_label,0);
        int cy = m_cgefPtr->m_centroids.at<double>(m_label,1);
        int area = m_cgefPtr->m_stats.at<int>(m_label, CC_STAT_AREA);

        cellbin *cptr = new cellbin(cx, cy, area, m_label, m_cgefPtr->m_block_size);
        uint64_t l_id = 0;
        for(int i=m_rect.x;i<m_rect.width;i++)
        {
            for(int j=m_rect.y;j<m_rect.height;j++)
            {
                if(m_cgefPtr->m_outimg.at<int>(j,i) == m_label)
                {
                    l_id = i;
                    l_id = (l_id<<32) | j;
                    auto itor = m_cgefPtr->m_hash_vecdnb.find(l_id);
                    if(itor != m_cgefPtr->m_hash_vecdnb.end())
                    {
                        cptr->add(itor->second);
                    }
                }
            }
        }

        if(cptr->m_genecnt)
        {
            getborder(cptr);
        }

        m_cgefPtr->m_cellqueuePtr->addqueue(cptr);
    }

    void getborder(cellbin *cptr)
    {
        int i=0;
        int sz = m_vecpoint.size();
        if(sz > BORDERCNT)
        {
            vector<Point> tmpborder;
            double epsilon = 0.01 * arcLength(m_vecpoint, true);
            approxPolyDP(m_vecpoint, tmpborder, epsilon, true);

            sz = tmpborder.size();
            for(;i<sz;i++)
            {
                cptr->m_vecborder.emplace_back(tmpborder[i].x - cx);
                cptr->m_vecborder.emplace_back(tmpborder[i].y - cy);
            }
        }
        else
        {
            for(;i<sz;i++)
            {
                cptr->m_vecborder.emplace_back(vborder[i].x - cx);
                cptr->m_vecborder.emplace_back(vborder[i].y - cy);
            }
        }
        
        for(;i<BORDERCNT;i++) //不足补0
        {
            cptr->m_vecborder.emplace_back(SHRT_MAX);
            cptr->m_vecborder.emplace_back(SHRT_MAX);
        }
    }

private:
    uint32_t m_label = 0;
    Rect &m_rect;
    vector<Point> &m_vecpoint;
    cgefCellgem *m_cgefPtr = nullptr;
};

#endif