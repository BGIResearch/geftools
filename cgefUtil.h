/*
 * @Author: zhaozijian
 * @Date: 2022-04-01 10:15:19
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-07 17:34:53
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFUTIL_H_
#define GEFTOOLS_CGEFUTIL_H_

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "utils.h"

class cgef_cell
{
public:
    cgef_cell(int label):m_celllabel(label){}
    cgef_cell(int label, char *ptr,int len):
    m_celllabel(label)
    {
        m_bsort = true;
        memcpy(m_type, ptr, len);
        m_type[len]='\0';
    }
    ~cgef_cell(){}
    bool add(std::string &gene, unsigned short umi)
    {
        dnbcnt++;
        expcnt+=umi;

        if(m_map_cellexp.find(gene) != m_map_cellexp.end())
        {
            m_map_cellexp[gene] += umi;
        }
        else
        {
            m_map_cellexp.emplace(gene, umi);
        }
        return true;
    }

    bool add(std::string &gene, unsigned short umi, int x, int y)
    {
        dnbcnt++;
        expcnt+=umi;
        m_vecPoint.emplace_back(x,y);

        if(m_map_cellexp.find(gene) != m_map_cellexp.end())
        {
            m_map_cellexp[gene] += umi;
        }
        else
        {
            m_map_cellexp.emplace(gene, umi);
        }
        return true;
    }

    bool merge(cgef_cell &cell)
    {
        dnbcnt += cell.dnbcnt;
        expcnt += cell.expcnt;

        auto itor = cell.m_map_cellexp.begin();
        for(;itor!= cell.m_map_cellexp.end();itor++)
        {
            if(m_map_cellexp.find(itor->first) != m_map_cellexp.end())
            {
                m_map_cellexp[itor->first] += itor->second;
            }
            else
            {
                m_map_cellexp.emplace(itor->first, itor->second);
            }
        }

        if(m_bsort)
        {
            // m_vecx.insert(m_vecx.end(), cell.m_vecx.begin(), cell.m_vecx.end());
            // m_vecy.insert(m_vecy.end(), cell.m_vecy.begin(), cell.m_vecy.end());
            m_vecPoint.insert(m_vecPoint.end(), cell.m_vecPoint.begin(), cell.m_vecPoint.end());
        }
        return true;
    }

    bool getCenter(unsigned int *block_size, int offx, int offy)
    {
        if(m_vecPoint.size() <3) return false;
        vector<Point> tmp, hull;
        convexHull(m_vecPoint, hull, true);
        m_vecPoint.swap(tmp);
        int sz = hull.size();
        // if(sz < 2)
        // {
        //     m_centerx = m_vecx[0] - offx;
        //     m_centery = m_vecy[0] - offy;
        // }
        // else
        // {
        //     std::sort(m_vecx.begin(), m_vecx.end(), [](int a, int b){return a<b;});
        //     std::sort(m_vecy.begin(), m_vecy.end(), [](int a, int b){return a<b;});

        //     int pos = ceil( (sz + 1)*1.0 /2 );
        //     m_centerx = ceil (m_vecx[pos-2]*0.5+m_vecx[pos-1]*0.5) - offx; //中位数
        //     m_centery = ceil (m_vecy[pos-2]*0.5+m_vecy[pos-1]*0.5) - offy;

        //     vector<int> tmp1,tmp2;
        //     m_vecx.swap(tmp1);
        //     m_vecy.swap(tmp2);
        // }
        if(sz <= 2)
        {
            return false;
        }
        else
        {
            if(sz > BORDERCNT)
            {
                double epsilon = 0.01 * arcLength(hull, true);
                approxPolyDP(hull, m_border, epsilon, true); 
            }
            else
            {
                m_border.swap(hull);
            }

            Moments mu = moments(m_border, true);
            if(mu.m00 == 0) return false;    
            m_center = Point(static_cast<int>(mu.m10/mu.m00), static_cast<int>(mu.m01/mu.m00));
            m_area = mu.m00;
            
        }

        m_blkid = m_center.x/block_size[0] + (m_center.y/block_size[1])*block_size[2];
        assert(m_blkid < block_size[2] * block_size[3]);
        return true;
    }

public:
    bool m_bsort = true;
    Point m_center;
    vector<Point> m_vecPoint;
    vector<Point> m_border;
    int m_blkid = 0;
    int m_celllabel = 0;
    unsigned short expcnt = 0; //细胞umi之和
    unsigned short dnbcnt = 0; //细胞所有坐标点个数 >= gene个数
    char m_type[32]={0};
    double m_area;
    std::map<std::string, unsigned short> m_map_cellexp;//包含的基因情况
};

class cgef_gene
{
public:
    cgef_gene(){}
    ~cgef_gene(){}
    bool add(int label, unsigned short umi)
    {
        expcnt += umi;
        if(m_map_geneexp.find(label) != m_map_geneexp.end())
        {
            m_map_geneexp[label] += umi;
        }
        else
        {
            m_map_geneexp.emplace(label, umi);
        }
        return true;
    }

    bool merge(cgef_gene &gene)
    {
        expcnt += gene.expcnt;
        auto itor = gene.m_map_geneexp.begin();
        for(;itor != gene.m_map_geneexp.end();itor++)
        {
            if(m_map_geneexp.find(itor->first) != m_map_geneexp.end())
            {
                m_map_geneexp[itor->first] += itor->second;
            }
            else
            {
                m_map_geneexp.emplace(itor->first, itor->second);
            }
        }
        return true;
    }
public:
    unsigned int expcnt = 0; //基因umi之和
    std::map<int, unsigned short> m_map_geneexp;//出现的细胞情况
};

class bgef_gene
{
public:
    bgef_gene(){};
    ~bgef_gene(){};
    void add(uint32_t x, uint32_t y, uint32_t midcnt)
    {
        m_vecExp.emplace_back(x,y,midcnt);
    }
    void merge(bgef_gene &other)
    {
        m_vecExp.insert(m_vecExp.end(), other.m_vecExp.begin(), other.m_vecExp.end());
    }
public:
    vector<Expression> m_vecExp;
};

class bgef_cell
{
public:
    bgef_cell(){};
    ~bgef_cell(){};
    void add(uint16_t gid, uint16_t exp, uint16_t dnb)
    {
        m_vecCexp.emplace_back(gid, exp);
        m_expcnt += exp;
        m_dnbcnt += dnb;
    }
public:
    uint16_t m_expcnt;
    uint16_t m_dnbcnt;
    vector<CellExpData> m_vecCexp;
};

#endif