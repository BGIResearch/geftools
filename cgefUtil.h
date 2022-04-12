/*
 * @Author: zhaozijian
 * @Date: 2022-04-01 10:15:19
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-11 17:35:39
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFUTIL_H_
#define GEFTOOLS_CGEFUTIL_H_

#include <map>
#include <vector>
#include <string>

class cgef_cell
{
public:
    cgef_cell(){}
    ~cgef_cell(){}
    bool add(std::string &gene, unsigned short umi)
    {
        dnbcnt++;
        expcnt+=umi;
        // if(raw)
        // {
        //     m_xtotal += x;
        //     m_ytotal += y;
        // }

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
        // m_xtotal += cell.m_xtotal;
        // m_ytotal += cell.m_ytotal;

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
        return true;
    }

public:
    //long m_xtotal = 0, m_ytotal = 0;
    unsigned short expcnt = 0; //细胞umi之和
    unsigned short dnbcnt = 0; //细胞所有坐标点个数 >= gene个数
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

#endif