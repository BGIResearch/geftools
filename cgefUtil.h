/*
 * @Author: zhaozijian
 * @Date: 2022-04-01 10:15:19
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-22 16:24:28
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
    cgef_cell(int label):m_celllabel(label){}
    cgef_cell(int label, char *ptr,int len):
    m_celllabel(label)
    {
        m_bsort = true;
        memcpy(m_type, ptr, len);
        m_type[len]='\0';
        m_vecx.reserve(1000);
        m_vecy.reserve(1000);
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
        m_vecx.push_back(x);
        m_vecy.push_back(y);

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
            m_vecx.insert(m_vecx.end(), cell.m_vecx.begin(), cell.m_vecx.end());
            m_vecy.insert(m_vecy.end(), cell.m_vecy.begin(), cell.m_vecy.end());
        }
        return true;
    }

    void getCenter(unsigned int *block_size, int offx, int offy)
    {
        int sz = m_vecx.size();
        if(sz < 2)
        {
            m_centerx = m_vecx[0] - offx;
            m_centery = m_vecy[0] - offy;
        }
        else
        {
            sort(m_vecx.begin(), m_vecx.end(), [](int a, int b){return a<b;});
            sort(m_vecy.begin(), m_vecy.end(), [](int a, int b){return a<b;});

            int pos = ceil( (sz + 1)*1.0 /2 );
            m_centerx = ceil (m_vecx[pos-2]*0.5+m_vecx[pos-1]*0.5) - offx; //中位数
            m_centery = ceil (m_vecy[pos-2]*0.5+m_vecy[pos-1]*0.5) - offy;

            vector<int> tmp1,tmp2;
            m_vecx.swap(tmp1);
            m_vecy.swap(tmp2);
        }

        m_blkid = m_centerx/block_size[0] + (m_centery/block_size[1])*block_size[2];
    }

public:
    bool m_bsort = false;
    long m_centerx = 0, m_centery = 0;
    vector<int> m_vecx;
    vector<int> m_vecy;
    int m_blkid = 0;
    int m_celllabel = 0;
    unsigned short expcnt = 0; //细胞umi之和
    unsigned short dnbcnt = 0; //细胞所有坐标点个数 >= gene个数
    char m_type[32]={0};
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