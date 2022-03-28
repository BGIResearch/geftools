/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:37
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-03-25 16:52:49
 * @Description: file content
 */

#include "readCellgemTask.h"
#include "cgefParam.h"

string readCellgemTask::m_leftstr;
mutex readCellgemTask::m_readmtx;
mutex readCellgemTask::m_mergemtx;

readCellgemTask::readCellgemTask()
{
    m_pbuf = new char[READLEN];
}

readCellgemTask::~readCellgemTask()
{
    delete[] m_pbuf;
}

void readCellgemTask::doTask()
{
    bool brun = true;
    while (brun)
    {
        brun = readbuf();
        getGeneInfo();
    }

    mergeGeneinfo();
}

bool readCellgemTask::readbuf()
{
    lock_guard<mutex> lock(m_readmtx);
    char *pbuf = m_pbuf;
    int leftlen = m_leftstr.length();
    memcpy(pbuf, m_leftstr.c_str(), leftlen);
    m_leftstr.clear();
    pbuf += leftlen;
    int readlen = READLEN-leftlen;
    int reallen = gzread(cgefParam::GetInstance()->m_infile, pbuf, readlen);

    m_buflen = reallen;
    if(reallen == readlen)
    {
        cuttail(m_pbuf);
    }
    else
    {
        if (m_buflen != 0)
            m_buflen += leftlen;
        return false;
    }
    return true;
}

int readCellgemTask::cuttail(char *pbuf)
{
    int i = READLEN-1;
    for(;i>0;i--)
    {
        if(pbuf[i] == '\n')
        {
            break;
        }
    }

    m_buflen = i+1;
    m_leftstr.append(&pbuf[m_buflen], READLEN-m_buflen);
    return 0;
}

int readCellgemTask::getCellInfo()
{
    int i = 0, k = 0, celllabel=0;
    char *ptr = m_pbuf;

    GeneExpData gcell;
    GeneExp expression{};
    int len = 0;
    for(;i<m_buflen;i++)
    {
        if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                len = &m_pbuf[i]-ptr;
                memcpy(expression.gname, ptr, len);
                expression.gname[len]='\0';
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 1:
                expression.x = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 2:
                expression.y = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 3:
                expression.umi = atoi(ptr);
                gcell.count = expression.umi;
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 4:
                celllabel = atoi(ptr);
                gcell.cell_id = celllabel;
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 5:
                if(memcmp(ptr, "raw", 3) == 0)
                {
                    expression.braw = true;
                }
                else
                {
                    expression.braw = false;
                }
                k = 0;
                
                m_map_cell[celllabel].push_back(expression);
                m_map_gene[expression.gname].push_back(gcell);
                ptr = &m_pbuf[i+1];
                break;
            default:
                break;
            }
        }
    }

    return m_map_cell.size();
}

int readCellgemTask::mergeCellinfo()
{
    lock_guard<mutex> lock(m_mergemtx);

    auto itor = m_map_cell.begin();
    auto &tmap_cell = cgefParam::GetInstance()->m_map_cell;
    for(;itor!=m_map_cell.end();itor++)
    {
        std::vector<GeneExp> &vec = tmap_cell[itor->first];
        vec.insert(vec.end(), itor->second.begin(), itor->second.end());
    }
    return 0;
}