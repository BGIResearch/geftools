/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:37
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 14:00:17
 * @Description: file content
 */

#include "readCellgemTask.h"
#include "cgefParam.h"
#include "BufPool.h"

string readCellgemTask::m_leftstr;
mutex readCellgemTask::m_readmtx;
mutex readCellgemTask::m_mergemtx;

readCellgemTask::readCellgemTask(int type):m_type(type)
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

    switch (m_type)
    {
    case 1: //geneID  x       y       UMICount        label   tag
    {
        while (brun)
        {
            brun = readbuf();
            getCellInfo();
        }

        mergeCellinfo();
    }
        break;
    case 2: //geneID  x       y       MIDCount
    {
        while (brun)
        {
            brun = readbuf();
            getInfo();
        }

        mergeinfo();
    }
        break;
    case 3: //geneID  xPos    yPos    UMICount        cellN   areaID
    {
        while (brun)
        {
            brun = readbuf();
            getInfo_celltype();
        }

        mergeCellinfo();
    }
        break;
    case 4:
    {
        while (brun)
        {
            brun = getInfo_celltype_new();
        }
        mergeCellinfo();
    }
    break;
    default:
        break;
    }
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

    //char gname[64]={0};
    string gname;
    int len = 0, x = 0, y=0, umi=0;
    for(;i<m_buflen;i++)
    {
        if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                len = &m_pbuf[i]-ptr;
                gname.clear();
                gname.append(ptr, len);
                // memcpy(gname, ptr, len);
                // gname[len]='\0';
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 1:
                //x = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 2:
                //y = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 3:
                umi = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 4:
                celllabel = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 5:
                k = 0;
                if(m_map_cell.find(celllabel) == m_map_cell.end())
                {
                    cgef_cell *cptr = new cgef_cell(celllabel);
                    m_map_cell.emplace(celllabel, cptr);
                }
                m_map_cell[celllabel]->add(gname,umi);

                if(m_map_gene.find(gname) == m_map_gene.end())
                {
                    cgef_gene *gptr = new cgef_gene();
                    m_map_gene.emplace(gname, gptr);
                }
                m_map_gene[gname]->add(celllabel, umi);
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

    // auto itor = m_map_cell.begin();
    // auto &tmap_cell = cgefParam::GetInstance()->m_map_cell;
    // for(;itor!=m_map_cell.end();itor++)
    // {
    //     if(tmap_cell.find(itor->first) != tmap_cell.end())
    //     {
    //         tmap_cell[itor->first]->merge(*(itor->second));
    //         delete itor->second;
    //     }
    //     else
    //     {
    //         tmap_cell.emplace(itor->first, itor->second);
    //     }
    // }

    // auto itor_g = m_map_gene.begin();
    // auto &tmp_gene = cgefParam::GetInstance()->m_map_gene;
    // for(;itor_g!=m_map_gene.end();itor_g++)
    // {
    //     if(tmp_gene.find(itor_g->first) != tmp_gene.end())
    //     {
    //         tmp_gene[itor_g->first]->merge(*(itor_g->second));
    //         delete itor_g->second;
    //     }
    //     else
    //     {
    //         tmp_gene.emplace(itor_g->first, itor_g->second);
    //     }
    // }

    //if(m_type == 3)
    {
        cgefParam::GetInstance()->m_min_x = std::min(cgefParam::GetInstance()->m_min_x, m_min_x);
        cgefParam::GetInstance()->m_max_x = std::max(cgefParam::GetInstance()->m_max_x, m_max_x);
        cgefParam::GetInstance()->m_min_y = std::min(cgefParam::GetInstance()->m_min_y, m_min_y);
        cgefParam::GetInstance()->m_max_y = std::max(cgefParam::GetInstance()->m_max_y, m_max_y);
    }

    return 0;
}

///////////////////////
int readCellgemTask::getInfo()
{
    int i = 0, k = 0, celllabel=0;
    char *ptr = m_pbuf;

    int len = 0, x = 0, y=0, umi=0;
    for(;i<m_buflen;i++)
    {
        if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 1:
                x = atoi(ptr);
                m_min_x = std::min(m_min_x, x);
                m_max_x = std::max(m_max_x, x);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 2:
                y = atoi(ptr);
                m_min_y = std::min(m_min_y, y);
                m_max_y = std::max(m_max_y, y);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 3:
                k++;
                cgefParam::GetInstance()->m_pdata[y*22747+x] = 1;
                ptr = &m_pbuf[i+1];
                break;
            default:
                break;
            }
        }
    }
}

int readCellgemTask::mergeinfo()
{
    lock_guard<mutex> lock(m_mergemtx);
    cgefParam::GetInstance()->m_min_x = std::min(cgefParam::GetInstance()->m_min_x, m_min_x);
    cgefParam::GetInstance()->m_min_y = std::min(cgefParam::GetInstance()->m_min_y, m_min_y);
    cgefParam::GetInstance()->m_max_x = std::max(cgefParam::GetInstance()->m_max_x, m_max_x);
    cgefParam::GetInstance()->m_max_y = std::max(cgefParam::GetInstance()->m_max_y, m_max_y);
}

int readCellgemTask::getInfo_celltype()
{
    int i = 0, k = 0, celllabel=0;
    char *ptr = m_pbuf;

    string gname;
    int len = 0, x = 0, y=0, umi=0;
    for(;i<m_buflen;i++)
    {
        if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                len = &m_pbuf[i]-ptr;
                gname.clear();
                gname.append(ptr, len);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 1:
                x = atoi(ptr);
                m_min_x = std::min(m_min_x, x);
                m_max_x = std::max(m_max_x, x);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 2:
                y = atoi(ptr);
                m_min_y = std::min(m_min_y, y);
                m_max_y = std::max(m_max_y, y);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 3:
                umi = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 4:
                celllabel = atoi(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 5:
                k = 0;
                if(celllabel > 0 && (memcmp(ptr,"NA",2)!=0) )
                {
                    if(m_map_cell.find(celllabel) == m_map_cell.end())
                    {
                        cgef_cell *cptr = new cgef_cell(celllabel, ptr, &m_pbuf[i]-ptr);
                        m_map_cell.emplace(celllabel, cptr);
                    }
                    m_map_cell[celllabel]->add(gname,umi, x, y);

                    if(m_map_gene.find(gname) == m_map_gene.end())
                    {
                        cgef_gene *gptr = new cgef_gene();
                        m_map_gene.emplace(gname, gptr);
                    }
                    m_map_gene[gname]->add(celllabel, umi);
                }
                ptr = &m_pbuf[i+1];
                break;
            default:
                break;
            }
        }
    }

    return m_map_cell.size();
}


int readCellgemTask::getInfo_celltype_new()
{
    Membuf *mptr = cgefParam::GetInstance()->m_bpPtr->getptr();
    if(mptr == nullptr) return 0;
    int i = 0, k = 0, celllabel=0;
    char *ptr = mptr->m_buf;

    string gname;
    int len = 0, x = 0, y=0, umi=0;
    for(;i<mptr->m_len;i++)
    {
        if(mptr->m_buf[i] == '\t' || mptr->m_buf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                len = mptr->m_buf+i-ptr;
                gname.clear();
                gname.append(ptr, len);
                k++;
                ptr = mptr->m_buf+i+1;
                break;
            case 1:
                x = atoi(ptr);
                m_min_x = std::min(m_min_x, x);
                m_max_x = std::max(m_max_x, x);
                k++;
                ptr = mptr->m_buf+i+1;
                break;
            case 2:
                y = atoi(ptr);
                m_min_y = std::min(m_min_y, y);
                m_max_y = std::max(m_max_y, y);
                k++;
                ptr = mptr->m_buf+i+1;
                break;
            case 3:
                umi = atoi(ptr);
                k++;
                ptr = mptr->m_buf+i+1;
                break;
            case 4:
                celllabel = atoi(ptr);
                k++;
                ptr = mptr->m_buf+i+1;;
                break;
            case 5:
                k = 0;
                if(celllabel > 0 && (memcmp(ptr,"NA",2)!=0) )
                {
                    if(m_map_cell.find(celllabel) == m_map_cell.end())
                    {
                        cgef_cell *cptr = new cgef_cell(celllabel, ptr, mptr->m_buf+i-ptr);
                        m_map_cell.emplace(celllabel, cptr);
                    }
                    m_map_cell[celllabel]->add(gname,umi, x, y);

                    if(m_map_gene.find(gname) == m_map_gene.end())
                    {
                        cgef_gene *gptr = new cgef_gene();
                        m_map_gene.emplace(gname, gptr);
                    }
                    m_map_gene[gname]->add(celllabel, umi);
                }
                ptr = mptr->m_buf+i+1;;
                break;
            default:
                break;
            }
        }
    }

    delete mptr;
    return m_map_cell.size();
}