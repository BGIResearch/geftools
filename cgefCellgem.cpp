/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:30
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-03-25 17:32:01
 * @Description: file content
 */

#include "cgefCellgem.h"
#include "cgefParam.h"
#include "readCellgemTask.h"
#include "mask.h"


cgefCellgem::cgefCellgem(/* args */)
{
    cgefParam::GetInstance()->m_infile = gzopen(opts->input_file_.c_str(), "r");
    gzbuffer(opts->infile_, READLEN);
}

cgefCellgem::~cgefCellgem()
{
    gzclose(cgefParam::GetInstance()->m_infile);
}

void cgefCellgem::readcellgem()
{
    std::string line;
    while (readline(cgefParam::GetInstance()->m_infile, line))
    {
        if (line.substr(0, 6) == "geneID") break;
    }

    ThreadPool thpool(cgefParam::GetInstance()->m_threadcnt);
    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        readCellgemTask *rtask = new readCellgemTask();
        thpool.addTask(rtask);
    }
    thpool.waitTaskDone();
}

struct cellinfo
{
    int celllabel;
    vector<GeneExp> *ptr;
};


void cgefCellgem::mapCell(Mask *maskptr)
{
    //把hash分成若干段，利用多线程去做映射
    int sz = cgefParam::GetInstance()->m_map_cell.size();
    int tsz = sz/cgefParam::GetInstance()->m_threadcnt + 1;

    vector<cellinfo> veccell;
    auto itor = cgefParam::GetInstance()->m_map_cell.begin();
    for(;itor != cgefParam::GetInstance()->m_map_cell.end();itor++)
    {
        cellinfo ci;
        ci.celllabel = itor->first;
        ci.ptr = &(itor->second);
        veccell.emplace(std::move(ci));
    }

    
}



