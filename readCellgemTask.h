/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-09 15:16:20
 * @Description: file content
 */
#ifndef GEFTOOLS_READCELLGEMTASK_H_
#define GEFTOOLS_READCELLGEMTASK_H_

#include <unordered_map>
#include "utils.h"
#include "thread_pool.h"
#include "gef.h"
#include "cgefUtil.h"
#include <sstream>
#include <sys/time.h>

class readCellgemTask:public ITask
{
public:
    readCellgemTask();
    ~readCellgemTask();
    void doTask();
protected:
    bool readbuf();
    int cuttail(char *pbuf);
    virtual int getInfo();
    virtual int mergeinfo();

protected:
    int m_buflen = 0;
    char *m_pbuf = nullptr;
    unordered_map<int, cgef_cell*> m_map_cell; 
    unordered_map<string, cgef_gene*> m_map_gene; //保存一个基因出现的细胞情况

    static string m_leftstr;
    static mutex m_readmtx; //读文件锁
    static mutex m_mergemtx; //合并锁

    int m_min_x = INT_MAX, m_min_y = INT_MAX, m_max_x = 0, m_max_y = 0;
};

//geneID  xPos    yPos    UMICount 
class readCellgemTask_raw:public readCellgemTask
{
public:
    readCellgemTask_raw(){};
    ~readCellgemTask_raw(){};
    int getInfo();
    int mergeinfo();
private:
    unordered_map<string, bgef_gene*> m_map_bgene;
};

//geneID  x       y       UMICount        label   tag(adjust/raw)
//Hnrnpdl 19414   15299   1       39054.0 adjust
class readCellgemTask_tag:public readCellgemTask
{
public:
    readCellgemTask_tag(){};
    ~readCellgemTask_tag(){};
    int getInfo();
    //int mergeinfo();
};

//geneID  xPos    yPos    UMICount        cellN   areaID
//NOC2L   30632   35146   2       15324   L-F1-l3
class readCellgemTask_areaID:public readCellgemTask
{
public:
    readCellgemTask_areaID(){};
    ~readCellgemTask_areaID(){};
    int getInfo();
    //int mergeinfo();
};

//geneID  x       y       MIDCount        cell
//Cr2     15653   20188   1       113231
class readCellgemTask_cell:public readCellgemTask
{
public:
    readCellgemTask_cell(){};
    ~readCellgemTask_cell(){};
    int getInfo();
    //int mergeinfo();
};



#endif