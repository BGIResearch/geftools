/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-02 17:24:14
 * @Description: file content
 */
#ifndef GEFTOOLS_READCELLGEMTASK_H_
#define GEFTOOLS_READCELLGEMTASK_H_

#include <unordered_map>
#include "utils.h"
#include "thread_pool.h"
#include "gef.h"
#include "cgefUtil.h"

class readCellgemTask:public ITask
{
public:
    readCellgemTask(int type);
    ~readCellgemTask();
    void doTask();
private:
    bool readbuf();
    int cuttail(char *pbuf);
    int getCellInfo();
    int mergeCellinfo();
    int getInfo();
    int mergeinfo();
private:
    int m_type;
    int m_buflen = 0;
    char *m_pbuf = nullptr;
    unordered_map<int, cgef_cell*> m_map_cell; 
    unordered_map<string, cgef_gene*> m_map_gene; //保存一个基因出现的细胞情况

    static string m_leftstr;
    static mutex m_readmtx; //读文件锁
    static mutex m_mergemtx; //合并锁

    int m_min_x = INT_MAX, m_min_y = INT_MAX;
};

#endif