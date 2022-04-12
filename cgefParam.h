/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:56:17
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-07 11:24:56
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFPARAM_H_
#define GEFTOOLS_CGEFPARAM_H_

#include <zlib.h>
#include <unordered_map>
#include <vector>
#include "gef.h"
#include "cgefUtil.h"

class cgefParam
{
public:
    static cgefParam *GetInstance()
    {
        static cgefParam instance;
        return &instance;
    }

    int m_threadcnt = 1;
    int m_block_size[2]={256,256};
    std::string m_rawgemstr;
    std::string m_maskstr;
    std::string m_cellgemstr;
    gzFile m_infile;
    std::unordered_map<int, cgef_cell*> m_map_cell;
    std::unordered_map<std::string, cgef_gene*> m_map_gene;
    int m_min_x = INT_MAX, m_min_y = INT_MAX;
private:
    cgefParam(/* args */){};
    ~cgefParam(){};

};

#endif