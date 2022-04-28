/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:56:17
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-28 10:04:25
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFPARAM_H_
#define GEFTOOLS_CGEFPARAM_H_

#include <zlib.h>
#include <unordered_map>
#include <vector>
#include "gef.h"
#include "cgefUtil.h"

enum InputType
{
    INPUTTYPE_BGEF = 0,
    INPUTTYPE_GEM,
    INPUTTYPE_GEM_ADJUST,
    INPUTTYPE_GEM_LABEL
};

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
    int m_sn = 500;
    std::string m_rawgemstr;
    std::string m_maskstr;
    std::string m_cellgemstr;
    gzFile m_infile;
    std::unordered_map<int, cgef_cell*> m_map_cell;
    std::unordered_map<std::string, cgef_gene*> m_map_gene;
    int m_min_x = INT_MAX, m_min_y = INT_MAX, m_max_x = 0, m_max_y = 0;
    InputType m_intype = INPUTTYPE_BGEF;
private:
    cgefParam(/* args */){};
    ~cgefParam(){};

};

#endif