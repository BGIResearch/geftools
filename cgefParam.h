/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:56:17
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-09 10:31:03
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFPARAM_H_
#define GEFTOOLS_CGEFPARAM_H_

#include <zlib.h>
#include <unordered_map>
#include <vector>
#include "gef.h"
#include "cgefUtil.h"
#include "BufPool.h"

enum InputType
{
    INPUTTYPE_BGEF = 0,
    INPUTTYPE_GEM_MASK,
    INPUTTYPE_GEM_TAGMASK,
    INPUTTYPE_GEM_AREAID,
    INPUTTYPE_GEM_CELL
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
    std::string m_rawgemstr;
    std::string m_maskstr;
    std::string m_cellgemstr;
    gzFile m_infile;
    std::unordered_map<int, cgef_cell*> m_map_cell;
    std::unordered_map<std::string, cgef_gene*> m_map_gene;
    int m_min_x = INT_MAX, m_min_y = INT_MAX, m_max_x = 0, m_max_y = 0;
    InputType m_intype = INPUTTYPE_BGEF;
    BufPool *m_bpPtr = nullptr;
    char *m_pdata = nullptr;

    unordered_map<string, bgef_gene*> m_map_bgene;
    uint32_t m_minExp = UINT_MAX, m_maxExp = 0, m_minCell = UINT_MAX, m_maxCell = 0, m_maxExp_gexp = 0;
private:
    cgefParam(/* args */){};
    ~cgefParam(){};

};

#endif