/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:56:17
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-03-25 17:25:52
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFPARAM_H_
#define GEFTOOLS_CGEFPARAM_H_

#include <zlib.h>
#include <unordered_map>

class cgefParam
{
public:
    static cgefParam *GetInstance()
    {
        static cgefParam instance;
        return &instance;
    }

    int m_threadcnt = 1;
    gzFile m_infile;
    std::unordered_map<int, std::vector<GeneExp>> m_map_cell;
private:
    cgefParam(/* args */){};
    ~cgefParam(){};

};

#endif