/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-16 11:28:09
 * @Description: file content
 */
#ifndef GEFTOOLS_GETCELLGENEDATA_H
#define GEFTOOLS_GETCELLGENEDATA_H

#include "gef.h"
#include <unordered_map>
#include "opencv2/opencv.hpp"
using namespace cv;

class getCellgeneData
{
public:
    getCellgeneData(string &bgef, string &cgef);
    ~getCellgeneData();
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
private:
    string m_bgef;
    string m_cgef;

    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    int m_min_x, m_min_y, m_max_x, m_max_y, m_resolution;
    unordered_map<uint64_t, vector<DnbExpression>> m_hash_vecdnb;
    unordered_map<uint32_t, Rect> m_hash_cellrect;
};


void getCellgeneData(string &bgef, string &cgef);


#endif