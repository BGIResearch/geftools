/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-17 13:46:14
 * @Description: file content
 */
#ifndef GEFTOOLS_CELLADJUST_H
#define GEFTOOLS_CELLADJUST_H

#include "gef.h"
#include <unordered_map>
#include "opencv2/opencv.hpp"
using namespace cv;

class cellAdjust
{
public:
    cellAdjust(string &bgef, string &cgef);
    ~cellAdjust();
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    void 
private:
    string m_bgef;
    string m_cgef;

    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    int m_min_x, m_min_y, m_max_x, m_max_y, m_resolution;
    unordered_map<uint64_t, vector<DnbExpression>> m_hash_vecdnb;
    unordered_map<uint32_t, Rect> m_hash_cellrect;
};


void cellAdjust(string &bgef, string &cgef);


#endif