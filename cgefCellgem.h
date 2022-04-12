/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:25
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-07 17:34:55
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFCELLGEM_H_
#define GEFTOOLS_CGEFCELLGEM_H_

#include "cgef_writer.h"
#include "thread_pool.h"
#include "opencv2/opencv.hpp"
using namespace cv;

struct celldata
{
    celldata(int cid, int lid):c_idx(cid),l_idx(lid){}
    int c_idx; //轮廓idx
    int l_idx; //连通域idx
};


class cgefCellgem
{
public:
    cgefCellgem(/* args */);
    ~cgefCellgem();
    void readmask();
    void readxy();
    void readcellgem();
    //void mapCell();
    void writeFile(CgefWriter *cwptr);
    void writeAttr();
    //void writeBorder();
    void writeCell();
    void writeGene();
    
private:
    unsigned int m_maskcellnum; //从mask文件获取的细胞个数
    unsigned int m_blocknum;
    unsigned int m_block_size[4] = {0};

    int m_min_x{INT_MAX}, m_max_x{0}, m_min_y{INT_MAX}, m_max_y{0}, m_rows{0}, m_cols{0};
    Mat m_stats, m_centroids;
    vector<vector<Point>> m_contours;
    vector<vector<celldata>> m_vec_veccell;

    CgefWriter *m_cgefwPtr = nullptr;
    ThreadPool *m_thpoolPtr = nullptr;
    unordered_map<int, int> m_hash_clabel2cid; //建立从label到cellid的映射
    unordered_map<string, int> m_hash_gname2gid;//gname到geneid的映射
};

#endif