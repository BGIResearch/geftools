/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:25
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-13 14:12:31
 * @Description: file content
 */
#ifndef GEFTOOLS_CGEFCELLGEM_H_
#define GEFTOOLS_CGEFCELLGEM_H_

#include "cgef_writer.h"
#include "thread_pool.h"
#include "opencv2/opencv.hpp"
#include "cgefUtil.h"
using namespace cv;

struct celldata
{
    celldata(int cid, int lid):c_idx(cid),l_idx(lid){}
    int c_idx; //轮廓idx
    int l_idx; //连通域idx
};

struct cellt
{
    uint16_t expcnt;
    uint16_t dnbcnt;
};


class cgefCellgem
{
public:
    cgefCellgem();
    ~cgefCellgem();
    void readmask(const string &strmask);
    void readxy(const string &strrawgem);
    void readcellgem(const string &strinput);
    void writeFile(CgefWriter *cwptr, const string &strmask, const string &strinput, const string &strrawgem);
    void writeAttr();
    void writeCell();
    void writeGene();
    void getCelldata();
    void getCelldata_celltype();
    void writeCell_celltype();
    void addCellborder(int cx, int cy, vector<short> &vec_border, uint32_t idx);

    void readcellgem_new();

    void clabeltocid();
    void writeGene_raw();
    void writeCell_raw();

    void readBgef(const string &strinput);
    void writeGene_bgef();
private:
    unsigned int m_maskcellnum = 0; //从mask文件获取的细胞个数
    unsigned int m_blocknum;
    unsigned int m_block_size[4] = {0};

    int m_min_x{INT_MAX}, m_max_x{0}, m_min_y{INT_MAX}, m_max_y{0}, m_rows{0}, m_cols{0};
    Mat m_stats, m_outimg, m_centroids;
    vector<vector<Point>> m_contours;
    vector<vector<celldata>> m_vec_veccell;

    CgefWriter *m_cgefwPtr = nullptr;
    ThreadPool *m_thpoolPtr = nullptr;
    unordered_map<uint32_t, uint32_t> m_hash_clabel2cid; //建立从label到cellid的映射
    unordered_map<string, int> m_hash_gname2gid;//gname到geneid的映射
    unordered_map<string, int> m_hash_celltype;//celltype到typeid的映射
    vector<string> m_vec_celltype;

    vector<uint32_t> m_vec_blkidx;
    vector<uint32_t> m_vec_cellLabel;
    //map<uint32_t, bgef_cell *> m_mapcellexp;
    vector<bgef_cell *> m_vec_cellexp;

    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
};

#endif