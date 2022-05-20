/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-19 16:22:41
 * @Description: file content
 */
#ifndef GEFTOOLS_CELLADJUST_H
#define GEFTOOLS_CELLADJUST_H

#include "gef.h"
#include "cgef_writer.h"
#include <unordered_map>
#include "opencv2/opencv.hpp"
using namespace cv;

struct Dnbs
{
    Dnbs(uint16_t gid, uint16_t cnt):geneid(gid),midcnt(cnt){};
    uint16_t geneid;
    uint16_t midcnt;
};

struct cellgem_label
{
    cellgem_label(uint32_t geneid, int x, int y, uint32_t midcnt, uint32_t cellid):
    geneid(geneid),x(x),y(y),midcnt(midcnt),cellid(cellid)
    {
    }
    uint32_t geneid;
    int x;
    int y;
    uint32_t midcnt;
    uint32_t cellid;
};



class cellAdjust
{
public:
    cellAdjust(string &outpath);
    ~cellAdjust();
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    uint32_t getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem);
    void writeCellAdjust(vector<Cell> &veccell, vector<DnbExpression> &vecDnb);
    bool addborder(unsigned int cid, vector<Point> &vecPoint, vector<Point> &border, vector<short> &vec_border);
    void writeCell(vector<Cell> &veccell, vector<DnbExpression> &vecDnb);
    void writeGene();
private:
    uint32_t m_genencnt;
    uint32_t m_geneexpcnt;
    uint32_t m_cellcnt;
    int m_offsetX, m_offsetY;
    vector<string> m_vecgenename;
    int m_min_x, m_min_y, m_max_x, m_max_y, m_resolution;
    unordered_map<uint64_t, vector<Dnbs>> m_hash_vecdnb;
    unordered_map<uint32_t, Rect> m_hash_cellrect;
    Mat m_fill_points;
    unsigned int *m_blkidxPtr = nullptr;
    unsigned int m_block_size[4];
    CellData *m_cell_arrayptr = nullptr;
    CgefWriter *m_cgefwPtr = nullptr;
    map<uint32_t, vector<GeneExpData>> m_map_gene;
};



#endif