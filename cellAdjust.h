/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:15
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 10:55:23
 * @Description: file content
 */
#ifndef GEFTOOLS_CELLADJUST_H
#define GEFTOOLS_CELLADJUST_H

#include "gef.h"
#include "cgef_writer.h"
#include <unordered_map>
#include "opencv2/opencv.hpp"
using namespace cv;


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
    cellAdjust();
    ~cellAdjust();
    void readBgef(const string &strinput);
    void readCgef(const string &strinput);
    uint32_t getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem);
    void writeCellAdjust(const string &outpath, Cell *cellptr, int cellcnt, DnbExpression *dnbptr, int dnbcnt);
    bool addborder(unsigned int cid, vector<Point> &vecPoint, vector<Point> &border, vector<short> &vec_border);
    void writeCell(Cell *cellptr, unsigned int cellcnt, DnbExpression *dnbptr, unsigned int dnbcnt);
    void writeGene();
    void cgeftogem(const string &strbgef, const string &strcgef, const string &strout);
    void cgeftogem_exon(const string &strbgef, const string &strcgef, const string &strout);
    bool bexon(){return m_bexon;}
    void createRegionGef(const string &strout);
    void getRegionGenedata(vector<vector<int>> &m_vecpos);
private:
    bool m_bexon = false;
    uint32_t m_genencnt;
    uint32_t m_geneexpcnt;
    uint32_t m_cellcnt;
    int m_offsetX, m_offsetY;
    vector<string> m_vecgenename;
    int m_min_x, m_min_y, m_max_x, m_max_y;
    uint32_t m_resolution;
    unordered_map<uint64_t, vector<Dnbs_exon>> m_hash_vecdnb_exon;
    unordered_map<uint32_t, Rect> m_hash_cellrect;
    Mat m_fill_points;
    unsigned int m_block_size[4];
    CellData *m_cell_arrayptr = nullptr;
    CgefWriter *m_cgefwPtr = nullptr;
    map<uint32_t, vector<GeneExpData>> m_map_gene;
    BgefOptions *m_bgefopts = nullptr;

    char m_szomics[32]={0};
    int m_maxx = 0, m_maxy = 0;
    hid_t m_bgeffile_id = 0;
};



#endif