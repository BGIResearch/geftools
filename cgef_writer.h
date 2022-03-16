/** @file cgef_writer.h
    @brief Declare a CgefWriter class for writing cell bin gef.

    Created by huangzhibo on 2021/12/27.
*/

#ifndef GEFTOOLS_CGEF_WRITER_H
#define GEFTOOLS_CGEF_WRITER_H

#include <vector>
#include <map>
#include <unordered_set>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "utils.h"
#include "gef.h"
#include "mask.h"
#include "bgef_reader.h"

struct block
{
    block(int off, int cnt):offset(off), count(cnt){};
    int offset;
    int count;
};

class CgefWriter {
  public:
    explicit CgefWriter(const string& output_cell_gef, bool verbose = false);
    ~CgefWriter();


    /**
     * @brief Add dnb expression info of one cell
     *
     * This method can only be used when the class is constructed in mode="w"
     * @param dnb_coordinates A vector of dnb coordinates inner one cell region
     * @param bin_gene_exp_map  A map abort geneID and expCount of the genes for each bin, key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     * @param center_point Center point of the cell polygon
     * @param area The polygon area of the cell
     */

    void addDnbExp(vector<Point> & dnb_coordinates,
                   map<unsigned long long int, pair<unsigned int, unsigned short>> & bin_gene_exp_map,
                   const DnbExpression *dnb_expression,
                   const Point& center_point,
                   unsigned short area);

    static unsigned short calcMaxCountOfGeneExp(vector<GeneExpData> & gene_exps);

    void storeAttr(CellBinAttr& cell_attr) const;
    void storeCell(unsigned int block_num, unsigned int *block_index, const unsigned int * block_size);
    void storeCellExp();
    void storeCellBorder(char* borderPath, unsigned int cell_num) const;
    void storeCellBorderWithAttr(char* borderPath, unsigned int cell_num, unsigned int* effective_rect) const;
    void storeCellTypeList();

    /**
     * @brief Writing the contents of geneData and geneExpData to GEF.
     * @param gene_name_list
     */
    void storeGeneAndGeneExp(const vector<string> &gene_name_list);

    /**
     * @brief Writing to cgef.
     * @return
     */
    int write(BgefReader& common_bin_gef, Mask& mask);

    bool isVerbose() const;

    void setVerbose(bool verbose);

    unsigned short getRandomCellTypeNum() const;

    void setRandomCellTypeNum(unsigned short random_cell_type_num);

    int addLevel();
    void getblkcelldata_top(int lev, int cnt, int left);
    void getblkcelldata_bottom(int lev, int cnt);
    void getblkcelldata(int lev, int cnt, int left);
    void createBlktype();
    void writeCelldata(int lev, vector<block> &blk, vector<int> &vecid);
    

  private:
    hid_t file_id_;
    hid_t group_id_;
    hid_t str32_type_;
    // string bin_gef_;
    // string mask_file_;
    map<unsigned short, vector<GeneExpData>> gene_exp_map_;
    vector<CellData> cell_list_;
    vector<CellExpData> cell_exp_list_;
    vector<S32> cell_type_list_;
    char *m_borderptr = nullptr;
    int m_x_len = 0;
    int m_y_len = 0;
    unordered_set<int> m_hash_cellid;
    hid_t m_level_gid;
    hid_t m_blk_memtype;
    hid_t m_blk_filetype;

    CellAttr cell_attr_ = {
        .average_gene_count=0.0,
        .average_exp_count=0.0,
        .average_dnb_count=0.0,
        .average_area=0.0,
        .median_gene_count=0.0,
        .median_exp_count=0.0,
        .median_dnb_count=0.0,
        .median_area=0.0,
        .min_x=USHRT_MAX,
        .min_y=USHRT_MAX,
        .min_gene_count=USHRT_MAX,
        .min_exp_count=USHRT_MAX,
        .min_dnb_count=USHRT_MAX,
        .min_area=USHRT_MAX,
        .max_x=0,
        .max_y=0,
        .max_gene_count=0,
        .max_exp_count=0,
        .max_dnb_count=0,
        .max_area=0};
    unsigned long long int exp_count_sum_ = 0;
    unsigned long long int dnb_count_sum_ = 0;
    unsigned long long int area_sum_ = 0;
    unsigned short gene_num_ = 0;
    unsigned int cell_num_ = 0;
    unsigned int expression_num_ = 0;
    unsigned short max_mid_count_ = 0;
    unsigned short random_cell_type_num_ = 0;
    bool verbose_ = false;

};

#endif //GEFTOOLS_CGEF_WRITER_H
