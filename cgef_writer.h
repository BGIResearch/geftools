/** @file cgef_writer.h
    @brief Declare a CgefWriter class for writing cell bin gef.

    Created by huangzhibo on 2021/12/27.
*/

#ifndef GEFTOOLS_CGEF_WRITER_H
#define GEFTOOLS_CGEF_WRITER_H

#include <vector>
#include <map>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "utils.h"
#include "gef.h"
#include "mask.h"
#include "common_bin.h"

class CgefWriter {
  public:
    explicit CgefWriter(const string& output_cell_gef);
    ~CgefWriter();


    /**
     * @brief Add dnb expression info of one cell
     *
     * This method can only be used when the class is constructed in mode="w"
     * @param dnb_coordinates A vector of dnb coordinates inner one cell region
     * @param bin_gene_exp_map  A map abort geneID and expCount of the genes for each bin, key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     * @param center_point Center point of cell polygon
     * @param area The polygon area of the cell
     */
    void addDnbExp(vector<Point> & dnb_coordinates,
    map<unsigned long long int, vector<CellExpData>> & bin_gene_exp_map,
    const Point& center_point,
    unsigned short area);

    static unsigned short calcMaxCountOfGeneExp(vector<GeneExpData> & gene_exps);

    void storeAttr(CellBinAttr& cell_attr) const;
    void storeCell();
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
    int write(CommonBin& common_bin_gef, Mask& mask);

  private:
    hid_t file_id_;
    hid_t group_id_;
    hid_t str32_type_;
    string bin_gef_;
    string mask_file_;
    map<unsigned short, vector<GeneExpData>> gene_exp_map_;
    vector<CellData> cell_list_;
    vector<CellExpData> cell_exp_list_;
    vector<S32> cell_type_list_;
    CellAttr cell_attr_ = {
        .min_gene_count=USHRT_MAX,
        .min_exp_count=USHRT_MAX,
        .min_dnb_count=USHRT_MAX,
        .min_area=USHRT_MAX,
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

};

#endif //GEFTOOLS_CGEF_WRITER_H
