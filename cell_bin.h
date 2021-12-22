/** @file cell_bin.h
    @brief Declare a CellBin class and several related structs.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__CELL_BIN_H_
#define GEFTOOLS__CELL_BIN_H_

#include <vector>
#include <map>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "utils.h"

using namespace std;
using namespace cv;

/**
 * @brief Describe the Cell dataset in the cell bin GEF file
 */
struct CellData {
  unsigned int x; ///< Coordinate X of center point in this cell
  unsigned int y; ///< Coordinate Y of center point in this cell
  unsigned int offset;  ///< Offset of current cell in cellExp, 0-based
  unsigned short gene_count; ///< The number of gene in this cell
  unsigned short exp_count; ///< The total expression count of all genes in this cell
  unsigned short dnb_count; ///< Dnb number in this cell
  unsigned short area; ///< The polygon area of this cell
  unsigned short cell_type_id; ///< Cell type ID to index the CellTypeList
};

struct CellAttr {
    float average_gene_count;
    float average_exp_count;
    float average_dnb_count;
    float average_area;
    unsigned short min_gene_count;
    unsigned short min_exp_count;
    unsigned short min_dnb_count;
    unsigned short min_area;
    unsigned short max_gene_count;
    unsigned short max_exp_count;
    unsigned short max_dnb_count;
    unsigned short max_area;
};

struct GeneData {
    GeneData(const char* g, unsigned int o, unsigned int c, unsigned m)
    {
        int i = 0;
        while (g[i] != '\0')
        {
            gene_name[i] = g[i];
            ++i;
        }
        offset = o;
        cell_count = c;
        max_mid_count = m;
    }
    char gene_name[32] = {0};
    unsigned int offset;  ///< Offset of current gene in geneExp, 0-based
    unsigned int cell_count;
    unsigned short max_mid_count;  ///< max MID count of current gene
};

struct CellExpData {
//    explicit CellExpData(unsigned int gene_exp){
//        geneID = gene_exp >> 16;
//        count = gene_exp  & 0xFFFF;
//    }
//    CellExpData(unsigned short g, unsigned short c){
//        geneID = g;
//        count = c;
//    }
    unsigned short geneID;
    unsigned short count;
};

struct GeneExpData {
    unsigned int cell_id;
    unsigned short count;
};

struct CellBinAttr
{
    unsigned int version;
    unsigned int resolution;
    unsigned int offsetX;
    unsigned int offsetY;
};


class CellBin{
  private:
    hid_t file_id_;
    hid_t group_id_;
    hid_t str32_type;

    map<unsigned short, vector<GeneExpData>> gene_exp_map_;
    vector<CellData> cell_list_;
    vector<CellExpData> cell_exp_list_;
    vector<S32> cell_type_list_;
    vector<GeneData> gene_list_;
    vector<string> gene_name_list_;

private:

    vector<GeneExpData> gene_exp_list_;
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
    unsigned int gene_num_ = 0;
    unsigned int cell_num_ = 0;
    unsigned int expression_num_ = 0;

  public:
    CellBin(const string& filepath,  const string& mode);
    ~CellBin();

    unsigned int getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;

    /**
     * @brief Add dnb expression info of one cell
     *
     * This method can only be used when the class is constructed in mode="w"
     * @param dnb_coordinates A vector of dnb coordinates inner one cell region
     * @param bin_gene_exp_map  A map abort geneID and expCount of the genes for each bin, key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     * @param center_point Center point of cell polygon
     * @param area The polygon area of the cell
     */
    void addDnbInCell(vector<Point> & dnb_coordinates,
                               map<unsigned long long int, vector<unsigned int>> & bin_gene_exp_map,
                               const Point& center_point,
                               unsigned short area);

    static unsigned short calcMaxCountOfGeneExp(vector<GeneExpData> & gene_exps);

    void storeAttr(CellBinAttr& cell_attr) const;
    void storeCell();
    void storeGeneList(vector<string>& geneList) const;
    void storeGeneList(GeneData * gene) const;
    void storeCellExp();
    void storeCellBorder(char* borderPath, unsigned int cell_num) const;
    void storeCellTypeList();
    void storeGeneExp();
    void storeGeneAndGeneExp(const vector<string> &gene_name_list);
};


#endif //GEFTOOLS__CELL_BIN_H_
