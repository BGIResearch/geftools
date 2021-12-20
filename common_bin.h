/** @file common_bin.h
    @brief Declare a CommonBin class and several related structs.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__COMMON_BIN_H_
#define GEFTOOLS__COMMON_BIN_H_

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "hdf5.h"
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

/**
 * @brief Expression struct
 */
struct Expression {
    unsigned int x; ///< dnb coordinates x
    unsigned int y; ///< dnb coordinates x
    unsigned int cnt; ///< expression count (MIDcount)
};

/**
 * @brief ExpressionAttr struct, record the attributes of the dataset named expressione
 */
struct ExpressionAttr
{
    unsigned int min_x; ///< Min X of dnb coordinate
    unsigned int min_y; ///< Min Y of dnb coordinate
    unsigned int max_x; ///< Max X of dnb coordinate
    unsigned int max_y; ///< Max Y of dnb coordinate
    unsigned int max_exp;  ///< Max expression count
    unsigned int resolution; ///< The resolution of stereo chip
};

struct Gene {
    Gene(const char* g, unsigned int o, unsigned c)
    {
        int i = 0;
        while (g[i] != '\0')
        {
            gene[i] = g[i];
            ++i;
        }
        offset = o;
        count = c;
    }
    char gene[32] = {0};
    unsigned int offset;
    unsigned int count;
};

// wholeExp Matrix
struct DnbStat {
    unsigned int MIDcount;
    unsigned short genecount;
};


class CommonBin {
  private:
    int bin_size_ = 0;
    unsigned int gene_num_ = 0;
    unsigned int cell_num_ = 0;
    unsigned long long expression_num_ = 0;
    ExpressionAttr expression_attr_{};
    bool expression_attr_init_ = false;
    unsigned int dnb_stat_matrix_shape_[2];
    Gene* genes_ = nullptr;
    Expression* expressions_ = nullptr;
    Mat whole_exp_matrix_t_;
    int version_;

    hid_t file_id_;
    hid_t exp_dataspace_id_{};
    hid_t exp_dataset_id_{};
    hid_t gene_dataspace_id_{};
    hid_t gene_dataset_id_{};
    hid_t whole_exp_dataspace_id_{};
    hid_t whole_exp_dataset_id_{};

    void openExpressionSpace();
    void openGeneSpace();
    void openWholeExpSpace();


  public:
    const unsigned int *getDnbStatMatrixShape() const;
    CommonBin(const string &filename, int bin_size);
    virtual ~CommonBin();
    int getVersion() const;
    int getBinSize() const;
    unsigned int getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;
    ExpressionAttr &getExpressionAttr();

    Gene *getGene();
    Gene *getGeneIndexes();
    Expression * getExpression();

    void cacheWholeExpMatrix();

    Mat getWholeExpMatrix(Rect roi);

    void getDnbStatMatrix(unsigned int offset_x,
                          unsigned int offset_y,
                          unsigned int rows,
                          unsigned int cols,
                          unsigned char * matrix) const;

    void getDnbStatMatrixT(unsigned int offset_x,
                          unsigned int offset_y,
                          unsigned int rows,
                          unsigned int cols,
                          unsigned char * matrix) const;

    //sparse matrix indexes
    vector<unsigned long long int> getSparseMatrixIndexesOfExp(unsigned int * cell_index, unsigned int * count);
    vector<string> getSparseMatrixIndexesOfGene(unsigned int * gene_index);

    /**
     * @brief Get geneID and expCount of this gene for each bin
     * @return key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     */
    map<unsigned long long int, vector<unsigned int>> getBinGeneExpMap();

    /**
     * @brief Free memory for cache variables
     */
    void clear();
};

#endif //GEFTOOLS__COMMON_BIN_H_
