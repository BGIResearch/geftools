//
// Created by huangzhibo on 2021/12/15.
//

#ifndef GEFTOOLS__COMMON_BIN_H_
#define GEFTOOLS__COMMON_BIN_H_

#include <string>
#include <vector>
#include <iostream>
#include "hdf5.h"

using namespace std;

struct Expression {
    unsigned int x;
    unsigned int y;
    unsigned int cnt;
};

struct ExpressionAttr
{
    unsigned int min_x;
    unsigned int min_y;
    unsigned int max_x;
    unsigned int max_y;
    unsigned int max_exp;
    unsigned int resolution;
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
    int getBinSize() const;
    unsigned int getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;
    ExpressionAttr &getExpressionAttr();

    Gene *getGene() const;
    Expression * getExpression() const;

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
    vector<unsigned long long> getSparseMatrixIndexesOfExp(unsigned int * cell_index, unsigned int * count);
    vector<string> getSparseMatrixIndexesOfGene(unsigned int * gene_index) const;
};

#endif //GEFTOOLS__COMMON_BIN_H_
