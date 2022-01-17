/** @file common_bin.h
    @brief Declare a BgefReader class for reading common bin GEF.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__COMMON_BIN_H_
#define GEFTOOLS__COMMON_BIN_H_

#include <numeric>
#include <vector>
#include <map>
#include "hdf5.h"
#include "utils.h"
#include "gef.h"

class BgefReader {
  private:
    int bin_size_ = 0;
    unsigned short gene_num_ = 0;
    unsigned int cell_num_ = 0;
    vector<Coordinate> cell_pos_;
    unsigned int * cell_indices_ = nullptr;
    unsigned int expression_num_ = 0;
    ExpressionAttr expression_attr_{};
    bool expression_attr_init_ = false;
    unsigned int whole_exp_matrix_shape_[2] = {0};
    Gene* genes_ = nullptr;
    Expression* expressions_ = nullptr;
    Mat whole_exp_matrix_t_;
    int version_{};
    int verbose_ = true;

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
    Gene *getGene();
    void buildCellInfo();
    void buildCellInfo2();

  public:
    BgefReader(const string &filename, int bin_size, bool verbose = false);
    virtual ~BgefReader();
    int getVersion() const;
    int getBinSize() const;
    unsigned short getGeneNum() const;
    unsigned int getCellNum();

    /**
     * @brief Get the number of expression.
     */
    unsigned int getExpressionNum() const;
    ExpressionAttr &getExpressionAttr();

    /**
     * @brief Get the shape of wholeExp matrix.
     * @return [rows, cols]
     */
    const unsigned int *getWholeExpMatrixShape() const;

    /**
     * @brief Get gene name list.
     *
     * @param gene_list
     */
    void getGeneNameList(vector<string> & gene_list);

    /**
     *  @brief Get cell name list.
     * @param cell_list
     */
    void getCellPosList(unsigned long long int * cell_list);

    Expression * getExpression();

    void cacheWholeExpMatrix();

    Mat getWholeExpMatrix(Rect roi);

    /**
     * @brief Read WholeExp data to matrix.
     * @param offset_x The starting position on the x-axis to be read.
     * @param offset_y The starting position on the y-axis to be read.
     * @param rows    Number of rows to read.
     * @param cols    Number of cols to read.
     * @param key     MIDcount or genecount.
     * @param matrix  When the value is greater than 255, it will be set to 255.
     */
    void readWholeExpMatrix(unsigned int offset_x,
                           unsigned int offset_y,
                           unsigned int rows,
                           unsigned int cols,
                           string & key,
                           unsigned char *matrix) const;

    /**
     * @brief Read WholeExp data to matrix.
     * @param key     MIDcount or genecount.
     * @param matrix  When the value is greater than 255, it will be set to 255.
     */
    void readWholeExpMatrix(string & key,
                            unsigned char *matrix) const;


    /**
     * Gets sparse matrix indexes
     * @param cell_ind
     * @param count
     * @return
     */
    vector<unsigned long long int> getSparseMatrixIndicesOfExp(unsigned int * cell_ind, unsigned int * count);

    /**
     *
     * @param gene_ind
     * @param gene_names
     */
    void getSparseMatrixIndicesOfGene(unsigned int * gene_ind, char * gene_names);
    //vector<string> getSparseMatrixIndicesOfGene(unsigned int * gene_index);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * Examples:
     * @code
     * # Python
     * from scipy import sparse
     * sparse.csr_matrix((data, indices, indptr), shape=(cell_num, gene_num))
     * @endcode
     * @param indices  CSR format index array of the matrix. Cell id array, the column indices,
     * is the same size as count.
     * @param indptr   CSR format index pointer array of the matrix. indptr length = gene_num_ + 1 .
     * @param count    CSR format data array of the matrix. Expression count.
     * @return
     */
    int getSparseMatrixIndices(unsigned int * indices, unsigned int * indptr, unsigned int * count);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * @param cell_ind     CSR format index array of the matrix. same size as count.
     * @param gene_ind     CSR format index array of the matrix. same size as count.
     * @param count        CSR format data array of the matrix. Expression count.
     */
    int getSparseMatrixIndices2(unsigned int * cell_ind, unsigned int * gene_ind, unsigned int * count);

    /**
     * @brief Get geneID and expCount of this gene for each bin
     * @return key is bin id: x << 32 | y, value is a vector of gene_exp compoud value (geneID << 16 | geneExpCount).
     */
    void getBinGeneExpMap(map<unsigned long long int, vector<CellExpData>>& bin_exp_map);

    /**
     * @brief Free memory for cache variables
     */
    void clear();

    static bool expressionComp(const Expression &p1, const Expression &p2);
};

#endif //GEFTOOLS__COMMON_BIN_H_
