/** @file cgef_reader.h
    @brief Declare a CgefReader class for reading cell bin GEFs.

    Created by huangzhibo on 2021/12/27.
*/

#ifndef GEFTOOLS_CGEF_READER_H
#define GEFTOOLS_CGEF_READER_H

#include <vector>
#include <map>
#include <unordered_map>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "utils.h"
#include "gef.h"

class CgefReader {
  public:
    explicit CgefReader(const string &filename, bool verbose = false);
    ~CgefReader();

    unordered_map<string,unsigned short> genename_to_id_;

    unsigned short getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;

    GeneData *getGene();

    CellData *getCell();

    string getGeneName(unsigned short gene_id);
    int getGeneId(string& gene_name);
    GeneData getGeneData(unsigned short gene_id);
    CellData getCellData(unsigned int cell_id);

    /**
     * @brief Gets gene name array, 32 bit for each string.
     * @param gene_list
     */
    void getGeneNameList(char *gene_list);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * Examples:
     * @code
     * # Python
     * from scipy import sparse
     * sparse.csr_matrix((data, indices, indptr), shape=(M, N))
     * @endcode
     * @param indices  CSR format index array of the matrix. Cell id or Gene id array, same size as count.
     * @param indptr   CSR format index pointer array of the matrix.
     * @param count    CSR format data array of the matrix. Expression count.
     * @param order    Order of count, "gene" or "cell".
     * @return
     */
    int getSparseMatrixIndicesOfExp(unsigned int * indices,
                                    unsigned int * indptr,
                                    unsigned int * count,
                                    const char * order);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * @param cell_ind     CSR format index array of the matrix. same size as count.
     * @param gene_ind     CSR format index array of the matrix. same size as count.
     * @param count        CSR format data array of the matrix. Expression count.
     */
    int getSparseMatrixIndicesOfExp2(unsigned int * cell_ind,
                                     unsigned int * gene_ind,
                                     unsigned int * count);

    /**
     * @brief Gets cellId and count from the geneExp dataset.
     * @param cell_id  Cell id.
     * @param count    Gene expression count.
     */
    void getCellIdAndCount(unsigned int *cell_id, unsigned int *count) const;

    /**
     * @brief Gets geneId and count from the cellExp dataset.
     * @param gene_id  Gene id.
     * @param count    Gene expression count.
     */
    void getGeneIdAndCount(unsigned int *gene_id, unsigned int *count) const;

    /**
     *
     * @param gene_name
     * @return
     */
    int getExpressionCountByGene(string& gene_name, GeneExpData* expressions);

    int getExpressionCountByGenes(vector<string> & gene_names, GeneExpData* expressions);

    unsigned int getExpressionCountByGeneId(unsigned short gene_id, GeneExpData* expressions);

    GeneData getGeneDataByGeneId(unsigned short gene_id);

    void getGeneExpByOffset(unsigned int offset, unsigned int cell_count, GeneExpData *expressions) const;

    /**
     * @brief Output gene expression info to gem format.
     *
     * GEM format:
     * # geneName  x  y  MIDcount  cellID
     * @param filename  The output filename. If it's "stdout", the contents will be printed to standard output.
     * @param gene_name_list List of genes to include.
     * @param force_genes Only warn about unknown subset genes, default: false.
     * @param exclude The genes in gene_name_list is to exclude, default: false.
     * @return Number of output entries.
     */
    unsigned int toGem(string & filename,
                       vector<string> & gene_name_list,
                       bool force_genes = false,
                       bool exclude = false);


  private:
    hid_t file_id_;
    hid_t group_id_;
    hid_t str32_type_;
    hid_t cell_dataset_id_;
    hid_t cell_exp_dataset_id_;
    hid_t gene_dataset_id_;
    hid_t gene_exp_dataset_id_;
    hid_t gene_exp_dataspace_id_;

    GeneData* gene_array_ = nullptr;
    CellData* cell_array_ = nullptr;
    unsigned short gene_num_ = 0;
    unsigned int cell_num_ = 0;
    unsigned int expression_num_ = 0;

    void openCellDataset();
    void openCellExpDataset();
    void openGeneDataset();
    void openGeneExpDataset();

    bool verbose_ = false;
  public:
    bool isVerbose() const;

    void setVerbose(bool verbose);

};

#endif //GEFTOOLS_CGEF_READER_H
