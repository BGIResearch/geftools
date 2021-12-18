//
// Created by huangzhibo on 2021/12/15.
//

#include "common_bin.h"
#include "khash.h"

KHASH_MAP_INIT_INT64(m64, unsigned int)

CommonBin::CommonBin(const string &filename, int bin_size) {
    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bin_size_ = bin_size;
    openExpressionSpace();
    openGeneSpace();
    openWholeExpSpace();
}

void CommonBin::openExpressionSpace() {
    hsize_t dims[1];
    // Read raw data
    char expName[128]={0};
    sprintf(expName, "/geneExp/bin%d/expression", bin_size_);
    exp_dataset_id_ = H5Dopen(file_id_, expName, H5P_DEFAULT);
    if (exp_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<expName<<endl;
        return;
    }
    exp_dataspace_id_ = H5Dget_space(exp_dataset_id_);
    H5Sget_simple_extent_dims(exp_dataspace_id_, dims, NULL);
    expression_num_ = dims[0];
}

void CommonBin::openGeneSpace() {
    hsize_t dims[1];

    // Read index
    char idxName[128]={0};
    sprintf(idxName, "/geneExp/bin%d/gene", bin_size_);
    gene_dataset_id_ = H5Dopen(file_id_, idxName, H5P_DEFAULT);
    if (gene_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<idxName<<endl;
        return;
    }
    gene_dataspace_id_ = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id_, dims, NULL);
    gene_num_ = dims[0];
}

void CommonBin::openWholeExpSpace() {
    hsize_t dims[2];

    // Read index
    char idxName[128]={0};
    sprintf(idxName, "/wholeExp/bin%d", bin_size_);
    whole_exp_dataset_id_ = H5Dopen(file_id_, idxName, H5P_DEFAULT);
    if (whole_exp_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<idxName<<endl;
        return;
    }
    whole_exp_dataspace_id_ = H5Dget_space(whole_exp_dataset_id_);
    H5Sget_simple_extent_dims(whole_exp_dataspace_id_, dims, NULL);

    dnb_stat_matrix_shape_[0] = dims[0];
    dnb_stat_matrix_shape_[1] = dims[1];
}

CommonBin::~CommonBin() {
    H5Fclose(file_id_);
    H5Dclose(exp_dataset_id_);
    H5Sclose(exp_dataspace_id_);
    H5Dclose(gene_dataset_id_);
    H5Sclose(gene_dataspace_id_);
    H5Dclose(whole_exp_dataset_id_);
    H5Sclose(whole_exp_dataspace_id_);
}

int CommonBin::getBinSize() const {
    return bin_size_;
}

unsigned int CommonBin::getGeneNum() const {
    return gene_num_;
}

unsigned int CommonBin::getCellNum() const {
    return cell_num_;
}

unsigned long long int CommonBin::getExpressionNum() const {
    return expression_num_;
}

ExpressionAttr &CommonBin::getExpressionAttr() {
    if(expression_attr_init_)
        return expression_attr_;
    hid_t attr;
    attr = H5Aopen(exp_dataset_id_, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.min_x));
    attr = H5Aopen(exp_dataset_id_, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.min_y));
    attr = H5Aopen(exp_dataset_id_, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_x));
    attr = H5Aopen(exp_dataset_id_, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_y));
    attr = H5Aopen(exp_dataset_id_, "maxExp", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_exp));
    attr = H5Aopen(exp_dataset_id_, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.resolution));
    expression_attr_init_ = true;

    H5Aclose(attr);
    return expression_attr_;
}

vector<unsigned long long> CommonBin::getSparseMatrixIndexesOfExp(unsigned int * cell_index, unsigned int * count)
{
    unsigned long long uniq_cell_id;
    Expression* expData = getExpression();

    vector<unsigned long long> uniq_cells;
    unsigned int index = 0;

    int absent, is_missing;
    khint_t k;
    khash_t(m64) *h = kh_init(m64);  // allocate a hash table

    for (int i = 0; i < expression_num_; ++i) {
        uniq_cell_id = expData[i].x;
        uniq_cell_id = uniq_cell_id << 32 | expData[i].y;

        k = kh_get(m64, h, uniq_cell_id);
        is_missing = (k == kh_end(h));
        if (!is_missing){
            cell_index[i] = kh_value(h, k);
        }else {
            cell_index[i] = index;
            uniq_cells.push_back(uniq_cell_id);
            k = kh_put(m64, h, uniq_cell_id, &absent);  // insert a key to the hash table
            kh_value(h, k) = index;
            ++index;
        }

        count[i] = expData[i].cnt;
    }

    cell_num_ = index;

    if (expData != nullptr)
        free(expData);

    kh_destroy(m64, h);

    return uniq_cells;
}

vector<string> CommonBin::getSparseMatrixIndexesOfGene(unsigned int *gene_index) const {
    Gene* geneData = getGene();

    vector<string> uniq_genes;
    unsigned long long exp_len_index = 0;
    for (unsigned int i = 0; i < gene_num_; ++i)
    {
        const char* gene = geneData[i].gene;
        uniq_genes.emplace_back(gene);
        unsigned int c = geneData[i].count;
        for (int j = 0; j < c; ++j)
        {
            gene_index[exp_len_index++] = i;
        }
    }

    assert(exp_len_index == expression_num_);

    if (geneData != nullptr)
        free(geneData);

    return uniq_genes;
}

Expression *CommonBin::getExpression() const {
    hid_t memtype;

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, cnt), H5T_NATIVE_UINT);

    Expression *expData;
    expData = (Expression *) malloc(expression_num_ * sizeof(Expression));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expData);

    H5Tclose(memtype);
    return expData;
}

Gene *CommonBin::getGene() const {
    hid_t memtype, strtype;

    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);

    Gene* geneData;
    geneData = (Gene*)malloc(gene_num_ * sizeof(Gene));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, geneData);
    H5Tclose(strtype);
    H5Tclose(memtype);
    return geneData;
}

void CommonBin::getDnbStatMatrix(unsigned int offset_x,
                                 unsigned int offset_y,
                                 unsigned int rows,
                                 unsigned int cols,
                                 unsigned char *matrix) const {
    hsize_t start[2] = {offset_x, offset_y},
            count[2] = {rows, cols},
            offset_out[2] = {0, 0};

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, 1);
    // genecount的值大于255将读取为255
    H5Tinsert(memtype, "genecount", 0, H5T_NATIVE_UCHAR);

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(2, count,NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL);

    H5Sselect_hyperslab (whole_exp_dataspace_id_, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dread (whole_exp_dataset_id_, memtype, memspace, whole_exp_dataspace_id_, H5P_DEFAULT, matrix);
}

const unsigned int *CommonBin::getDnbStatMatrixShape() const {
    return dnb_stat_matrix_shape_;
}

void CommonBin::getDnbStatMatrixT(unsigned int offset_x,
                                  unsigned int offset_y,
                                  unsigned int rows,
                                  unsigned int cols,
                                  unsigned char *matrix) const {
    getDnbStatMatrix(offset_y, offset_x, cols, rows, matrix);
}

