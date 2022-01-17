//
// Created by huangzhibo on 2021/12/15.
//

#include "bgef_reader.h"
#include "khash.h"

KHASH_MAP_INIT_INT64(m64, unsigned int)

BgefReader::BgefReader(const string &filename, int bin_size, bool verbose) {
    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bin_size_ = bin_size;
    verbose_ = verbose;

    hid_t attr;
    attr = H5Aopen(file_id_, "version", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &version_);
    H5Aclose(attr);

    openExpressionSpace();
    openGeneSpace();
    openWholeExpSpace();
}

BgefReader::~BgefReader() {
    if(genes_ != nullptr)
        free(genes_);
    if(cell_indices_ != nullptr)
        free(cell_indices_);
    if(expressions_ != nullptr)
        free(expressions_);
    H5Fclose(file_id_);
    H5Dclose(exp_dataset_id_);
    H5Sclose(exp_dataspace_id_);
    H5Dclose(gene_dataset_id_);
    H5Sclose(gene_dataspace_id_);
    H5Dclose(whole_exp_dataset_id_);
    H5Sclose(whole_exp_dataspace_id_);
}

void BgefReader::openExpressionSpace() {
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
    H5Sget_simple_extent_dims(exp_dataspace_id_, dims, nullptr);
    expression_num_ = dims[0];
}

void BgefReader::openGeneSpace() {
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
    H5Sget_simple_extent_dims(gene_dataspace_id_, dims, nullptr);
    gene_num_ = dims[0];
}

void BgefReader::openWholeExpSpace() {
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
    H5Sget_simple_extent_dims(whole_exp_dataspace_id_, dims, nullptr);

    whole_exp_matrix_shape_[0] = dims[0];
    whole_exp_matrix_shape_[1] = dims[1];
}


void BgefReader::buildCellInfo() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return;

    hid_t memtype;

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Coordinate));
    H5Tinsert(memtype, "x", HOFFSET(Coordinate, pos[0]), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Coordinate, pos[1]), H5T_NATIVE_UINT);

    Coordinate * xy_id;
    xy_id = (Coordinate *) malloc(expression_num_ * sizeof(Coordinate));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, xy_id);

    Coordinate uniq_cell_id{};
    unsigned int index = 0;

    cell_indices_ = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));

    int absent, is_missing;
    khint_t k;
    khash_t(m64) *h = kh_init(m64);  // allocate a hash table

    for (int i = 0; i < expression_num_; ++i) {
        uniq_cell_id = xy_id[i];
        k = kh_get(m64, h, uniq_cell_id.pos_id);
        is_missing = (k == kh_end(h));
        if (!is_missing){
            cell_indices_[i] = kh_value(h, k);
        }else {
            cell_indices_[i] = index;
            cell_pos_.emplace_back(uniq_cell_id);
            k = kh_put(m64, h, uniq_cell_id.pos_id, &absent);  // insert a key to the hash table
            kh_value(h, k) = index;
            ++index;
        }
    }

    cell_num_ = index;
    kh_destroy(m64, h);
    H5Tclose(memtype);
    free(xy_id);
    if(verbose_) printCpuTime(cprev, "buildCellInfo");
}

void BgefReader::buildCellInfo2() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return;

    hid_t memtype;
    Coordinate * xy_id;

    unsigned long cprev2=clock();
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Coordinate));
    H5Tinsert(memtype, "x", HOFFSET(Coordinate, pos[0]), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Coordinate, pos[1]), H5T_NATIVE_UINT);

    xy_id = (Coordinate *) malloc(expression_num_ * sizeof(Coordinate));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, xy_id);
    if(verbose_) printCpuTime(cprev2, "read");


    cell_indices_ = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));

    auto * exp_index = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));
    iota(exp_index, exp_index+expression_num_, 0);
    sort(exp_index, exp_index+expression_num_,
         [xy_id](int a,int b){return xy_id[a].pos_id < xy_id[b].pos_id; });

    Coordinate uniq_cell_id{}, pre_xy = xy_id[exp_index[0]];
    unsigned int index = 0;
    cell_indices_[0] = 0;
    for (int i = 1; i < expression_num_; ++i) {
        uniq_cell_id = xy_id[exp_index[i]];
        if (uniq_cell_id.pos_id != pre_xy.pos_id){
            cell_pos_.emplace_back(pre_xy);
            ++index;
            pre_xy = uniq_cell_id;
        }
        cell_indices_[i] = index;
    }

    cell_pos_.emplace_back(uniq_cell_id);
    ++index;

    cell_num_ = index;
    H5Tclose(memtype);
    free(xy_id);
    if(verbose_) printCpuTime(cprev, "buildCellInfo2");
}


int BgefReader::getBinSize() const {
    return bin_size_;
}

unsigned short BgefReader::getGeneNum() const {
    return gene_num_;
}

unsigned int BgefReader::getCellNum() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return cell_num_;

    buildCellInfo2();
    if(verbose_) printCpuTime(cprev, "getCellNum");
    return cell_num_;
}

unsigned int BgefReader::getExpressionNum() const {
    return expression_num_;
}

ExpressionAttr &BgefReader::getExpressionAttr() {
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

vector<unsigned long long> BgefReader::getSparseMatrixIndicesOfExp(unsigned int * cell_ind, unsigned int * count)
{
    unsigned long cprev=clock();
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
            cell_ind[i] = kh_value(h, k);
        }else {
            cell_ind[i] = index;
            uniq_cells.push_back(uniq_cell_id);
            k = kh_put(m64, h, uniq_cell_id, &absent);  // insert a key to the hash table
            kh_value(h, k) = index;
            ++index;
        }

        count[i] = expData[i].cnt;
    }

    cell_num_ = index;
    kh_destroy(m64, h);
    if(verbose_) printCpuTime(cprev, "getSparseMatrixIndicesOfExp");
    return uniq_cells;
}

//vector<string> BgefReader::getSparseMatrixIndicesOfGene(unsigned int *gene_index) {
//    Gene* gene_data = getGene();
//
//    vector<string> uniq_genes;
//    unsigned long long exp_len_index = 0;
//    for (unsigned short i = 0; i < gene_num_; ++i)
//    {
//        const char* gene = gene_data[i].gene;
//        uniq_genes.emplace_back(gene);
//        unsigned int c = gene_data[i].count;
//        for (int j = 0; j < c; ++j)
//        {
//            gene_index[exp_len_index++] = i;
//        }
//    }
//
//    assert(exp_len_index == expression_num_);
//
//    return uniq_genes;
//}


void BgefReader::getSparseMatrixIndicesOfGene(unsigned int *gene_ind, char * gene_names) {
    Gene* gene_data = getGene();

    unsigned long long exp_len_index = 0;
    for (unsigned short i = 0; i < gene_num_; ++i)
    {
        memcpy(&gene_names[i*32], gene_data[i].gene, 32);
        unsigned int c = gene_data[i].count;
        for (int j = 0; j < c; ++j)
        {
            gene_ind[exp_len_index++] = i;
        }
    }
}

Expression *BgefReader::getExpression() {
    if(expressions_ != nullptr)
        return expressions_;

    hid_t memtype;

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, cnt), H5T_NATIVE_UINT);

    expressions_ = (Expression *) malloc(expression_num_ * sizeof(Expression));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expressions_);

    H5Tclose(memtype);
    return expressions_;
}

Gene *BgefReader::getGene() {
    if(genes_ != nullptr)
        return genes_;

    hid_t memtype, strtype;

    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);

    genes_ = (Gene*)malloc(gene_num_ * sizeof(Gene));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genes_);
    H5Tclose(strtype);
    H5Tclose(memtype);
    return genes_;
}

Mat BgefReader::getWholeExpMatrix(Rect roi){
    if(whole_exp_matrix_t_.empty())
        cacheWholeExpMatrix();
    return whole_exp_matrix_t_(roi);
}

void BgefReader::getBinGeneExpMap(map<unsigned long long int, vector<CellExpData>>& bin_exp_map) {
    Gene* gene_data = getGene();
    Expression* expression_data = getExpression();
//    auto * exp_index = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));
//    iota(exp_index, exp_index+expression_num_, 0);
//
//    sort(exp_index, exp_index+expression_num_,
//         [expression_data](int a,int b){
//        return expression_data[a].x < expression_data[b].x ||
//        (expression_data[a].x == expression_data[b].x && expression_data[a].y < expression_data[b].y); });

    unsigned long long int bin_id, exp_len_index = 0;

    for (unsigned short i = 0; i < gene_num_; ++i)
    {
        assert(exp_len_index + gene_data[i].count <= expression_num_);
        for (unsigned int j = 0; j < gene_data[i].count; ++j)
        {
            Expression exp = expression_data[exp_len_index];
            bin_id = static_cast<unsigned long long int>(exp.x);
            bin_id = bin_id << 32 | exp.y;

            CellExpData gene_exp = {i, static_cast<unsigned short>(exp.cnt)};

            exp_len_index++;

            auto iter = bin_exp_map.find(bin_id);
            if(iter != bin_exp_map.end()){
                iter->second.emplace_back(gene_exp);
            }else{
                vector<CellExpData> gene_exp_list;
                gene_exp_list.emplace_back(gene_exp);
                bin_exp_map.insert(map<unsigned long long, vector<CellExpData>>::value_type (bin_id, gene_exp_list));
            }
        }
    }
}

void BgefReader::clear() {
    if(genes_ != nullptr)
        free(genes_);
    if(expressions_ != nullptr)
        free(expressions_);
    whole_exp_matrix_t_.release();
}

void BgefReader::cacheWholeExpMatrix() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, 1);
    // genecount的值大于255将读取为255
    whole_exp_matrix_t_ = Mat::zeros(
            (int)whole_exp_matrix_shape_[0], (int)whole_exp_matrix_shape_[1], CV_8UC1);
    H5Tinsert(memtype, "genecount", 0, H5T_NATIVE_UCHAR);
    H5Dread(whole_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, whole_exp_matrix_t_.data);
    whole_exp_matrix_t_ = whole_exp_matrix_t_.t();
    H5Tclose(memtype);
}

void BgefReader::readWholeExpMatrix(string &key, unsigned char *matrix) const {
    readWholeExpMatrix(0,
                       0,
                       whole_exp_matrix_shape_[0],
                       whole_exp_matrix_shape_[1],
                       key,
                       matrix);
}

void BgefReader::readWholeExpMatrix(unsigned int offset_x, unsigned int offset_y, unsigned int rows, unsigned int cols,
                                    string &key, unsigned char *matrix) const {

    hsize_t start[2] = {offset_x, offset_y},
            count[2] = {rows, cols},
            offset_out[2] = {0, 0};

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, 1);
    // genecount的值大于255将读取为255
    H5Tinsert(memtype, key.c_str(), 0, H5T_NATIVE_UCHAR);

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(2, count,nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(whole_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(whole_exp_dataset_id_, memtype, memspace, whole_exp_dataspace_id_, H5P_DEFAULT, matrix);

    H5Tclose(memtype);
    H5Sclose(memspace);
}

int BgefReader::getVersion() const {
    return version_;
}

void BgefReader::getGeneNameList(vector<string> & gene_list) {
    Gene * genes = getGene();
    for(unsigned int i = 0; i < gene_num_; i++){
        string name = genes[i].gene;
        gene_list.emplace_back(name);
    }
}

//void BgefReader::getGeneNameList(char *gene_list) {
//    Gene * genes = getGene();
//    for(unsigned int i = 0; i < gene_num_; i++){
//        memcpy(&gene_list[i], genes[i].gene, 32);
//    }
//}


bool BgefReader::expressionComp(const Expression& p1, const Expression& p2) {
    return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
}

const unsigned int *BgefReader::getWholeExpMatrixShape() const {
    return whole_exp_matrix_shape_;
}


int BgefReader::getSparseMatrixIndices(unsigned int *indices, unsigned int *indptr, unsigned int *count) {
    unsigned long cprev=clock();

    if(cell_indices_ == nullptr) buildCellInfo2();
    memcpy(indices, cell_indices_, expression_num_ * sizeof(unsigned int));

    Gene * gene_data = getGene();

    indptr[0] = 0;
    for(unsigned int i = 1; i < gene_num_; i++){
        indptr[i] = gene_data[i].offset;
    }
    indptr[cell_num_] = gene_data[gene_num_-1].offset + gene_data[gene_num_-1].count;

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "count", 0, H5T_NATIVE_UINT);
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);

    if(verbose_) printCpuTime(cprev, "getSparseMatrixIndices");
    return 0;
}

int BgefReader::getSparseMatrixIndices2(unsigned int * cell_ind, unsigned int * gene_ind, unsigned int * count){
    unsigned long cprev=clock();
    Gene * gene_data = getGene();

    if(cell_indices_ == nullptr) buildCellInfo2();
    memcpy(cell_ind, cell_indices_, expression_num_ * sizeof(unsigned int));

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "count", 0, H5T_NATIVE_UINT);
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);

    unsigned int n = 0;
    for(unsigned short i = 0; i < gene_num_; i++){
        for(unsigned int j = 0; j < gene_data[i].count; j++){
            gene_ind[n++] = i;
        }
    }

    H5Tclose(memtype);
    if(verbose_) printCpuTime(cprev, "getSparseMatrixIndices2");
    return 0;
}

void BgefReader::getCellPosList(unsigned long long int *cell_list) {
    memcpy(cell_list, cell_pos_.data(), cell_num_ * sizeof(unsigned long long int));
}



