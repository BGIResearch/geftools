
#include "cgef_reader.h"

CgefReader::CgefReader(const string &filename, bool verbose) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);
    verbose_ = verbose;

    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id_ = H5Gopen(file_id_, "/cellBin", H5P_DEFAULT);
    cell_dataset_id_ = openCellDataset(group_id_);
    cell_exp_dataset_id_ = openCellExpDataset(group_id_);
    gene_dataset_id_ = openGeneDataset(group_id_);
    gene_exp_dataset_id_ = openGeneExpDataset(group_id_);
    gene_exp_dataspace_id_ = H5Dget_space(gene_exp_dataset_id_);

    hsize_t dims[1];
    cell_exp_dataspace_id_ = H5Dget_space(cell_exp_dataset_id_);
    H5Sget_simple_extent_dims(cell_exp_dataspace_id_, dims, nullptr);
    expression_num_ = dims[0];
    expression_num_current_ = dims[0];

    cell_dataspace_id_ = H5Dget_space(cell_dataset_id_);
    H5Sget_simple_extent_dims(cell_dataspace_id_, dims, nullptr);
    cell_num_ = dims[0];
    cell_num_current_ = dims[0];

    gene_array_ = loadGene();
}

CgefReader::~CgefReader() {
    H5Tclose(str32_type_);
    H5Fclose(file_id_);
    H5Gclose(group_id_);
    H5Dclose(cell_dataset_id_);
    H5Dclose(gene_dataset_id_);
    H5Dclose(cell_exp_dataset_id_);
    H5Dclose(gene_exp_dataset_id_);
    H5Sclose(cell_dataspace_id_);
    H5Sclose(cell_exp_dataspace_id_);
    H5Sclose(gene_exp_dataspace_id_);
    free(gene_array_);
    if (cell_array_ != nullptr)
        free(cell_array_);
    if (cell_array_current_ != nullptr) free(cell_array_current_);
    if (cell_id_array_current_ != nullptr) free(cell_id_array_current_);
}

hid_t CgefReader::openCellDataset(hid_t group_id) {
    cell_dataset_id_ = H5Dopen(group_id, "cell", H5P_DEFAULT);
    if (cell_dataset_id_ < 0) {
        cerr << "failed open dataset: cell" << endl;
        exit(3);
    }

    hsize_t dims_attr[1];
    hid_t attr, attr_dataspace;
    attr = H5Aopen(cell_dataset_id_, "blockIndex", H5P_DEFAULT);
    attr_dataspace = H5Aget_space(attr);
    H5Sget_simple_extent_dims(attr_dataspace, dims_attr, nullptr);

    block_index_ = static_cast<unsigned int *>(
            malloc(dims_attr[0] * sizeof(unsigned int)));

    H5Aread(attr, H5T_NATIVE_UINT32, block_index_);

    attr = H5Aopen(cell_dataset_id_, "blockSize", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT32, block_size_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);

    return cell_dataset_id_;
}

hid_t CgefReader::openGeneDataset(hid_t group_id) {
    hsize_t dims[1];
    gene_dataset_id_ = H5Dopen(group_id, "gene", H5P_DEFAULT);
    if (gene_dataset_id_ < 0) {
        cerr << "failed open dataset: gene" << endl;
        return gene_dataset_id_;
    }
    hid_t gene_dataspace_id = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id, dims, nullptr);
    gene_num_ = dims[0];
    gene_num_current_ = dims[0];
    H5Sclose(gene_dataspace_id);
    return gene_dataset_id_;
}

hid_t CgefReader::openCellExpDataset(hid_t group_id) {
    cell_exp_dataset_id_ = H5Dopen(group_id, "cellExp", H5P_DEFAULT);
    if (cell_exp_dataset_id_ < 0) {
        cerr << "failed open dataset: cellExp" << endl;
        exit(3);
    }
    return cell_exp_dataset_id_;
}

hid_t CgefReader::openGeneExpDataset(hid_t group_id) {
    gene_exp_dataset_id_ = H5Dopen(group_id, "geneExp", H5P_DEFAULT);
    if (gene_exp_dataset_id_ < 0) {
        cerr << "failed open dataset: geneExp" << endl;
        return gene_exp_dataset_id_;
    }
    return gene_exp_dataset_id_;
}

unsigned int CgefReader::getCellNum() const {
    return cell_num_current_;
}

unsigned short CgefReader::getGeneNum() const {
    return gene_num_current_;
}

unsigned long long int CgefReader::getExpressionNum() const {
    return expression_num_current_;
}

GeneData *CgefReader::loadGene(bool reload) {
    unsigned long cprev = clock();
    if (gene_array_ != nullptr) {
        if (reload) free(gene_array_);
        else return gene_array_;
    }

    hid_t memtype = getMemtypeOfGeneData();
    gene_array_ = (GeneData *) malloc(gene_num_ * sizeof(GeneData));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_array_);

    for (unsigned short i = 0; i < gene_num_; i++) {
        genename_to_id_[gene_array_[i].gene_name] = i;
    }

    H5Tclose(memtype);
    if (verbose_) printCpuTime(cprev, "loadGene");
    return gene_array_;
}

CellData *CgefReader::loadCell(bool reload) {
    unsigned long cprev = clock();
    if (cell_array_ != nullptr) {
        if (reload) free(cell_array_);
        else return cell_array_;
    }

    hid_t memtype = getMemtypeOfCellData();
    cell_array_ = (CellData *) malloc(cell_num_ * sizeof(CellData));
    H5Dread(cell_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array_);
    H5Tclose(memtype);
    if (verbose_) printCpuTime(cprev, "getCell");
    return cell_array_;
}

void CgefReader::getGeneNameList(vector<string> &gene_list) {
    for (unsigned int i = 0; i < gene_num_; i++) {
        gene_list.emplace_back(gene_array_[i].gene_name);
//        memcpy(&gene_list[i*32], genes[i].gene_name, 32);
    }
}


void CgefReader::getGeneName(char *gene_list) {
    if (gene_id_set_current_.empty()) {
        for (unsigned int i = 0; i < gene_num_; i++) {
            memcpy(&gene_list[i * 32], gene_array_[i].gene_name, 32);
        }
    } else {
        int i = 0;
        for (auto gene_id: gene_id_set_current_) {
            memcpy(&gene_list[i * 32], gene_array_[gene_id].gene_name, 32);
            i++;
        }
    }
}

void CgefReader::getCellNameList(unsigned long long int *cell_pos_list) {
    if (restrict_region_) {
        for (unsigned int i = 0; i < cell_num_current_; i++) {
            cell_pos_list[i] = cell_array_current_[i].x;
            cell_pos_list[i] = cell_pos_list[i] << 32 | cell_array_current_[i].y;
        }
    } else {
        CellData *cells = loadCell();
        for (unsigned int i = 0; i < cell_num_; i++) {
            cell_pos_list[i] = cells[i].x;
            cell_pos_list[i] = cell_pos_list[i] << 32 | cells[i].y;
        }
    }
}

//indptr length = gene_num_ + 1
int CgefReader::getSparseMatrixIndices(unsigned int *indices, unsigned int *indptr, unsigned int *count,
                                       const char *order) {
    if (order[0] == 'g') {
        if (restrict_region_) {
            cerr << "This method is inefficient when using restrictRegion and (order == gene),"
                    " please use order == cell." << endl;
            exit(2);

//            unordered_map<unsigned int, bool> cellid_map;
//            for(unsigned int i = 0; i < cell_num_current_; i++) {
//                cellid_map[cell_id_array_current_[i]] = true;
//            }
//
//            unsigned int n = 0;
//            indptr[0] = 0;
//
//            auto *gene_exp_data = static_cast<GeneExpData *>(malloc(cell_num_ * sizeof(GeneExpData)));
//            for(unsigned int i = 0; i < gene_num_current_; i++) {
//                unsigned int gene_id = restrict_gene_ ? gene_id_array_current_[i] : i;
//                GeneData gene_data = gene_array_[gene_id];
//                selectGeneExp(gene_data.offset, gene_data.cell_count, gene_exp_data);
//
//                unsigned int c_count = 0;
//                for (unsigned int j = 0; j < gene_data.cell_count; j++) {
//                    unsigned int cid = gene_exp_data[j].cell_id;
//                    if (cellid_map.find(cid) == cellid_map.end()) continue;
//                    indices[n] = cid;
//                    count[n] = gene_exp_data[j].count;
//                    n++;
//                    c_count++;
//                }
//                indptr[i+1] = indptr[i] + c_count;
//            }
//            free(gene_exp_data);
        } else {
            hid_t memtype;
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
            H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "cellID", 0, H5T_NATIVE_UINT);
            H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
            for (unsigned int i = 0; i < gene_num_; i++) {
                indptr[i] = gene_array_[i].offset;
            }
            indptr[gene_num_] = gene_array_[gene_num_ - 1].offset + gene_array_[gene_num_ - 1].cell_count;
            H5Tclose(memtype);
        }

    } else if (order[0] == 'c') {
        if (restrict_region_) {
            auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
            for (unsigned int i = 0; i < cell_num_current_; i++) {
                CellData cell = cell_array_current_[i];
                selectCellExp(cell.offset, cell.gene_count, cell_exp_data);

                indptr[i] = cell.offset;
                for (unsigned int j = 0; j < cell.gene_count; j++) {
                    indices[i + j] = cell_exp_data[j].gene_id;
                    count[i + j] = cell_exp_data[j].count;
                }
            }
            indptr[cell_num_current_] = cell_array_current_[cell_num_current_].offset
                                        + cell_array_current_[cell_num_current_].gene_count;
            free(cell_exp_data);
        } else {
            hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "geneID", 0, H5T_NATIVE_USHORT);
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);

            CellData *cell_data = loadCell();
            indptr[0] = 0;
            for (unsigned int i = 1; i < cell_num_; i++) {
                indptr[i] = cell_data[i].offset;
            }
            indptr[cell_num_] = cell_data[cell_num_ - 1].offset + cell_data[cell_num_ - 1].gene_count;
            H5Tclose(memtype);
        }
    } else {
        return -1;
    }

    return 0;
}

int CgefReader::getSparseMatrixIndices2(unsigned int *cell_ind, unsigned int *gene_ind, unsigned int *count) {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "cellID", 0, H5T_NATIVE_UINT);
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_ind);

    unsigned int n = 0;
    for (unsigned short i = 0; i < gene_num_; i++) {
        unsigned int cell_count = gene_array_[i].cell_count;
        for (unsigned int j = 0; j < cell_count; j++) {
            gene_ind[n++] = i;
        }
    }

    H5Tclose(memtype);
    return 0;
}

void CgefReader::getCellIdAndCount(unsigned int *cell_id, unsigned short *count) const{
    hid_t memtype = getMemtypeOfGeneExpData();
    GeneExpData* gene_exp_data;
    gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);
    for(unsigned int i = 0; i < expression_num_; i++){
        cell_id[i] = gene_exp_data->cell_id;
        count[i] = gene_exp_data->count;
    }
    free(gene_exp_data);
    H5Tclose(memtype);
}

void CgefReader::getGeneIdAndCount(unsigned short *gene_id, unsigned short *count) const{
    hid_t memtype = getMemtypeOfCellExpData();
    CellExpData* cell_exp_data;
    cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
    H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);
    for(unsigned int i = 0; i < expression_num_; i++){
        gene_id[i] = cell_exp_data->gene_id;
        count[i] = cell_exp_data->count;
    }
    free(cell_exp_data);
    H5Tclose(memtype);
}

int CgefReader::getExpressionCountByGene(string &gene_name, GeneExpData *expressions) {
    int gene_id = getGeneId(gene_name);
    if (gene_id < 0) {
        cerr << "Gene ID < 0 : " << gene_id << endl;
        exit(2);
    }
    getExpressionCountByGeneId(gene_id, expressions);
    return 0;
}

unsigned int CgefReader::getExpressionCountByGeneId(unsigned short gene_id, GeneExpData *expressions) {
    selectGeneExp(gene_array_[gene_id].offset, gene_array_[gene_id].cell_count, expressions);
    return gene_array_[gene_id].cell_count;
}

GeneData CgefReader::getGeneDataByGeneId(unsigned short gene_id) {
    return gene_array_[gene_id];
}

CellData CgefReader::getCell(unsigned int cell_id) const {
    CellData cell_data{};
    selectCells(cell_id, 1, &cell_data);
    return cell_data;
}

CellData *CgefReader::getCell() {
    if (cell_array_current_ != nullptr) {
        return cell_array_current_;
    }
    return loadCell();
}

string CgefReader::getGeneName(unsigned short gene_id) {
    return gene_array_[gene_id].gene_name;
}

GeneData CgefReader::getGene(unsigned short gene_id) const {
    return gene_array_[gene_id];
}

GeneData *CgefReader::getGene() {
    if (gene_array_current_ != nullptr) {
        return gene_array_current_;
    } else if (!gene_id_set_current_.empty()) {
        gene_array_current_ = static_cast<GeneData *>(malloc(gene_num_current_ * sizeof(GeneData)));
        int i = 0;
        for (auto gene_id: gene_id_set_current_) {
            memcpy(&gene_array_current_[i], &gene_array_[gene_id], sizeof(GeneData));
        }
        return gene_array_current_;
    }
    return gene_array_;
}

int CgefReader::getGeneId(string &gene_name) {
    auto iter = genename_to_id_.find(gene_name);
    if (iter != genename_to_id_.end()) {
        return iter->second;
    }
    return -1;
}

unsigned int CgefReader::toGem(string &filename, const vector<string> &gene_name_list, bool force_genes, bool exclude) {
    unsigned long cprev = clock();
    unsigned short gene_ids[gene_num_];
    unsigned short n = 0;
    unsigned int cell_num = 0;
    unordered_map<string, bool> genename_map;
    for (const auto &gene_name: gene_name_list) {
        genename_map[gene_name] = true;
    }
    if (exclude) {
        for (unsigned short i = 0; i < gene_num_; i++) {
            if (genename_map.find(gene_array_[i].gene_name) == genename_map.end()) {
                gene_ids[n++] = i;
            }
        }
    } else {
        for (const auto &gene_name: gene_name_list) {
            if (genename_to_id_.find(gene_name) == genename_to_id_.end()) {
                cerr << "Gene ( " << gene_name << " ) does not exist." << endl;
                if (!force_genes) exit(2);
                continue;
            }
            gene_ids[n++] = genename_to_id_[gene_name];
        }
    }

    ofstream fout;

    bool to_file = true;
    if (filename != "stdout") {
        fout.open(filename);
        if (!fout.is_open()) cerr << "Fail to open file : " << filename << endl;
        fout << "#geneName\tx\ty\tcount\tcellID" << endl;
    } else {
        to_file = false;
        cout << "#geneName\tx\ty\tcount\tcellID" << endl;
    }

    if (restrict_region_) {
        if (verbose_) cerr << "toGem restrict_region_ true" << endl;
        auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
        for (unsigned int i = 0; i < cell_num_current_; i++) {
            unsigned int cell_id = cell_id_array_current_[i];
            CellData cell = cell_array_current_[i];
            selectCellExp(cell.offset, cell.gene_count, cell_exp_data);

            for (unsigned int j = 0; j < cell.gene_count; j++) {
                string gene_name = gene_array_[cell_exp_data[j].gene_id].gene_name;

                if (!gene_name_list.empty() && (
                        (exclude && genename_map.find(gene_name) != genename_map.end())
                        || (!exclude && genename_map.find(gene_name) == genename_map.end())
                ))
                    continue;

                cell_num++;
                if (to_file) {
                    fout << gene_name << "\t" << cell.x << "\t" << cell.y
                         << "\t" << cell_exp_data[j].count << "\t" << cell_id << endl;
                } else {
                    cout << gene_name << "\t" << cell.x << "\t" << cell.y
                         << "\t" << cell_exp_data[j].count << "\t" << cell_id << endl;
                }
            }
        }

        free(cell_exp_data);
    } else {
        unsigned int max_malloc = 1000;
        auto *expression = static_cast<GeneExpData *>(malloc(max_malloc * sizeof(GeneExpData)));

        for (unsigned short i = 0; i < n; i++) {
            unsigned short gene_id = gene_ids[i];
            GeneData gene_data = getGeneDataByGeneId(gene_id);
            unsigned int cell_count = gene_data.cell_count;

            cell_num += cell_count;

            if (cell_count > max_malloc) {
                expression = static_cast<GeneExpData *>(realloc(expression, cell_count));
                max_malloc = cell_count + 1000;
            }

            selectGeneExp(gene_data.offset, cell_count, expression);

            for (unsigned int j = 0; j < cell_count; j++) {
                CellData cell_data = getCell(expression[j].cell_id);
                if (to_file) {
                    fout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
                         << "\t" << expression[j].count << "\t" << expression[j].cell_id << endl;
                } else {
                    cout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
                         << "\t" << expression[j].count << "\t" << expression[j].cell_id << endl;
                }
            }
        }
        free(expression);
    }

    if (verbose_) printCpuTime(cprev, "toGem");

    if (to_file) fout.close();
    return cell_num;
}

bool CgefReader::isVerbose() const {
    return verbose_;
}

void CgefReader::setVerbose(bool verbose) {
    verbose_ = verbose;
}

void CgefReader::selectCells(unsigned int offset,
                             unsigned int cell_count,
                             CellData *cell) const {
    hsize_t start[1] = {offset},
            count[1] = {cell_count},
            offset_out[1] = {0};

    hid_t memtype;
    memtype = getMemtypeOfCellData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_dataset_id_, memtype, memspace, cell_dataspace_id_, H5P_DEFAULT, cell);
}

void CgefReader::selectCellExp(unsigned int offset,
                               unsigned int gene_count,
                               CellExpData *cell_exp_data) const {
    hsize_t start[1] = {offset},
            count[1] = {gene_count},
            offset_out[1] = {0};

    hid_t memtype;
    memtype = getMemtypeOfCellExpData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_exp_dataset_id_, memtype, memspace, cell_exp_dataspace_id_, H5P_DEFAULT, cell_exp_data);
}

void CgefReader::selectGeneExp(unsigned int offset,
                               unsigned int cell_count,
                               GeneExpData *gene_exp_data) const {
    hsize_t start[1] = {offset},
            count[1] = {cell_count},
            offset_out[1] = {0};

    hid_t memtype = getMemtypeOfGeneExpData();

    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(gene_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(gene_exp_dataset_id_, memtype, memspace, gene_exp_dataspace_id_, H5P_DEFAULT, gene_exp_data);
}

bool CgefReader::isRestrictRegion() const {
    return restrict_region_;
}

bool CgefReader::isRestrictGene() const {
    return restrict_gene_;
}

void CgefReader::restrictRegion(unsigned int min_x, unsigned int max_x, unsigned int min_y, unsigned int max_y) {
    unsigned long cprev = clock();
    if (restrict_gene_) {
        cerr << "Please call restrictRegion before restrictGene function." << endl;
        exit(2);
    }

    restrict_region_ = true;
    unsigned int x_block_num = block_size_[2];
    unsigned int y_block_num = block_size_[3];
    unsigned int min_block_x = min_x / block_size_[0];
    unsigned int max_block_x = max_x / block_size_[0];
    unsigned int min_block_y = min_y / block_size_[1];
    unsigned int max_block_y = max_y / block_size_[1];

    max_block_x = max_block_x > x_block_num ? x_block_num : max_block_x;
    max_block_y = max_block_y > y_block_num ? y_block_num : max_block_y;

    unsigned int block_id_y, offset, cell_num_tmp = 0;
    for (unsigned int y = min_block_y; y <= max_block_y; y++) {
        block_id_y = y * x_block_num;
        cell_num_tmp += block_index_[max_block_x + block_id_y + 1] - block_index_[min_block_x + block_id_y];
    }

    // support rerun this method
    cell_num_current_ = 0;
    expression_num_current_ = 0;
    if (cell_array_current_ != nullptr) free(cell_array_current_);
    if (cell_id_array_current_ != nullptr) free(cell_id_array_current_);
    if (!gene_id_set_current_.empty()) gene_id_set_current_.clear();

    cell_array_current_ = static_cast<CellData *>(malloc(cell_num_tmp * sizeof(CellData)));
    cell_id_array_current_ = static_cast<unsigned int *>(malloc(cell_num_tmp * sizeof(unsigned int)));
    auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));

    for (unsigned int y = min_block_y; y <= max_block_y; y++) {
        block_id_y = y * x_block_num;
        offset = block_index_[min_block_x + block_id_y];
        cell_num_tmp = block_index_[max_block_x + block_id_y + 1] - offset;
        selectCells(offset, cell_num_tmp, &cell_array_current_[cell_num_current_]);

        for (unsigned int j = 0; j < cell_num_tmp; j++) {
            CellData cell = cell_array_current_[cell_num_current_ + j];
            if (cell.x < min_x || cell.x > max_x || cell.y < min_y || cell.y > max_y)
                continue;
            memmove(&(cell_array_current_[cell_num_current_]), &cell, sizeof(CellData));
            cell_id_array_current_[cell_num_current_] = offset + j;
            cell_num_current_++;
            expression_num_current_ += cell.gene_count;

            selectCellExp(cell.offset, cell.gene_count, cell_exp_data);
            for (unsigned int m = 0; m < cell.gene_count; m++) {
                gene_id_set_current_.insert(cell_exp_data[m].gene_id);
            }
        }
    }

    gene_num_current_ = gene_id_set_current_.size();
    free(cell_exp_data);
    if (verbose_) printCpuTime(cprev, "restrictRegion");
}

void CgefReader::restrictGene(vector<string> &gene_list, bool exclude) {
    restrict_gene_ = true;
    if (restrict_region_) {
        for (const auto &gene_name: gene_list) {
            if (genename_to_id_.find(gene_name) != genename_to_id_.end()) {
                if (exclude) {
                    auto iter = gene_id_set_current_.find(genename_to_id_[gene_name]);
                    if (iter != gene_id_set_current_.end())
                        gene_id_set_current_.erase(iter);
                } else {
                    gene_id_set_current_.insert(genename_to_id_[gene_name]);
                }
            }
        }
    } else {
        unordered_set<string> genename_set;
        for (const auto &gene_name: gene_list) {
            genename_set.insert(gene_name);
        }

//        if(gene_id_array_current_ == nullptr)
//            gene_id_array_current_ = (unsigned short*)malloc(gene_num_ * sizeof(unsigned short));
//
//
//        gene_num_current_ = 0;
        for (unsigned short i = 0; i < gene_num_; i++) {
            if ((exclude && genename_set.find(gene_array_[i].gene_name) == genename_set.end())
                || (!exclude && genename_set.find(gene_array_[i].gene_name) != genename_set.end())
                    ) {
                gene_id_set_current_.insert(i);
//                gene_id_array_current_[gene_num_current_] = i;
                //            memmove(&gene_array_[gene_num_current_], &gene_array_[i], sizeof(GeneData));
//                gene_num_current_++;
            }
        }
    }
    gene_num_current_ = gene_id_set_current_.size();
}

void CgefReader::updateGeneInfo() {
    auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
    for (unsigned int i = 0; i < cell_num_current_; i++) {
        CellData cell = cell_array_current_[i];
        selectCellExp(cell.offset, cell.gene_count, cell_exp_data);
        for (unsigned int j = 0; j < cell.gene_count; j++) {
            gene_id_set_current_.insert(cell_exp_data[j].gene_id);
        }
    }
    gene_num_current_ = gene_id_set_current_.size();
    free(cell_exp_data);
}

unsigned int CgefReader::getCellCount(string &gene_name) {
    auto iter = genename_to_id_.find(gene_name);
    if (iter == genename_to_id_.end())
        return 0;
    return gene_array_[iter->second].cell_count;
}

unsigned int CgefReader::getCellCount(unsigned short gene_id) {
    if (gene_id >= gene_num_)
        return 0;
    return gene_array_[gene_id].cell_count;
}

unsigned short CgefReader::getGeneCount(unsigned int cell_id) const {
    if (cell_id >= cell_num_)
        return 0;
    return getCell(cell_id).gene_count;
}

