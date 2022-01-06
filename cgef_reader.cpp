
#include "cgef_reader.h"


CgefReader::CgefReader(const string &filename, bool verbose) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);
    verbose_ = verbose;

    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id_ = H5Gopen(file_id_, "/CellBin", H5P_DEFAULT);
    openCellDataset();
    openCellExpDataset();
    openGeneDataset();
    openGeneExpDataset();

    getGene();
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
    if(gene_array_ != nullptr)
        free(gene_array_);
    if(cell_array_ != nullptr)
        free(cell_array_);
}

void CgefReader::openCellDataset() {
    hsize_t dims[1];
    cell_dataset_id_ = H5Dopen(group_id_, "cell", H5P_DEFAULT);
    if (cell_dataset_id_ < 0){
        cerr<<"failed open dataset: cell" <<endl;
        return;
    }
    cell_dataspace_id_ = H5Dget_space(cell_dataset_id_);
    H5Sget_simple_extent_dims(cell_dataspace_id_, dims, nullptr);
    cell_num_ = dims[0];

    hsize_t dims_attr[2];
    hid_t attr, attr_dataspace;
    attr = H5Aopen (cell_dataset_id_, "blockIndex", H5P_DEFAULT);
    attr_dataspace = H5Aget_space(attr);
    H5Sget_simple_extent_dims (attr_dataspace, dims_attr, nullptr);

    block_num_ = dims_attr[0];
    block_index_ = static_cast<unsigned int *>(
            malloc(dims_attr[0] * dims_attr[1] * sizeof(unsigned int)));

    H5Aread (attr, H5T_NATIVE_UINT32, block_index_);

    attr = H5Aopen (cell_dataset_id_, "blockSize", H5P_DEFAULT);
    H5Aread (attr, H5T_NATIVE_UINT32, block_size_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
}

void CgefReader::openGeneDataset() {
    hsize_t dims[1];
    gene_dataset_id_ = H5Dopen(group_id_, "gene", H5P_DEFAULT);
    if (gene_dataset_id_ < 0){
        cerr<<"failed open dataset: gene" <<endl;
        return;
    }
    hid_t gene_dataspace_id = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id, dims, nullptr);
    gene_num_ = dims[0];
    H5Sclose(gene_dataspace_id);
}


void CgefReader::openCellExpDataset() {
    hsize_t dims[1];
    cell_exp_dataset_id_ = H5Dopen(group_id_, "cellExp", H5P_DEFAULT);
    if (cell_exp_dataset_id_ < 0){
        cerr<<"failed open dataset: cellExp" <<endl;
        return;
    }
    cell_exp_dataspace_id_ = H5Dget_space(cell_exp_dataset_id_);
    H5Sget_simple_extent_dims(cell_exp_dataspace_id_, dims, nullptr);
    expression_num_ = dims[0];
}

void CgefReader::openGeneExpDataset() {
    hsize_t dims[1];
    gene_exp_dataset_id_ = H5Dopen(group_id_, "geneExp", H5P_DEFAULT);
    if (gene_exp_dataset_id_ < 0){
        cerr<<"failed open dataset: geneExp" <<endl;
        return;
    }
    gene_exp_dataspace_id_ = H5Dget_space(gene_exp_dataset_id_);
    H5Sget_simple_extent_dims(gene_exp_dataspace_id_, dims, nullptr);
}

unsigned int CgefReader::getCellNum() const {
    return cell_num_;
}

unsigned short CgefReader::getGeneNum() const {
    return gene_num_;
}

unsigned long long int CgefReader::getExpressionNum() const {
    return expression_num_;
}

GeneData *CgefReader::getGene() {
    unsigned long cprev=clock();
    if(gene_array_ != nullptr)
        return gene_array_;

    hid_t memtype = getMemtypeOfGeneData();
    gene_array_ = (GeneData*)malloc(gene_num_ * sizeof(GeneData));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_array_);


    for(unsigned short i=0; i < gene_num_; i++){
        genename_to_id_[gene_array_[i].gene_name] = i;
    }

    H5Tclose(memtype);
    if(verbose_) printCpuTime(cprev, "getGene");
    return gene_array_;
}

CellData *CgefReader::getCell() {
    if(cell_array_ != nullptr)
        return cell_array_;

    hid_t memtype = getMemtypeOfCellData();
    cell_array_ = (CellData*)malloc(cell_num_ * sizeof(CellData));
    H5Dread(cell_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array_);
    H5Tclose(memtype);
    return cell_array_;
}

void CgefReader::getGeneNameList(char *gene_list) {
    GeneData * genes = getGene();
    for(unsigned int i = 0; i < gene_num_; i++){
        memcpy(&gene_list[i], genes[i].gene_name, 32);
    }
}

int CgefReader::getSparseMatrixIndicesOfExp(unsigned int *indices, unsigned int *indptr, unsigned int *count,
                                         const char *order) {
    if(order[0] == 'g'){
        getCellIdAndCount(indices, count);
        GeneData * gene_data = getGene();
        //indptr length = gene_num_ + 1
        indptr[0] = 0;
        for(unsigned int i = 1; i < cell_num_; i++){
            indptr[i] = gene_data->cell_count;
        }
    }else if(order[0] == 'c'){
        getGeneIdAndCount(indices, count);
        //indptr length = gene_num_ + 1
        CellData * cell_data = getCell();
        indptr[0] = 0;
        for(unsigned int i = 1; i < cell_num_; i++){
            indptr[i] = cell_data->gene_count;
        }
    }else {
        return -1;
    }
    return 0;
}

int CgefReader::getSparseMatrixIndicesOfExp2(unsigned int *cell_ind,
                                          unsigned int *gene_ind,
                                          unsigned int *count) {
    getCellIdAndCount(cell_ind, count);
    GeneData * gene_data = getGene();
    unsigned int n = 0;
    for(unsigned short i = 0; i < gene_num_; i++){
        unsigned int cell_count = gene_data->cell_count;
        for(unsigned int j = 0; j < cell_count; j++){
            gene_ind[n++] = i;
        }
    }
    return 0;
}

void CgefReader::getCellIdAndCount(unsigned int *cell_id, unsigned int *count) const{
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

void CgefReader::getGeneIdAndCount(unsigned int *gene_id, unsigned int *count) const{
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


    return 0;
}

unsigned int CgefReader::getExpressionCountByGeneId(unsigned short gene_id, GeneExpData *expressions) {
    GeneData * genes = getGene();
    expressions = static_cast<GeneExpData *>(malloc(genes[gene_id].cell_count * sizeof(GeneExpData)));

    getGeneExpByOffset(genes[gene_id].offset, genes[gene_id].cell_count, expressions);

    return genes[gene_id].cell_count;
}

void CgefReader::getGeneExpByOffset(unsigned int offset, unsigned int cell_count, GeneExpData *expressions) const {
    hsize_t start[1] = {offset},
            count[1] = {cell_count},
            offset_out[2] = {0, 0};

    hid_t memtype = getMemtypeOfGeneExpData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count,nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab (gene_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread (gene_exp_dataset_id_, memtype, memspace, gene_exp_dataspace_id_, H5P_DEFAULT, expressions);
}

GeneData CgefReader::getGeneDataByGeneId(unsigned short gene_id) {
    GeneData * genes = getGene();
    return genes[gene_id];
}

CellData CgefReader::getCellData(unsigned int cell_id) {
    return getCell()[cell_id];
}

string CgefReader::getGeneName(unsigned short gene_id) {
    GeneData gene_data = getGene()[gene_id];
    return gene_data.gene_name;
}

GeneData CgefReader::getGeneData(unsigned short gene_id) {
    return getGene()[gene_id];
}

int CgefReader::getGeneId(string &gene_name) {
    if (genename_to_id_.find(gene_name) != genename_to_id_.end()){
        return genename_to_id_[gene_name];
    }
    return -1;
}

unsigned int CgefReader::toGem(string &filename, const vector<string> &gene_name_list, bool force_genes, bool exclude) {
    unsigned long cprev=clock();
    unsigned short gene_ids[gene_num_];
    unsigned short n = 0;
    unsigned int cell_num = 0;
    unordered_map<string,bool> genename_map;
    for (const auto& gene_name: gene_name_list) {
        genename_map[gene_name] = true;
    }
    if(exclude){
        for(unsigned short i = 0; i < gene_num_; i++){
            if(genename_map.find(gene_array_[i].gene_name) == genename_map.end()){
                gene_ids[n++] = i;
            }
        }
    }else {
        for (const auto& gene_name: gene_name_list) {
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
    if(filename != "stdout"){
        fout.open(filename);
        if(!fout.is_open()) cerr << "Fail to open file : " << filename << endl;
        fout << "#geneName\tx\ty\tcount\tcellID" << endl;
    }else {
        to_file = false;
        cout << "#geneName\tx\ty\tcount\tcellID" << endl;
    }

    if(use_region_){
        if(verbose_) cerr << "toGem use_region_ true" << endl;
        auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
        for(unsigned int i = 0; i < cell_num_; i++){
            CellData cell = cell_array_[i];
            unsigned int cell_id = cell_id_array_[i];
            selectCellExp(cell.offset, cell.gene_count, cell_exp_data);

            for(unsigned int j = 0; j < cell.gene_count; j++){
                string gene_name = gene_array_[cell_exp_data[j].gene_id].gene_name;

                if ( !gene_name_list.empty() && (
                        (exclude && genename_map.find(gene_name) != genename_map.end())
                        || (!exclude && genename_map.find(gene_name) == genename_map.end())
                        ) ) continue;

                cell_num ++;
                if(to_file){
                    fout << gene_name << "\t" << cell.x << "\t" << cell.y
                         << "\t" << cell_exp_data[j].count << "\t" << cell_id <<endl;
                }else {
                    cout << gene_name << "\t" << cell.x << "\t" << cell.y
                         << "\t" << cell_exp_data[j].count << "\t" << cell_id <<endl;
                }
            }
        }

        free(cell_exp_data);
    }else {
        unsigned int max_malloc = 1000;
        auto *expression = static_cast<GeneExpData *>(malloc(max_malloc * sizeof(GeneExpData)));

        for (unsigned short i = 0; i < n; i++) {
            unsigned short gene_id = gene_ids[i];
            GeneData gene_data = getGeneDataByGeneId(gene_id);
            unsigned int cell_count = gene_data.cell_count;

            cell_num += cell_count;

            if(cell_count > max_malloc){
                expression = static_cast<GeneExpData *>(realloc(expression, cell_count));
                max_malloc = cell_count + 1000;
            }

            getGeneExpByOffset(gene_data.offset, cell_count, expression);

            for(unsigned int j = 0; j < cell_count; j++){
                CellData cell_data = getCellData(expression[j].cell_id);
                if(to_file){
                    fout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
                         << "\t" << expression[j].count << "\t" << expression[j].cell_id <<endl;
                }else {
                    cout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
                         << "\t" << expression[j].count << "\t" << expression[j].cell_id <<endl;
                }
            }
        }
        free(expression);
    }

    if(verbose_) printCpuTime(cprev, "toGem");

    if(to_file) fout.close();
    return cell_num;
}

bool CgefReader::isVerbose() const {
    return verbose_;
}

void CgefReader::setVerbose(bool verbose) {
    CgefReader::verbose_ = verbose;
}

void CgefReader::useRegion(unsigned int min_x, unsigned int max_x, unsigned int min_y, unsigned int max_y) {
    unsigned long cprev=clock();
    use_region_ = true;
    unsigned int x_block_num = block_size_[2];
    unsigned int min_block_x = min_x/block_size_[0];
    unsigned int max_block_x = max_x/block_size_[0];
    unsigned int min_block_y = min_y/block_size_[1];
    unsigned int max_block_y = max_y/block_size_[1];

    active_block_ids_num_ = (max_block_y - min_block_y + 1) * (max_block_x - min_block_x + 1);
    active_block_ids_ = static_cast<unsigned int *>(malloc(active_block_ids_num_ * sizeof(unsigned int)));

    auto b_arr = (int (*)[2])block_index_;

    unsigned int cell_num_tmp = 0;
    int index = 0;
    unsigned int block_id, block_id_y;
    for(unsigned int y = min_block_y; y <= max_block_y; y++) {
        block_id_y = y * x_block_num;
        for (unsigned int x = min_block_x; x <= max_block_x; x++) {
            block_id = x + block_id_y;
            active_block_ids_[index++] = block_id;
            cell_num_tmp += b_arr[block_id][1];
        }
    }

    cell_num_ = 0;
    auto * cell_data = static_cast<CellData *>(malloc(cell_num_tmp * sizeof(CellData)));
    cell_array_ = static_cast<CellData *>(malloc(cell_num_tmp * sizeof(CellData)));
    cell_id_array_ = static_cast<unsigned int *>(malloc(cell_num_tmp * sizeof(unsigned int)));

    unsigned int offset, count;
    for(unsigned int i = 0; i < index; i++){
        block_id = active_block_ids_[i];
        offset = b_arr[block_id][0];
        count = b_arr[block_id][1];

        selectCells(offset, count, cell_data);

        for(unsigned int j = 0; j < count; j++){
            CellData cell = cell_data[j];
            if(cell.x < min_x || cell.x > max_x || cell.y < min_y || cell.y > max_y)
                continue;
            memcpy(&(cell_array_[cell_num_]), &cell, sizeof(CellData));
            cell_id_array_[cell_num_] = offset + j ;
            cell_num_++;
        }
    }

    free(cell_data);
    if(verbose_) printCpuTime(cprev, "useRegion");
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
    hid_t memspace = H5Screate_simple(1, count,nullptr);
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
    hid_t memspace = H5Screate_simple(1, count,nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_exp_dataset_id_, memtype, memspace, cell_exp_dataspace_id_, H5P_DEFAULT, cell_exp_data);
}

unsigned int CgefReader::getActiveBlockId() {
    return 0;
}

bool CgefReader::isUseRegion() const {
    return use_region_;
}
