//
// Created by huangzhibo on 2021/12/15.
//

#include "cell_bin.h"

CellBin::CellBin(const string &filename, const string &mode) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);

    mode_ = mode;

    if(mode == "r"){
        file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        group_id_ = H5Gopen(file_id_, "/cellBin", H5P_DEFAULT);
        openCellDataset();
        openCellExpDataset();
        openGeneDataset();
        openGeneExpDataset();
    } else if(mode == "w"){
        cerr << "create h5 file: " <<  filename << endl;
        file_id_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        group_id_ = H5Gcreate(file_id_, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
}

CellBin::~CellBin() {
    H5Tclose(str32_type_);
    H5Fclose(file_id_);
    H5Gclose(group_id_);
    if(mode_ == "r"){
        H5Dclose(cell_dataset_id_);
        H5Dclose(gene_dataset_id_);
        H5Dclose(cell_exp_dataset_id_);
        H5Dclose(gene_exp_dataset_id_);
    }
    if(gene_array_ != nullptr)
        free(gene_array_);
    if(cell_array_ != nullptr)
        free(cell_array_);
}

void CellBin::openCellDataset() {
    hsize_t dims[1];
    cell_dataset_id_ = H5Dopen(group_id_, "cell", H5P_DEFAULT);
    if (cell_dataset_id_ < 0){
        cerr<<"failed open dataset: cell" <<endl;
        return;
    }
    hid_t cell_dataspace_id = H5Dget_space(cell_dataset_id_);
    H5Sget_simple_extent_dims(cell_dataspace_id, dims, NULL);
    cell_num_ = dims[0];
    H5Sclose(cell_dataspace_id);
}

void CellBin::openGeneDataset() {
    hsize_t dims[1];
    gene_dataset_id_ = H5Dopen(group_id_, "gene", H5P_DEFAULT);
    if (gene_dataset_id_ < 0){
        cerr<<"failed open dataset: gene" <<endl;
        return;
    }
    hid_t gene_dataspace_id = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id, dims, NULL);
    gene_num_ = dims[0];
    H5Sclose(gene_dataspace_id);
}


void CellBin::openCellExpDataset() {
    hsize_t dims[1];
    cell_exp_dataset_id_ = H5Dopen(group_id_, "cellExp", H5P_DEFAULT);
    if (cell_exp_dataset_id_ < 0){
        cerr<<"failed open dataset: cellExp" <<endl;
        return;
    }
    hid_t cell_exp_dataspace_id = H5Dget_space(cell_exp_dataset_id_);
    H5Sget_simple_extent_dims(cell_exp_dataspace_id, dims, NULL);
    expression_num_ = dims[0];
    H5Sclose(cell_exp_dataspace_id);
}

void CellBin::openGeneExpDataset() {
    hsize_t dims[1];
    gene_exp_dataset_id_ = H5Dopen(group_id_, "geneExp", H5P_DEFAULT);
    if (gene_exp_dataset_id_ < 0){
        cerr<<"failed open dataset: geneExp" <<endl;
        return;
    }
}

unsigned int CellBin::getCellNum() const {
    return cell_num_;
}

unsigned short CellBin::getGeneNum() const {
    return gene_num_;
}

unsigned long long int CellBin::getExpressionNum() const {
    return expression_num_;
}

GeneData *CellBin::getGene() {
    if(gene_array_ != nullptr)
        return gene_array_;

    hid_t memtype = getMemtypeOfGeneData();
    gene_array_ = (GeneData*)malloc(gene_num_ * sizeof(GeneData));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_array_);
    H5Tclose(memtype);
    return gene_array_;
}

CellData *CellBin::getCell() {
    if(cell_array_ != nullptr)
        return cell_array_;

    hid_t memtype = getMemtypeOfCellData();
    cell_array_ = (CellData*)malloc(cell_num_ * sizeof(CellData));
    H5Dread(cell_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array_);
    H5Tclose(memtype);
    return cell_array_;
}

//GeneExpData * CellBin::getGeneExp(){
//
//}

hid_t CellBin::getMemtypeOfGeneData() const{
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneData));
    H5Tinsert(memtype, "geneName", HOFFSET(GeneData, gene_name), str32_type_);
    H5Tinsert(memtype, "offset", HOFFSET(GeneData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "cellCount", HOFFSET(GeneData, cell_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "maxMIDcount", HOFFSET(GeneData, max_mid_count), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t CellBin::getMemtypeOfGeneExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneExpData));
    H5Tinsert(memtype, "cellID", HOFFSET(GeneExpData, cell_id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(GeneExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}


hid_t CellBin::getMemtypeOfCellData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(memtype, "x", HOFFSET(CellData, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(CellData, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "offset", HOFFSET(CellData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "geneCount", HOFFSET(CellData, gene_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "expCount", HOFFSET(CellData, exp_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "dnbCount", HOFFSET(CellData, dnb_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "area", HOFFSET(CellData, area), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "cellTypeID", HOFFSET(CellData, cell_type_id), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t CellBin::getMemtypeOfCellExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}

void CellBin::storeCellBorder(char* borderPath, unsigned int cell_num) const {
    hsize_t dims[3];
    dims[0] = cell_num;
    dims[1] = 16;
    dims[2] = 2;

    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "cellBorder", H5T_STD_I8LE, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderPath);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}


void CellBin::storeCellBorderWithAttr(char *borderPath, unsigned int cell_num, unsigned int *effective_rect) const {
    storeCellBorder(borderPath, cell_num);

    hid_t dataset_id = H5Dopen(group_id_, "cellBorder", H5P_DEFAULT);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, NULL);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, effective_rect);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &effective_rect[1]);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &effective_rect[2]);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &effective_rect[3]);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Dclose(dataset_id);
}

void CellBin::storeCellExp() {

    hsize_t dims[1] = {cell_exp_list_.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);

    filetype = H5Tcreate(H5T_COMPOUND, 4);
    H5Tinsert(filetype, "geneID", 0, H5T_STD_U16LE);
    H5Tinsert(filetype, "count", 2, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_exp_list_[0]);


    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, NULL);
    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CellBin::storeCell() {
    hsize_t dims[1] = {(hsize_t)cell_num_};

    hid_t memtype, filetype;
    memtype = getMemtypeOfCellData();

    filetype = H5Tcreate(H5T_COMPOUND, 22);
    H5Tinsert(filetype, "x", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "y", 4, H5T_STD_U32LE);
    H5Tinsert(filetype, "offset", 8, H5T_STD_U32LE);
    H5Tinsert(filetype, "geneCount", 12, H5T_STD_U16LE);
    H5Tinsert(filetype, "expCount", 14, H5T_STD_U16LE);
    H5Tinsert(filetype, "dnbCount", 16, H5T_STD_U16LE);
    H5Tinsert(filetype, "area", 18, H5T_STD_U16LE);
    H5Tinsert(filetype, "cellTypeID", 20, H5T_STD_U16LE);
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

    hid_t dataset_id = H5Dcreate(group_id_, "cell", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_list_[0]);

    // Create cell attribute
    cell_attr_.average_gene_count = static_cast<float>(expression_num_)/  static_cast<float>(cell_num_);
    cell_attr_.average_exp_count = static_cast<float>(exp_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_dnb_count = static_cast<float>(dnb_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_area = static_cast<float>(area_sum_) /  static_cast<float>(cell_num_);

    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, NULL);
    attr = H5Acreate(dataset_id, "averageGeneCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_gene_count);
    attr = H5Acreate(dataset_id, "averageExpCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_exp_count);
    attr = H5Acreate(dataset_id, "averageDnbCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_dnb_count);
    attr = H5Acreate(dataset_id, "averageArea", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_area);
    attr = H5Acreate(dataset_id, "minGeneCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_gene_count);
    attr = H5Acreate(dataset_id, "minExpCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_exp_count);
    attr = H5Acreate(dataset_id, "minDnbCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_dnb_count);
    attr = H5Acreate(dataset_id, "minArea", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_area);
    attr = H5Acreate(dataset_id, "maxGeneCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_gene_count);
    attr = H5Acreate(dataset_id, "maxExpCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_exp_count);
    attr = H5Acreate(dataset_id, "maxDnbCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_dnb_count);
    attr = H5Acreate(dataset_id, "maxArea", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_area);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CellBin::addDnbExp(vector<Point> & dnb_coordinates,
                           map<unsigned long long int, vector<CellExpData>> & bin_gene_exp_map,
                           const Point& center_point,
                           unsigned short area) {
    unsigned long long int bin_id;
    map<unsigned short, unsigned short> gene_count_in_cell;
    unsigned short gene_count = 0;
    unsigned short exp_count = 0;

    for (auto & dnb_coordinate : dnb_coordinates) {
        bin_id = static_cast<unsigned long long int>(dnb_coordinate.x);
        bin_id = bin_id << 32 | static_cast<unsigned int>(dnb_coordinate.y);
        auto iter = bin_gene_exp_map.find(bin_id);
        if(iter != bin_gene_exp_map.end()){
            vector<CellExpData> gene_exps = iter->second;
            for (auto cell_gene_exp : gene_exps) {
                exp_count += cell_gene_exp.count;
                auto iter_gene = gene_count_in_cell.find(cell_gene_exp.gene_id);
                if(iter_gene != gene_count_in_cell.end()){
                    iter_gene->second += cell_gene_exp.count;
                } else{
                    gene_count_in_cell.insert(
                            map<unsigned short, unsigned short>::value_type(
                                    cell_gene_exp.gene_id, cell_gene_exp.count));
                    gene_count++;
                }
            }
        }
    }

    CellData cell = {
        static_cast<unsigned int>(center_point.x),
        static_cast<unsigned int>(center_point.y),
        expression_num_,
        static_cast<unsigned short>(gene_count),
        static_cast<unsigned short>(exp_count),
        static_cast<unsigned short>(dnb_coordinates.size()),
        area,
        0
    };

    cell_attr_.min_area = area < cell_attr_.min_area ? area : cell_attr_.min_area;
    cell_attr_.max_area = area > cell_attr_.max_area ? area : cell_attr_.max_area;
    cell_attr_.min_gene_count = gene_count < cell_attr_.min_gene_count ? gene_count : cell_attr_.min_gene_count;
    cell_attr_.max_gene_count = gene_count > cell_attr_.max_gene_count ? gene_count : cell_attr_.max_gene_count;
    cell_attr_.min_exp_count = exp_count < cell_attr_.min_exp_count ? exp_count : cell_attr_.min_exp_count;
    cell_attr_.max_exp_count = exp_count > cell_attr_.max_exp_count ? exp_count : cell_attr_.max_exp_count;

    unsigned int dnb_count = dnb_coordinates.size();
    cell_attr_.min_dnb_count = dnb_count < cell_attr_.min_dnb_count ? dnb_count : cell_attr_.min_dnb_count;
    cell_attr_.max_dnb_count = dnb_count > cell_attr_.max_dnb_count ? dnb_count : cell_attr_.max_dnb_count;

    exp_count_sum_ += exp_count;
    dnb_count_sum_ += dnb_count;
    area_sum_ += area;

    cell_list_.emplace_back(cell);

    expression_num_ += gene_count;

    map<unsigned short, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
//        gene_array_[iter_m->first].
        unsigned short gene_id = iter_m->first;
        unsigned short count = iter_m->second;
        max_mid_count_ = count > max_mid_count_ ? count : max_mid_count_;

        // 用于生成geneExp dataset的数据
        GeneExpData gene_exp_tmp = {cell_num_, count};
        auto iter_gene_exp_map = gene_exp_map_.find(gene_id);
        if(iter_gene_exp_map != gene_exp_map_.end()){
            iter_gene_exp_map->second.emplace_back(gene_exp_tmp);
        } else{
            vector<GeneExpData> v_gene_exp_data;
            v_gene_exp_data.emplace_back(gene_exp_tmp);
            gene_exp_map_.insert(map<unsigned short, vector<GeneExpData>>::value_type(gene_id, v_gene_exp_data));
        }

        CellExpData cexp_tmp = {gene_id, count};
        cell_exp_list_.emplace_back(cexp_tmp);
        ++iter_m;
    }

    cell_num_ += 1;
}

void CellBin::storeAttr(CellBinAttr & cell_bin_attr) const {
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, NULL);
    attr = H5Acreate(file_id_, "version", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.version);
    attr = H5Acreate(file_id_, "resolution", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.resolution);
    attr = H5Acreate(file_id_, "offsetX", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.offsetX);
    attr = H5Acreate(file_id_, "offsetY", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.offsetY);

    //Write createTime into cell bin gef
//    S32 time_str = getStrfTime();
//    attr = H5Acreate(file_id_, "createTime", str32_type_, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
//    H5Awrite(attr, str32_type_, &time_str);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
}

void CellBin::storeGeneAndGeneExp(const vector<string> &gene_name_list) {
    gene_num_ = gene_name_list.size();
    hsize_t dims[1] = {gene_num_};

    GeneData* gene_data_list;
    gene_data_list = static_cast<GeneData *>(malloc(gene_num_ * sizeof(GeneData)));

    unsigned int offset = 0;
    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(expression_num_);
    for(unsigned int i = 0; i < gene_num_; i++){
        auto iter_gene_exp_map = gene_exp_map_.find(i);
        if(iter_gene_exp_map != gene_exp_map_.end()){
            vector<GeneExpData> tmp = iter_gene_exp_map->second;
            gene_exp_list.insert(gene_exp_list.end(), tmp.begin(), tmp.end());

            gene_data_list[i] = {
                    gene_name_list[i].c_str(),
                    offset,
                    static_cast<unsigned int>(tmp.size()),
                    calcMaxCountOfGeneExp(tmp)};

            offset += tmp.size();
        }
    }

    hid_t memtype, filetype;
    memtype = getMemtypeOfGeneData();
    filetype = H5Tcreate(H5T_COMPOUND, 42);
    H5Tinsert(filetype, "geneName", 0, str32_type_);
    H5Tinsert(filetype, "offset", 32, H5T_STD_U32LE);
    H5Tinsert(filetype, "cellCount", 36, H5T_STD_U32LE);
    H5Tinsert(filetype, "maxMIDcount", 40, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "gene", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_data_list);

    memtype = getMemtypeOfGeneExpData();
    filetype = H5Tcreate(H5T_COMPOUND, 6);
    H5Tinsert(filetype, "cellID", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 4, H5T_STD_U16LE);

    hsize_t dims_exp[1] = {expression_num_};
    dataspace_id = H5Screate_simple(1, dims_exp, NULL);
    dataset_id = H5Dcreate(group_id_, "geneExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gene_exp_list[0]);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, NULL);
    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CellBin::storeCellTypeList() {
    hsize_t dims[1] = {1};

    S32 cell_type = S32("default");
    cell_type_list_.emplace_back(cell_type);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "cellTypeList", str32_type_, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, str32_type_, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_type_list_[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

unsigned short CellBin::calcMaxCountOfGeneExp(vector<GeneExpData> &gene_exps) {
    unsigned max = 0;
    for (auto gene_exp :gene_exps) {
        max = gene_exp.count > max ? gene_exp.count : max;
    }
    return max;
}

void CellBin::getGeneNameList(char *gene_list) {
    GeneData * genes = getGene();
    for(unsigned int i = 0; i < gene_num_; i++){
        memcpy(&gene_list[i], genes[i].gene_name, 32);
    }
}

int CellBin::getSparseMatrixIndicesOfExp(unsigned int *indices, unsigned int *indptr, unsigned int *count,
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

int CellBin::getSparseMatrixIndicesOfExp2(unsigned int *cell_ind,
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

void CellBin::getCellIdAndCount(unsigned int *cell_id, unsigned int *count) const{
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

void CellBin::getGeneIdAndCount(unsigned int *gene_id, unsigned int *count) const{
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
