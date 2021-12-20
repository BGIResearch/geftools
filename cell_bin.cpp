//
// Created by huangzhibo on 2021/12/15.
//

#include "cell_bin.h"

CellBin::CellBin(const string &filepath, const string &mode) {
    //    bool r = copyFile(inPath, outPath);
    //    if(!r) cerr << "Error writing to file (" << outPath << ") is failed!" << endl;

    printf("create h5 file: %s\n", filepath.c_str());
    file_id_ = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group_id_ = H5Gcreate(file_id_, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

CellBin::~CellBin() {
    H5Fclose(file_id_);
    H5Gclose(group_id_);
}

unsigned int CellBin::getGeneNum() const {
    return gene_num_;
}

unsigned int CellBin::getCellNum() const {
    return cell_num_;
}

unsigned long long int CellBin::getExpressionNum() const {
    return expression_num_;
}

void CellBin::storeGeneList(vector<string> & geneList) const {
    hsize_t dims[1] = {geneList.size()};
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "geneList", strtype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &geneList[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
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

void CellBin::storeCellExp() {

    hsize_t dims[1] = {cell_exp_list_.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, geneID), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);

    filetype = H5Tcreate(H5T_COMPOUND, 4);
    H5Tinsert(filetype, "geneID", 0, H5T_STD_U16LE);
    H5Tinsert(filetype, "count", 2, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(group_id_, "cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_exp_list_[0]);

    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CellBin::storeCell() {
    hsize_t dims[1] = {(hsize_t)cell_num_};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(memtype, "x", HOFFSET(CellData, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(CellData, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "offset", HOFFSET(CellData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "geneCount", HOFFSET(CellData, gene_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "expCount", HOFFSET(CellData, exp_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "dnbCount", HOFFSET(CellData, dnb_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "area", HOFFSET(CellData, area), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "cellTypeID", HOFFSET(CellData, cell_type_id), H5T_NATIVE_USHORT);

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
    cell_attr_.average_gene_count = static_cast<float>(gene_count_sum_)/  static_cast<float>(cell_num_);
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

void CellBin::addDnbInCell(vector<Point> & dnb_coordinates,
                           map<unsigned long long int, vector<unsigned int>> & bin_gene_exp_map,
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
            vector<unsigned int> gene_exps = iter->second;
            for (auto gene_exp : gene_exps) {
//                CellExpData cell_gene_exp = CellExpData(gene_exp);
                CellExpData cell_gene_exp{};
                memcpy(&cell_gene_exp, &gene_exp, sizeof(unsigned int));

                exp_count += cell_gene_exp.count;
                auto iter_gene = gene_count_in_cell.find(cell_gene_exp.geneID);
                if(iter_gene != gene_count_in_cell.end()){
                    iter_gene->second += cell_gene_exp.count;
                } else{
                    gene_count_in_cell.insert(
                            map<unsigned short, unsigned short>::value_type(
                                    cell_gene_exp.geneID, cell_gene_exp.count));
                    gene_count++;
                }
            }
        }
    }

    CellData cell = {
        static_cast<unsigned int>(center_point.x),
        static_cast<unsigned int>(center_point.y),
        current_cell_offset_,
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

    gene_count_sum_ += gene_count;
    exp_count_sum_ += exp_count;
    dnb_count_sum_ += dnb_count;
    area_sum_ += area;

    cell_list_.emplace_back(cell);
    cell_num_ += 1;

    current_cell_offset_ += gene_count;

    map<unsigned short, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
        CellExpData cexp_tmp = {iter_m->first, iter_m->second};
        cell_exp_list_.emplace_back(cexp_tmp);
        ++iter_m;
    }
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

}

