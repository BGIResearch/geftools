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

void CellBin::storeVersion() {
    int version = 1;
    hsize_t dimsAttr[1] = {1};
    hid_t data_space = H5Screate_simple(1, dimsAttr, NULL);
    hid_t attr = H5Acreate(file_id_, "version", H5T_STD_U32LE, data_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &version);
    H5Sclose(data_space);
    H5Aclose(attr);
}

