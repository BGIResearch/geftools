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


//void CellBin::addDnbInCell(unsigned int *dnb_coordinates, unsigned int size) {
void CellBin::addDnbInCell(vector<Point> & dnb_coordinates, Point center_point, unsigned short area) {
    unsigned long long int bin_id;
    map<unsigned short, unsigned short> gene_count_in_cell;
    unsigned short gene_count = 0;
    unsigned short exp_count = 0;



    for (unsigned int i = 0; i < dnb_coordinates.size(); ++i) {
        bin_id = static_cast<unsigned long long int>(dnb_coordinates[i].x);
        bin_id = bin_id << 32 | static_cast<unsigned int>(dnb_coordinates[i].y);
        auto iter = gene_exp_map_.find(bin_id);
        if(iter != gene_exp_map_.end()){
            vector<CellExpData> cxp = iter->second;
            auto it = cxp.begin();
            while (it != cxp.end()) {
                exp_count += it->count;
                auto iter_gene = gene_count_in_cell.find(it->geneID);
                if(iter_gene != gene_count_in_cell.end()){
                    iter_gene->second += it->count;
                } else{
                    gene_count_in_cell.insert(map<unsigned short, unsigned short>::value_type(it->geneID, it->count));
                    gene_count++;
                }
                ++it;
            }
        }
    }

    CellData cell = {
        static_cast<unsigned int>(center_point.x),
        static_cast<unsigned int>(center_point.y),
        0,
        static_cast<unsigned short>(gene_count),
        static_cast<unsigned short>(exp_count),
        static_cast<unsigned short>(dnb_coordinates.size()),
        area,
        0
    };

    map<unsigned short, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
        CellExpData cexp_tmp = {iter_m->first, iter_m->second};
        cell_exp_list.emplace_back(cexp_tmp);
        ++iter_m;
    }



    cell_exp_count_list.emplace_back(exp_count);

    if(cell_gene_exp_list.empty()){
        cell_gene_exp_list.emplace_back(0);
    }else {
        cell_gene_exp_list.emplace_back(cell_gene_exp_list.back()+cell_gene_count_list.back());
    }
    cell_gene_count_list.emplace_back(gene_count);
}

