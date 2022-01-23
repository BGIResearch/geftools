#include "cgef_writer.h"

CgefWriter::CgefWriter(const string& output_cell_gef, bool verbose) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);
    verbose_ = verbose;

    cerr << "create h5 file: " <<  output_cell_gef << endl;
    hid_t fpid = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fpid, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    file_id_ = H5Fcreate(output_cell_gef.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fpid);
    group_id_ = H5Gcreate(file_id_, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(fpid);
}

CgefWriter::~CgefWriter() {
    H5Tclose(str32_type_);
    H5Fclose(file_id_);
    H5Gclose(group_id_);
}

void CgefWriter::storeCellBorder(char* borderPath, unsigned int cell_num) const {
    hsize_t dims[3];
    dims[0] = cell_num;
    dims[1] = 16;
    dims[2] = 2;

    hid_t dataspace_id = H5Screate_simple(3, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellBorder", H5T_STD_I8LE, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderPath);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}


void CgefWriter::storeCellBorderWithAttr(char *borderPath, unsigned int cell_num, unsigned int *effective_rect) const {
    storeCellBorder(borderPath, cell_num);

    hid_t dataset_id = H5Dopen(group_id_, "cellBorder", H5P_DEFAULT);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
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

void CgefWriter::storeCellExp() {

    hsize_t dims[1] = {cell_exp_list_.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);

    filetype = H5Tcreate(H5T_COMPOUND, 4);
    H5Tinsert(filetype, "geneID", 0, H5T_STD_U16LE);
    H5Tinsert(filetype, "count", 2, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_exp_list_[0]);


    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CgefWriter::storeCell(unsigned int block_num, unsigned int * block_index, const unsigned int *block_size) {
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
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);


    hid_t dpid = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_attr_phase_change(dpid, 0, 0);
    hid_t dataset_id = H5Dcreate(group_id_, "cell", filetype, dataspace_id, H5P_DEFAULT,
                                 dpid, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_list_[0]);

    // Create cell attribute
    cell_attr_.average_gene_count = static_cast<float>(expression_num_)/  static_cast<float>(cell_num_);
    cell_attr_.average_exp_count = static_cast<float>(exp_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_dnb_count = static_cast<float>(dnb_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_area = static_cast<float>(area_sum_) /  static_cast<float>(cell_num_);

    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "averageGeneCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_gene_count);
    attr = H5Acreate(dataset_id, "averageExpCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_exp_count);
    attr = H5Acreate(dataset_id, "averageDnbCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_dnb_count);
    attr = H5Acreate(dataset_id, "averageArea", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_area);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_attr_.min_x);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_attr_.max_x);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_attr_.min_y);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_attr_.max_y);
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

    // write block index
    dimsAttr[0] = block_num + 1;
    attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "blockIndex", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, block_index);

    dimsAttr[0] = 4;
    attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "blockSize", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, block_size);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(attr_dataspace);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}


void CgefWriter::addDnbExp(vector<Point> & dnb_coordinates,
                           map<unsigned long long int, pair<unsigned int, unsigned short>> & bin_gene_exp_map,
                           const DnbExpression *dnb_expression,
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
            pair<unsigned int, unsigned short> gene_info = iter->second;
            unsigned int end = gene_info.first + gene_info.second;
            for(unsigned int i = gene_info.first; i < end; i++){
                exp_count += dnb_expression[i].count;
                auto iter_gene = gene_count_in_cell.find(dnb_expression[i].gene_id);
                if(iter_gene != gene_count_in_cell.end()){
                    iter_gene->second += dnb_expression[i].count;
                } else{
                    gene_count_in_cell.insert(
                            map<unsigned short, unsigned short>::value_type(
                                    dnb_expression[i].gene_id, dnb_expression[i].count));
                    gene_count++;
                }
            }
        }
    }

    unsigned short cell_type_id = random_cell_type_num_ == 0 ? 0 : rand()%(random_cell_type_num_ + 1);

    CellData cell = {
            static_cast<unsigned int>(center_point.x),
            static_cast<unsigned int>(center_point.y),
            expression_num_, //offset
            static_cast<unsigned short>(gene_count),
            static_cast<unsigned short>(exp_count),
            static_cast<unsigned short>(dnb_coordinates.size()),
            area,
            cell_type_id
    };
    expression_num_ += gene_count;

    cell_attr_.min_x = cell.x < cell_attr_.min_x ? cell.x : cell_attr_.min_x;
    cell_attr_.max_x = cell.x > cell_attr_.max_x ? cell.x : cell_attr_.max_x;
    cell_attr_.min_y = cell.y < cell_attr_.min_y ? cell.y : cell_attr_.min_y;
    cell_attr_.max_y = cell.y > cell_attr_.max_y ? cell.y : cell_attr_.max_y;

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

    map<unsigned short, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
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


void CgefWriter::addDnbExp(vector<Point> & dnb_coordinates,
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

    unsigned short cell_type_id = random_cell_type_num_ == 0 ? 0 : rand()%(random_cell_type_num_ + 1);

    CellData cell = {
            static_cast<unsigned int>(center_point.x),
            static_cast<unsigned int>(center_point.y),
            expression_num_, //offset
            static_cast<unsigned short>(gene_count),
            static_cast<unsigned short>(exp_count),
            static_cast<unsigned short>(dnb_coordinates.size()),
            area,
            cell_type_id
    };
    expression_num_ += gene_count;

    cell_attr_.min_x = cell.x < cell_attr_.min_x ? cell.x : cell_attr_.min_x;
    cell_attr_.max_x = cell.x > cell_attr_.max_x ? cell.x : cell_attr_.max_x;
    cell_attr_.min_y = cell.y < cell_attr_.min_y ? cell.y : cell_attr_.min_y;
    cell_attr_.max_y = cell.y > cell_attr_.max_y ? cell.y : cell_attr_.max_y;

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

    map<unsigned short, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
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

void CgefWriter::storeAttr(CellBinAttr & cell_bin_attr) const {
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
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

void CgefWriter::storeGeneAndGeneExp(const vector<string> &gene_name_list) {
    gene_num_ = gene_name_list.size();
    hsize_t dims[1] = {gene_num_};

    GeneData* gene_data_list;
    gene_data_list = static_cast<GeneData *>(malloc(gene_num_ * sizeof(GeneData)));

    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count;

    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(expression_num_);
    for(unsigned int i = 0; i < gene_num_; i++){
        auto iter_gene_exp_map = gene_exp_map_.find(i);
        if(iter_gene_exp_map != gene_exp_map_.end()){
            vector<GeneExpData> tmp = iter_gene_exp_map->second;
            gene_exp_list.insert(gene_exp_list.end(), tmp.begin(), tmp.end());

            cell_count = static_cast<unsigned int>(tmp.size());
            max_MID_count = 0;
            exp_count = 0;
            for (auto gene_exp :tmp) {
                exp_count += gene_exp.count;
                max_MID_count = gene_exp.count > max_MID_count ? gene_exp.count : max_MID_count;
            }

            min_exp_count = min_exp_count < exp_count ? min_exp_count : exp_count;
            max_exp_count = max_exp_count > exp_count ? max_exp_count : exp_count;
            min_cell_count = min_cell_count < cell_count ? min_cell_count : cell_count;
            max_cell_count = max_cell_count > cell_count ? max_cell_count : cell_count;

            gene_data_list[i] = GeneData(
                    gene_name_list[i].c_str(),
                    offset,
                    static_cast<unsigned int>(tmp.size()),
                    exp_count,
                    max_MID_count);

            offset += tmp.size();
        }else{
            gene_data_list[i] = GeneData(
                    gene_name_list[i].c_str(),
                    offset,
                    0,
                    0,
                    0);
        }
    }

    hid_t memtype, filetype;
    memtype = getMemtypeOfGeneData();
    filetype = H5Tcreate(H5T_COMPOUND, 46);
    H5Tinsert(filetype, "geneName", 0, str32_type_);
    H5Tinsert(filetype, "offset", 32, H5T_STD_U32LE);
    H5Tinsert(filetype, "cellCount", 36, H5T_STD_U32LE);
    H5Tinsert(filetype, "expCount", 40, H5T_STD_U32LE);
    H5Tinsert(filetype, "maxMIDcount", 44, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "gene", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_data_list);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    attr = H5Acreate(dataset_id, "minExpCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &min_exp_count);
    attr = H5Acreate(dataset_id, "maxExpCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &max_exp_count);
    attr = H5Acreate(dataset_id, "minCellCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &min_cell_count);
    attr = H5Acreate(dataset_id, "maxCellCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &max_cell_count);

    memtype = getMemtypeOfGeneExpData();
    filetype = H5Tcreate(H5T_COMPOUND, 6);
    H5Tinsert(filetype, "cellID", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 4, H5T_STD_U16LE);

    hsize_t dims_exp[1] = {expression_num_};
    dataspace_id = H5Screate_simple(1, dims_exp, nullptr);
    dataset_id = H5Dcreate(group_id_, "geneExp", filetype, dataspace_id, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gene_exp_list[0]);

    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CgefWriter::storeCellTypeList() {
    S32 cell_type = S32("default");
    cell_type_list_.emplace_back(cell_type);

    int i = 0;
    while(i < random_cell_type_num_){
        i++;
        cell_type = S32();
        sprintf(cell_type.value, "type%d", i);
        cell_type_list_.emplace_back(cell_type);
    }

    hsize_t dims[1];
    dims[0] = random_cell_type_num_ + 1;
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellTypeList", str32_type_, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, str32_type_, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_type_list_[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

unsigned short CgefWriter::calcMaxCountOfGeneExp(vector<GeneExpData> &gene_exps) {
    unsigned max = 0;
    for (auto gene_exp :gene_exps) {
        max = gene_exp.count > max ? gene_exp.count : max;
    }
    return max;
}

int CgefWriter::write(BgefReader &common_bin_gef, Mask &mask) {
    map<unsigned long long int, pair<unsigned int, unsigned short>> bin_gene_exp_map;
    auto * dnb_exp_info = (DnbExpression *) malloc(common_bin_gef.getExpressionNum()  * sizeof(DnbExpression));
    common_bin_gef.getBinGeneExpMap(bin_gene_exp_map, dnb_exp_info);

    const vector<Polygon>& polygons = mask.getPolygons();

    unsigned long cprev=clock();
    for(unsigned int i = 0; i < mask.getCellNum(); i++){
        Polygon p = polygons[i];
        Rect roi = Rect(p.getMinX(),p.getMinY(),p.getCols(),p.getRows());
        Mat roi_mat = common_bin_gef.getWholeExpMatrix(roi);
        Mat fill_points = p.getFillPolyMat();
        roi_mat = roi_mat.mul(fill_points);

        vector<Point> non_zero_coordinates, non_zero_coordinates_offset;
        findNonZero(roi_mat,non_zero_coordinates);
        Point offset = Point(-p.getMinX(), -p.getMinY());
        offsetCoordinates(non_zero_coordinates, non_zero_coordinates_offset, offset);

        addDnbExp(
            non_zero_coordinates_offset,
            bin_gene_exp_map,
            dnb_exp_info,
            p.getCenter(),
            p.getAreaUshort());
    }

    if(verbose_) printCpuTime(cprev, "addDnbExp");

    char* borders = static_cast<char *>(malloc(mask.getCellNum() * 16 * 2 * sizeof(char)));
    mask.getBorders(borders);

    ExpressionAttr expression_attr = common_bin_gef.getExpressionAttr();
    CellBinAttr cell_bin_attr = {
            .version = 1,
            .resolution = expression_attr.resolution,
            .offsetX = expression_attr.min_x,
            .offsetY = expression_attr.min_y
    };

    storeAttr(cell_bin_attr);

    unsigned int effective_rect[4];
    mask.getEffectiveRectangle(effective_rect);
    storeCellBorderWithAttr(borders, mask.getCellNum(), effective_rect);
    storeCell(mask.getBlockNum(), mask.getBlockIndex(), mask.getBlockSize());
    storeCellExp();
    storeCellTypeList();

    vector<string> gene_name_list;
    gene_name_list.reserve(common_bin_gef.getGeneNum());
    common_bin_gef.getGeneNameList(gene_name_list);
    storeGeneAndGeneExp(gene_name_list);

    free(dnb_exp_info);
    return 0;
}

unsigned short CgefWriter::getRandomCellTypeNum() const {
    return random_cell_type_num_;
}

void CgefWriter::setRandomCellTypeNum(unsigned short random_cell_type_num) {
    random_cell_type_num_ = random_cell_type_num;
}

bool CgefWriter::isVerbose() const {
    return verbose_;
}

void CgefWriter::setVerbose(bool verbose) {
    verbose_ = verbose;
}
