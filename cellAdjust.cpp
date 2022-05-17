/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:23
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-17 14:46:57
 * @Description: file content
 */
#include "cellAdjust.h"


cellAdjust::cellAdjust(string &bgef, string &cgef):
    m_bgef(bgef), m_cgef(cgef)
{
}

cellAdjust::~cellAdjust()
{
}


void cellAdjust::readBgef(const string &strinput)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(file_id, "/geneExp/bin1/gene", H5P_DEFAULT);
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);


    m_genePtr = (Gene*)malloc(dims[0] * sizeof(Gene));
    hid_t genememtype, strtype;
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    genememtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(genememtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(genememtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(genememtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    H5Dread(gene_did, genememtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_genePtr);
    H5Tclose(genememtype);
    H5Tclose(strtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);


    hid_t exp_did = H5Dopen(file_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_cgefwPtr->expression_num_ = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    m_expPtr = (Expression *) malloc(dims[0] * sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_expPtr);

    for(int i=0;i<)

    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    // attr = H5Aopen(exp_did, "maxExp", H5P_DEFAULT);
    // H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_exp));
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);
    H5Fclose(file_id);
}

void cellAdjust::readCgef(const string &strinput)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);

    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);

    hid_t memtype = getMemtypeOfCellData();
    CellData *cell_array = new CellData[dims[0]];
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array);

    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);

    hsize_t dims[3];
    hid_t border_id = H5Dopen(file_id, "/cellBin/cellBorder", H5P_DEFAULT);
    hid_t border_sid = H5Dget_space(border_id);
    H5Sget_simple_extent_dims(border_sid, dims, nullptr);

    short *borderdataPtr = (short*)calloc(dims[0]*dims[1]*dims[2], 2);
    H5Dread(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderdataPtr);

    for(int i=0;i<dims[0];i++)
    {
        for(int j=0;j<dims[1];j++)
        {
            if(borderdataPtr[j]==0)
            {
                break;
            }
            borderdataPtr[j] += cell_array[i].x;
            borderdataPtr[j] += cell_array[i].y;
        }
    }

    int min_x, min_y, max_x, max_y;
    hid_t attr = H5Aopen(border_id, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_x);
    attr = H5Aopen(border_id, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_y);
    attr = H5Aopen(border_id, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_x);
    attr = H5Aopen(border_id, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_y);

    H5Aclose(attr);
    H5Sclose(border_sid);
    H5Dclose(border_id);
    H5Fclose(file_id);
}

void cellAdjust(string &bgef, string &cgef)
{
    cellAdjust cellgene(bgef, cgef);
    cellgene.readBgef();
    cellgene.readCgef();
}