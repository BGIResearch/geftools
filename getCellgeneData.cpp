/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:23
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-16 11:15:21
 * @Description: file content
 */
#include "getCellgeneData.h"


getCellgeneData::getCellgeneData(string &bgef, string &cgef):
    m_bgef(bgef), m_cgef(cgef)
{
}

getCellgeneData::~getCellgeneData()
{
}


void getCellgeneData::readBgef(const string &strinput)
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

void getCellgeneData::readCgef(const string &strinput)
{

}

void getCellgeneData(string &bgef, string &cgef)
{
    getCellgeneData cellgene(bgef, cgef);
    cellgene.readBgef();
    cellgene.readCgef();
}