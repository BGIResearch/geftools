/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:23
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 11:35:58
 * @Description: file content
 */
#include "cellAdjust.h"
#include "timer.h"


cellAdjust::cellAdjust()
{
}

cellAdjust::~cellAdjust()
{
    delete []m_cell_arrayptr;
    free(m_blkidxPtr);
}

void cellAdjust::readBgef(const string &strinput)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(file_id, "/geneExp/bin1/gene", H5P_DEFAULT);
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    m_genencnt = dims[0];
    Gene *genePtr = (Gene*)malloc(dims[0] * sizeof(Gene));
    hid_t genememtype, strtype;
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    genememtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(genememtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(genememtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(genememtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    H5Dread(gene_did, genememtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genePtr);
    H5Tclose(genememtype);
    H5Tclose(strtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);


    hid_t exp_did = H5Dopen(file_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneexpcnt = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    Expression *expPtr = (Expression *) malloc(dims[0] * sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expPtr);

    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", m_min_x, m_min_y, m_max_x, m_max_y);
    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);
    H5Fclose(file_id);

    uint64_t l_id = 0;
    for(int i=0;i<m_genencnt;i++)
    {
        m_vecgenename.emplace_back(genePtr[i].gene);
        Expression *ptr = expPtr + genePtr[i].offset;
        for(int j=0;j<genePtr[i].count;j++)
        {
            l_id = ptr[j].x;
            l_id = (l_id<<32) | ptr[j].y;
            
            if(m_hash_vecdnb.find(l_id) == m_hash_vecdnb.end())
            {
                vector<Dnbs> tvec;
                m_hash_vecdnb.emplace(l_id, tvec);
            }
            m_hash_vecdnb[l_id].emplace_back(i, ptr[j].count);
        }
    }
    free(genePtr);
    free(expPtr);
    printf("gene:%d geneexp:%d hashcnt:%d\n", m_genencnt, m_geneexpcnt, m_hash_vecdnb.size());
}

void cellAdjust::readCgef(const string &strinput)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);

    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);

    m_cellcnt = dims[0];
    hid_t memtype = getMemtypeOfCellData();
    m_cell_arrayptr = new CellData[dims[0]];
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cell_arrayptr);

    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);

    hid_t d_id = H5Dopen(file_id, "/cellBin/blockIndex", H5P_DEFAULT);
    hid_t s_id = H5Dget_space(d_id);
    H5Sget_simple_extent_dims(s_id, dims, nullptr);
    m_blkidxPtr = static_cast<unsigned int *>(calloc(dims[0], sizeof(unsigned int)));
    H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_blkidxPtr);

    H5Sclose(s_id);
    H5Dclose(d_id);

    d_id = H5Dopen(file_id, "/cellBin/blockSize", H5P_DEFAULT);
    H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_block_size);
    H5Dclose(d_id);

    hsize_t bdims[3];
    hid_t border_id = H5Dopen(file_id, "/cellBin/cellBorder", H5P_DEFAULT);
    hid_t border_sid = H5Dget_space(border_id);
    H5Sget_simple_extent_dims(border_sid, bdims, nullptr);

    short *borderdataPtr = (short*)calloc(bdims[0]*bdims[1]*bdims[2], 2);
    H5Dread(border_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderdataPtr);

    int x,y;
    int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
    vector<Point> vecborder;
    //m_fill_points = Mat::zeros(m_max_y-m_min_y+1, m_max_x-m_min_x+1, CV_8UC1);
    m_fill_points = Mat::zeros(20501, 21801, CV_8UC1);
    short *ptmp = borderdataPtr;
    for(int i=0;i<bdims[0];i++)
    {
        minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
        vecborder.clear();
        for(int j=0;j<bdims[1];j++)
        {
            x = ptmp[2*j];
            y = ptmp[2*j+1];
            if(x==SHRT_MAX && y==SHRT_MAX)
            {
                break;
            }
            x += m_cell_arrayptr[i].x;
            y += m_cell_arrayptr[i].y;
            vecborder.emplace_back(x,y);
            minx = std::min(minx, x);
            maxx = std::max(maxx, x);
            miny = std::min(miny, y);
            maxy = std::max(maxy, y);

        }
        if(!vecborder.empty())
        {
            Rect rtmp(minx, miny, maxx-minx+1, maxy-miny+1);
            m_hash_cellrect.emplace(i, rtmp);
            fillPoly(m_fill_points, vecborder, 1);
        }
        else
        {
            printf("cid %d\n",i);
        }
        ptmp += 64;
    }

    free(borderdataPtr);
    printf("cellcnt:%d hashrect %d\n", m_cellcnt, m_hash_cellrect.size());


    int min_x, min_y, max_x, max_y;
    hid_t attr = H5Aopen(border_id, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_x);
    attr = H5Aopen(border_id, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_y);
    attr = H5Aopen(border_id, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_x);
    attr = H5Aopen(border_id, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_y);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", min_x, min_y, max_x, max_y);

    // attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    // H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);

    // attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    // H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    H5Aclose(attr);
    H5Sclose(border_sid);
    H5Dclose(border_id);
    H5Fclose(file_id);

}

uint32_t cellAdjust::getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem)
{
    genename.reserve(m_vecgenename.size());
    vecCellgem.reserve(m_geneexpcnt);

    uint64_t l_id = 0;
    vector<Point> vecpoint;
    int x,y;
    int t_f = 0, t_n = 0;
    auto itor = m_hash_cellrect.begin();
    for(;itor != m_hash_cellrect.end();itor++)
    {
        vecpoint.clear();
        Mat t = m_fill_points(itor->second);
        findNonZero(t,vecpoint);

        for(Point &pt : vecpoint)
        {
            x = pt.x+itor->second.x;
            y = pt.y+itor->second.y;
            l_id = x;
            l_id = (l_id << 32) | y;
            auto dnb_itor = m_hash_vecdnb.find(l_id);
            if(dnb_itor!= m_hash_vecdnb.end())
            {
                t_f++;
                for(Dnbs &dnbs : dnb_itor->second)
                {
                    vecCellgem.emplace_back(dnbs.geneid, x, y, dnbs.midcnt, itor->first+1);
                }
                m_hash_vecdnb.erase(l_id);
            }
        }
    }
    printf("%d\n", m_hash_vecdnb.size());

    auto itor_s = m_hash_vecdnb.begin();
    for(;itor_s != m_hash_vecdnb.end();itor_s++)
    {
        t_n++;
        x = (itor_s->first) >> 32;
        y = (itor_s->first) | 0xFFFFFFFF;
        for(Dnbs &dnbs : itor_s->second)
        {
            vecCellgem.emplace_back(dnbs.geneid, x, y, dnbs.midcnt, 0);
        }
    }

    genename.insert(genename.end(), m_vecgenename.begin(), m_vecgenename.end());
    return vecCellgem.size();
}

bool cellAdjust::addborder(unsigned int cid, vector<Point> &vecPoint, vector<Point> &border, vector<short> &vec_border)
{
    convexHull(vecPoint, border, true);
    if(border.size() < 3)
    {
        printf("err cid=%d\n", cid);
        return false;
    }

    int bsz = border.size();
    int i=0;

    if(bsz > BORDERCNT)
    {
        vector<Point> tmpborder;
        double epsilon = 0.01 * arcLength(border, true);
        approxPolyDP(border, tmpborder, epsilon, true);

        bsz = tmpborder.size();
        for(;i<bsz;i++)
        {
            vec_border.emplace_back(tmpborder[i].x - m_cell_arrayptr[cid].x);
            vec_border.emplace_back(tmpborder[i].y - m_cell_arrayptr[cid].y);
        }
    }
    else
    {
        for(;i<bsz;i++)
        {
            vec_border.emplace_back(border[i].x - m_cell_arrayptr[cid].x);
            vec_border.emplace_back(border[i].y - m_cell_arrayptr[cid].y);
        }
    }

    for(;i<BORDERCNT;i++) //不足补SHRT_MAX
    {
        vec_border.emplace_back(SHRT_MAX);
        vec_border.emplace_back(SHRT_MAX);
    }

    return true;
}

//void cellAdjust::writeCell(vector<Cell> &veccell, vector<DnbExpression> &vecDnb)
void cellAdjust::writeCell(Cell *cellptr, int cellcnt, DnbExpression *dnbptr, int dnbcnt)
{
    uint16_t gene_count, exp_count, dnb_count, area, cell_type_id;
    map<uint16_t, uint16_t> map_gene_cnt; //gid midcnt
    uint32_t offset = 0;
    vector<Point> vecPoint;
    vector<Point> border;

    vector<short> vec_border;
    vec_border.reserve(cellcnt*2*BORDERCNT);

    m_cgefwPtr->cell_num_ = cellcnt;
    //for(Cell &ce : veccell)
    for(int ci=0;ci<cellcnt;ci++)
    {
        Cell &ce = cellptr[ci];
        vecPoint.clear();
        border.clear();
        map_gene_cnt.clear();
        for(int i=0;i<ce.count;i++)
        {
            DnbExpression &dnb = dnbptr[ce.offset+i];
            if(map_gene_cnt.find(dnb.gene_id) != map_gene_cnt.end())
            {
                map_gene_cnt[dnb.gene_id] += dnb.count;
            }
            else
            {
                gene_count++;
                map_gene_cnt.emplace(dnb.gene_id, dnb.count);
            }
            exp_count += dnb.count;
            vecPoint.emplace_back(dnb.x,dnb.y);
        }

        bool ret = addborder(ce.cellid, vecPoint, border, vec_border);
        Moments mu = moments(border, true);
        area = mu.m00;

        auto itor = map_gene_cnt.begin();
        for(;itor != map_gene_cnt.end();itor++)
        {
            m_cgefwPtr->cell_exp_list_.emplace_back(itor->first, itor->second);
            if(m_map_gene.find(itor->first) == m_map_gene.end())
            {
                vector<GeneExpData> tvec;
                m_map_gene.emplace(itor->first, tvec);
            }
            m_map_gene[itor->first].emplace_back(ce.cellid, itor->second);
        }

        cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
        CellData cell = {
                ce.cellid,
                m_cell_arrayptr[ce.cellid].x,
                m_cell_arrayptr[ce.cellid].y,
                offset,
                gene_count,
                exp_count,
                ce.count,
                area,
                cell_type_id
        };
        offset += gene_count;

        m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
        m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

        m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
        m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

        m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
        m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

        m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
        m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

        m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
        m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

        m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, dnb_count);
        m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, dnb_count);

        m_cgefwPtr->expression_num_ += gene_count;
        m_cgefwPtr->exp_count_sum_ += exp_count;
        m_cgefwPtr->dnb_count_sum_ += dnb_count;
        m_cgefwPtr->area_sum_ += area;

        m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
    }

    int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), m_cgefwPtr->cell_num_, effective_rect);
    m_cgefwPtr->storeCell(m_block_size[2]*m_block_size[3], m_blkidxPtr, m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();
}

void cellAdjust::writeGene()
{
    timer st(__FUNCTION__);
    m_cgefwPtr->gene_num_ = m_map_gene.size();
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_ , sizeof(GeneData)));
    
    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count = 0;
    vector<GeneExpData> gene_exp_list;

    int i = 0;
    auto itor = m_map_gene.begin();
    for(;itor != m_map_gene.end();itor++,i++)
    {
        for(GeneExpData &gexp : itor->second)
        {
            gene_exp_list.emplace_back(gexp);
            max_MID_count = std::max(max_MID_count, gexp.count);
            m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, gexp.count);
            exp_count += gexp.count;
        }

        cell_count = itor->second.size();
        gene_data_list[i].cell_count = cell_count;
        gene_data_list[i].exp_count = exp_count;
        string &strgene = m_vecgenename[itor->first];
        memcpy(gene_data_list[i].gene_name, strgene.c_str(), strgene.length());
        gene_data_list[i].max_mid_count = max_MID_count;
        gene_data_list[i].offset = offset;
        offset += cell_count;

        min_exp_count = std::min(min_exp_count, exp_count);
        max_exp_count = std::max(max_exp_count, exp_count);
        min_cell_count = std::min(min_cell_count, cell_count);
        max_cell_count = std::max(max_cell_count, cell_count);
    }

    m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count,
                                    gene_data_list, gene_exp_list);
    free(gene_data_list);
}

void cellAdjust::writeCellAdjust(const string &outpath, Cell *cellptr, int cellcnt, DnbExpression *dnbptr, int dnbcnt)
{
    m_cgefwPtr = new CgefWriter();
    m_cgefwPtr->setOutput(outpath);
    writeCell(cellptr, cellcnt, dnbptr, dnbcnt);
    writeGene();
    delete m_cgefwPtr;
}
