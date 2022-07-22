/*
 * @Author: zhaozijian
 * @Date: 2022-05-16 11:02:23
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 11:35:58
 * @Description: file content
 */
#include "cellAdjust.h"
#include "timer.h"
#include <sstream>
#include <fstream>
#include "bgef_writer.h"
#include "dnb_merge_task.h"
#include "bin_task.h"


cellAdjust::cellAdjust()
{
}

cellAdjust::~cellAdjust()
{
    delete [] m_cell_arrayptr;
    H5Fclose(m_bgeffile_id);
}

void cellAdjust::readBgef(const string &strinput)
{
    timer st(__FUNCTION__);
    m_bgeffile_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/gene", H5P_DEFAULT);
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
    H5Sclose(gene_sid);
    H5Dclose(gene_did);


    hid_t exp_did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneexpcnt = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    Expression *expPtr = (Expression *) calloc(dims[0], sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expPtr);

    if(H5Lexists(m_bgeffile_id, "/geneExp/bin1/exon", H5P_DEFAULT)>0)
    {
        m_bexon = true;
        hsize_t edims[1];
        hid_t did = H5Dopen(m_bgeffile_id, "/geneExp/bin1/exon", H5P_DEFAULT);
        hid_t sid = H5Dget_space(did);
        H5Sget_simple_extent_dims(sid, edims, nullptr);
        assert(edims[0] == m_geneexpcnt);
        unsigned int *exonPtr = new unsigned int[edims[0]];
        H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, exonPtr);
        H5Sclose(sid);
        H5Dclose(did);
        for(int i=0;i<m_geneexpcnt;i++)
        {
            expPtr[i].exon = exonPtr[i];
        }
        delete []exonPtr;
    }

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

    if(H5Aexists(m_bgeffile_id, "omics"))
    {
        hid_t fattr = H5Aopen(m_bgeffile_id, "omics", H5P_DEFAULT);
        H5Aread(fattr, strtype, m_szomics);
    }
    H5Tclose(strtype);

    uint64_t l_id = 0;
    for(int i=0;i<m_genencnt;i++)
    {
        m_vecgenename.emplace_back(genePtr[i].gene);
        Expression *ptr = expPtr + genePtr[i].offset;
        for(int j=0;j<genePtr[i].count;j++)
        {
            l_id = ptr[j].x;
            l_id = (l_id<<32) | ptr[j].y;
            
            if(m_hash_vecdnb_exon.find(l_id) == m_hash_vecdnb_exon.end())
            {
                vector<Dnbs_exon> tvec;
                m_hash_vecdnb_exon.emplace(l_id, tvec);
            }
            m_hash_vecdnb_exon[l_id].emplace_back(i, ptr[j].count, ptr[j].exon);
        }
    }
    printf("gene:%d geneexp:%d hashcnt:%d\n", m_genencnt, m_geneexpcnt, m_hash_vecdnb_exon.size());
    free(genePtr);
    free(expPtr);
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

    // hid_t d_id = H5Dopen(file_id, "/cellBin/blockIndex", H5P_DEFAULT);
    // hid_t s_id = H5Dget_space(d_id);
    // H5Sget_simple_extent_dims(s_id, dims, nullptr);
    // m_blkidxPtr = static_cast<unsigned int *>(calloc(dims[0], sizeof(unsigned int)));
    // H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_blkidxPtr);

    // H5Sclose(s_id);
    // H5Dclose(d_id);

    hid_t d_id = H5Dopen(file_id, "/cellBin/blockSize", H5P_DEFAULT);
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
    m_fill_points = Mat::zeros(m_max_y-m_min_y+1, m_max_x-m_min_x+1, CV_8UC1);
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
            printf("empty cid %d\n",i);
        }
        ptmp += BORDERCNT*2;
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

    attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);

    attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    printf("offsetx:%d offsety:%d\n", m_offsetX, m_offsetY);
    H5Aclose(attr);
    H5Sclose(border_sid);
    H5Dclose(border_id);
    H5Fclose(file_id);

}

uint32_t cellAdjust::getCellLabelgem(vector<string> &genename, vector<cellgem_label> &vecCellgem)
{
    timer st(__FUNCTION__);
    genename.reserve(m_vecgenename.size());
    vecCellgem.reserve(m_geneexpcnt);

    uint64_t l_id = 0;
    vector<Point> vecpoint;
    int x,y;
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
            auto dnb_itor = m_hash_vecdnb_exon.find(l_id);
            if(dnb_itor!= m_hash_vecdnb_exon.end())
            {
                for(Dnbs_exon &dnbs : dnb_itor->second)
                {
                    vecCellgem.emplace_back(dnbs.geneid, x, y, dnbs.midcnt, itor->first+1);
                }
                m_hash_vecdnb_exon.erase(l_id);
            }
        }
    }
    printf("%d\n", m_hash_vecdnb_exon.size());

    auto itor_s = m_hash_vecdnb_exon.begin();
    for(;itor_s != m_hash_vecdnb_exon.end();itor_s++)
    {
        x = (itor_s->first) >> 32;
        y = (itor_s->first) & 0xFFFFFFFF;
        for(Dnbs_exon &dnbs : itor_s->second)
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
        //printf("borderwarn cid=%d pot=%d bor=%d dnb=%d\n", cid, vecPoint.size(), border.size(), m_cell_arrayptr[cid].dnb_count);
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

void cellAdjust::writeCell(Cell *cellptr, unsigned int cellcnt, DnbExpression *dnbptr, unsigned int dnbcnt)
{
    timer st(__FUNCTION__);
    uint16_t gene_count, exp_count, dnb_count, area, cell_type_id, maxExpmid = 0;
    map<uint16_t, uint16_t> map_gene_cnt; //gid midcnt
    uint32_t offset = 0;
    vector<Point> vecPoint;
    vector<Point> border;
    printf("rawcellcnt:%d newcellcnt:%d dnbcnt:%d\n", m_cellcnt, cellcnt, dnbcnt);

    unsigned int blocknum = m_block_size[2] * m_block_size[3];
    vector<vector<Cell>> vec_vec_cell;
    for(int i=0;i<blocknum;i++)
    {
        vector<Cell> vectmp;
        vec_vec_cell.emplace_back(std::move(vectmp));
    }
    unsigned int blkid = 0, cid=0;
    for(unsigned int ci=0;ci<cellcnt;ci++)
    {
        cid = cellptr[ci].cellid-1;
        CellData &cd = m_cell_arrayptr[cid];
        blkid = cd.x/m_block_size[0] + (cd.y/m_block_size[1])*m_block_size[2];
        vec_vec_cell[blkid].emplace_back(cellptr[ci]);
    }

    vector<uint32_t> vec_blkidx;
    vec_blkidx.reserve(blocknum+1);
    vector<short> vec_border;
    vec_border.reserve(cellcnt*2*BORDERCNT);

    uint32_t newcid = 0, offcnt = 0, blkcnt = 0;
    int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
    for(vector<Cell> &vecC : vec_vec_cell)
    {
        blkcnt = 0;
        for(Cell &ce : vecC)
        {
            cid = ce.cellid-1;
            vecPoint.clear();
            border.clear();
            gene_count = 0;
            exp_count = 0;
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

            bool ret = addborder(cid, vecPoint, border, vec_border);
            if(!ret) continue;
            
            Moments mu = moments(border, true);
            Rect trect = boundingRect(border);

            minx = std::min(minx, trect.x);
            maxx = std::max(maxx, trect.x+trect.width);
            miny = std::min(miny, trect.y);
            maxy = std::max(maxy, trect.y+trect.height);

            area = mu.m00;
            auto itor = map_gene_cnt.begin();
            for(;itor != map_gene_cnt.end();itor++)
            {
                m_cgefwPtr->cell_exp_list_.emplace_back(itor->first, itor->second);
                maxExpmid = std::max(maxExpmid, itor->second);
                if(m_map_gene.find(itor->first) == m_map_gene.end())
                {
                    vector<GeneExpData> tvec;
                    m_map_gene.emplace(itor->first, tvec);
                }
                m_map_gene[itor->first].emplace_back(newcid, itor->second);
            }

            cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
            CellData cell = {
                    newcid++,
                    m_cell_arrayptr[cid].x,
                    m_cell_arrayptr[cid].y,
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

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, cell.dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, cell.dnb_count);

            m_cgefwPtr->expression_num_ += gene_count;
            m_cgefwPtr->exp_count_sum_ += exp_count;
            m_cgefwPtr->dnb_count_sum_ += cell.dnb_count;
            m_cgefwPtr->area_sum_ += area;
            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
            blkcnt++;
        }

        vec_blkidx.emplace_back(offcnt);
        offcnt += blkcnt;
    }

    vec_blkidx.emplace_back(newcid);
    m_cgefwPtr->cell_num_ = newcid;
    m_cgefwPtr->max_mid_count_ = maxExpmid;

    int effective_rect[4] ={minx, miny, maxx, maxy};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), m_cgefwPtr->cell_num_, effective_rect);
    m_cgefwPtr->storeCell(m_block_size[2]*m_block_size[3], vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();
}

void cellAdjust::writeGene()
{
    timer st(__FUNCTION__);
    printf("genecnt:%d hashcnt:%d geneexpcnt:%d\n", m_genencnt, m_map_gene.size(), m_cgefwPtr->expression_num_);
    m_cgefwPtr->gene_num_ = m_genencnt;
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_ , sizeof(GeneData)));
    
    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count = 0;
    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(m_cgefwPtr->expression_num_);

    m_cgefwPtr->max_mid_count_ = 0;
    for(int i=0;i<m_genencnt;i++)
    {
        exp_count = 0;
        max_MID_count = 0;
        auto itor = m_map_gene.find(i);
        string &strgene = m_vecgenename[i];
        if(itor != m_map_gene.end())
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
            memcpy(gene_data_list[i].gene_name, strgene.c_str(), strgene.length());
            gene_data_list[i].max_mid_count = max_MID_count;
            gene_data_list[i].offset = offset;
            offset += cell_count;
        }
        else
        {
            memcpy(gene_data_list[i].gene_name, strgene.c_str(), strgene.length());
            gene_data_list[i].cell_count = 0;
            gene_data_list[i].exp_count = 0;
            gene_data_list[i].max_mid_count = 0; 
            gene_data_list[i].offset = 0;
        }

        m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, max_MID_count);
        min_exp_count = std::min(min_exp_count, exp_count);
        max_exp_count = std::max(max_exp_count, exp_count);
        min_cell_count = std::min(min_cell_count, cell_count);
        max_cell_count = std::max(max_cell_count, cell_count);
    }

    m_cgefwPtr->expression_num_ = gene_exp_list.size();
    m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count,
                                    gene_data_list, gene_exp_list);
    free(gene_data_list);
}

void cellAdjust::writeCellAdjust(const string &outpath, Cell *cellptr, int cellcnt, DnbExpression *dnbptr, int dnbcnt)
{
    m_cgefwPtr = new CgefWriter();
    m_cgefwPtr->setOutput(outpath);
    CellBinAttr cell_bin_attr = {
            .version = 2,
            .resolution = m_resolution,
            .offsetX = m_min_x,
            .offsetY = m_min_y
    };
    m_cgefwPtr->storeAttr(cell_bin_attr);
    writeCell(cellptr, cellcnt, dnbptr, dnbcnt);
    writeGene();
    delete m_cgefwPtr;
}

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    auto group_names = reinterpret_cast<std::vector<std::string>*>(opdata);
    group_names->push_back(name);
    return 0;
}

void cellAdjust::createRegionGef(const string &out)
{
    timer st(__FUNCTION__);
    hid_t gid = H5Gopen(m_bgeffile_id,"/geneExp", H5P_DEFAULT);
    std::vector<std::string> group_names;
    herr_t idx = H5Literate(gid, H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);
    H5Gclose(gid);

    for(string &str : group_names)
    {
        int bin = std::stoi(str.substr(3));
        m_bgefopts->bin_sizes_.push_back(bin);
    }

    m_bgefopts->m_genes_queue.init(m_bgefopts->map_gene_exp_.size());
    ThreadPool thpool(m_bgefopts->thread_ * 2);

    m_bgefopts->m_stromics.append(m_szomics);
    BgefWriter bgef_writer(out, false, m_bexon, m_bgefopts->m_stromics);
    bgef_writer.setResolution(m_resolution);

    int genecnt = 0;
    for(unsigned int bin : m_bgefopts->bin_sizes_)
    {
        auto& dnb_matrix = m_bgefopts->dnbmatrix_;
        auto& dnbAttr = m_bgefopts->dnbmatrix_.dnb_attr;

        dnbAttr.min_x = (m_min_x / bin) * bin;
        dnbAttr.len_x = (m_maxx)/bin + 1;
        dnbAttr.min_y = (m_min_y / bin) * bin;
        dnbAttr.len_y = (m_maxy)/bin + 1;
        dnbAttr.max_gene = 0;
        dnbAttr.max_mid = 0;
        dnbAttr.number = 0;
        unsigned long matrix_len = (unsigned long)(dnbAttr.len_x) * dnbAttr.len_y;
        printf("bin %d matrix: min_x=%d len_x=%d min_y=%d len_y=%d matrix_len=%lu\n",
               bin, dnbAttr.min_x, dnbAttr.len_x, dnbAttr.min_y, dnbAttr.len_y, matrix_len);
        if (bin == 1)
        {
            dnb_matrix.pmatrix_us = (BinStatUS*)calloc(matrix_len, sizeof(BinStatUS));
            assert(dnb_matrix.pmatrix_us);
            if(m_bexon)
            {
                dnb_matrix.pexon16 = (unsigned short*)calloc(matrix_len, 2);
                assert(dnb_matrix.pexon16);
            }
        }
        else
        {
            dnb_matrix.pmatrix = (BinStat*)calloc(matrix_len, sizeof(BinStat));
            assert(dnb_matrix.pmatrix);
            if(m_bexon)
            {
                dnb_matrix.pexon32 = (unsigned int*)calloc(matrix_len, 4);
                assert(dnb_matrix.pexon32);
            }
        }


        for(int i=0; i < m_bgefopts->thread_; i++)
        {
            auto *task = new DnbMergeTask(m_bgefopts->map_gene_exp_.size(), i, bin);
            thpool.addTask(task);
        }

        auto itor = m_bgefopts->map_gene_exp_.begin();
        for(;itor != m_bgefopts->map_gene_exp_.end();itor++)
        {
            auto *task = new BinTask(bin, itor->first.c_str());
            thpool.addTask(task);
        }

        unsigned int offset = 0;
        unsigned int maxexp = 0;
        unsigned int maxexon = 0;
        genecnt = 0;
        while (true)
        {
            GeneInfo *pgeneinfo = m_bgefopts->m_geneinfo_queue.getPtr();
            if (bin == 1){
                m_bgefopts->expressions_.insert(m_bgefopts->expressions_.end(), pgeneinfo->vecptr->begin(), pgeneinfo->vecptr->end());
            }
            else
            {
                for (auto g : *pgeneinfo->vecptr)
                {
                    g.x *= bin;
                    g.y *= bin;
                    m_bgefopts->expressions_.push_back(std::move(g));
                }
            }

            m_bgefopts->genes_.emplace_back(pgeneinfo->geneid, offset, static_cast<unsigned int>(pgeneinfo->vecptr->size()));
            offset += pgeneinfo->vecptr->size();
            maxexp = std::max(maxexp, pgeneinfo->maxexp);
            maxexon = std::max(maxexon, pgeneinfo->maxexon);
            
            if(bin == 100)
            {
                m_bgefopts->m_vec_bin100.emplace_back(pgeneinfo->geneid, pgeneinfo->umicnt, pgeneinfo->e10);
            }
            delete pgeneinfo;
            genecnt++;
            if(genecnt == m_bgefopts->map_gene_exp_.size())
            {
                break;
            }
        }

        bgef_writer.storeGene(m_bgefopts->expressions_, m_bgefopts->genes_, dnb_matrix.dnb_attr, maxexp, bin);
        bgef_writer.storeGeneExon(m_bgefopts->expressions_, maxexon, bin);
        m_bgefopts->expressions_.clear();
        m_bgefopts->genes_.clear();
        
        thpool.waitTaskDone();
        m_bgefopts->m_genes_queue.clear(bin);
        //write dnb
        if(bin == 100)
        {
            vector<GeneStat> &geneStat = m_bgefopts->m_vec_bin100;
            std::sort(geneStat.begin(), geneStat.end(), [](const GeneStat& p1, const GeneStat& p2){
                if (p1.mid_count > p2.mid_count)
                    return true;
                else if (p1.mid_count == p2.mid_count)
                {
                    int ret = strcmp(p1.gene, p2.gene);
                    return ret < 0;
                }
                else
                    return false;
            });
            bgef_writer.storeStat(geneStat);
        }

        vector<unsigned int> vec_mid;
        unsigned long number = 0;
        
        if (bin == 1)
        {
            for(unsigned long i=0;i<matrix_len;i++)
            {
                if(dnb_matrix.pmatrix_us[i].gene_count)
                {
                    ++number;
                    vec_mid.push_back(dnb_matrix.pmatrix_us[i].mid_count);
                }
            }
        }
        else
        {
            for(unsigned long i=0;i<matrix_len;i++)
            {
                if(dnb_matrix.pmatrix[i].gene_count)
                {
                    ++number;
                    vec_mid.push_back(dnb_matrix.pmatrix[i].mid_count);
                }
            }
        }

        int sz = vec_mid.size();
        sort(vec_mid.begin(), vec_mid.end(), [](const unsigned int a, const unsigned int b){return a<b;});
        if(bin > 50)
        {
            dnbAttr.max_mid = vec_mid[sz-1];
        }
        else
        {
            int limit = sz*0.999;
            dnbAttr.max_mid = vec_mid[limit];
        }
        
        dnbAttr.number = number;
        bgef_writer.storeDnb(dnb_matrix, bin);
        bgef_writer.storeWholeExon(dnb_matrix, bin);

        if (bin == 1)
        {
            if (dnb_matrix.pmatrix_us != nullptr)
            {
                free(dnb_matrix.pmatrix_us);
                dnb_matrix.pmatrix_us = nullptr;
                if(m_bexon)
                {
                    free(dnb_matrix.pexon16);
                    dnb_matrix.pexon16 = nullptr;
                }
            }
        }
        else
        {
            if (dnb_matrix.pmatrix != nullptr)
            {
                free(dnb_matrix.pmatrix);
                dnb_matrix.pmatrix = nullptr;
                if(m_bexon)
                {
                    free(dnb_matrix.pexon32);
                    dnb_matrix.pexon32 = nullptr;
                }
            }
        }
    }
}

void cellAdjust::getRegionGenedata(vector<vector<int>> &m_vecpos)
{
    timer st(__FUNCTION__);
    m_bgefopts = BgefOptions::GetInstance();
    int num = m_vecpos.size();
    uint64_t l_id = 0;
    vector<Point> non_zerovecpoint;
    vector<Point> relativepoint;
    int x,y;
    unsigned long totalSize = 0;

    for(int i=0;i<num;i++)
    {
        relativepoint.clear();
        non_zerovecpoint.clear();

        int cnt = m_vecpos[i].size();
        int *pos = m_vecpos[i].data();
        int minx = INT_MAX, miny = INT_MAX, maxx = 0, maxy = 0;
        for(int j=0;j<cnt;j+=2)
        {
            minx = std::min(minx, pos[j]);
            maxx = std::max(maxx, pos[j]);
            miny = std::min(miny, pos[j+1]);
            maxy = std::max(maxy, pos[j+1]);
        }
        m_maxx = std::max(m_maxx, maxx);
        m_maxy = std::max(m_maxy, maxy);

        for(int j=0;j<cnt;j+=2)
        {
            relativepoint.emplace_back(pos[j]-minx, pos[j+1]-miny);
        }

        int rows = maxy-miny+1;
        int cols = maxx-minx+1;
        Mat fill_points = Mat::zeros(rows, cols, CV_8UC1);
        fillPoly(fill_points, relativepoint, 1);
        findNonZero(fill_points, non_zerovecpoint);

        for(Point &pt : non_zerovecpoint)
        {
            x = pt.x+minx;
            y = pt.y+miny;
            l_id = x;
            l_id = (l_id << 32) | y;
            auto dnb_itor = m_hash_vecdnb_exon.find(l_id);
            if(dnb_itor!= m_hash_vecdnb_exon.end())
            {
                for(Dnbs_exon &dnbs : dnb_itor->second)
                {
                    string str(m_vecgenename[dnbs.geneid]);
                    auto itor_t = m_bgefopts->map_gene_exp_.find(str);
                    if(itor_t == m_bgefopts->map_gene_exp_.end())
                    {
                        vector<Expression> tvec;
                        m_bgefopts->map_gene_exp_.emplace(str, std::move(tvec));
                    }
                    m_bgefopts->map_gene_exp_[str].emplace_back(x,y,dnbs.midcnt, dnbs.exon);
                }
                m_hash_vecdnb_exon.erase(l_id);
                totalSize += dnb_itor->second.size();
            }
        }
    }

    m_bgefopts->expressions_.reserve(totalSize);
    m_bgefopts->genes_.reserve(m_bgefopts->map_gene_exp_.size());
}