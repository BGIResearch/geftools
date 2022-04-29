/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:30
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 15:35:18
 * @Description: file content
 */

#include "cgefCellgem.h"
#include "cgefParam.h"
#include "readCellgemTask.h"
#include "timer.h"
#include <functional>
#include "FileReader.h"

cgefCellgem::cgefCellgem(/* args */)
{
    m_thpoolPtr = new ThreadPool(cgefParam::GetInstance()->m_threadcnt);
}

cgefCellgem::~cgefCellgem()
{
    delete m_thpoolPtr;
}

size_t Rect_hash(const Rect &rect)
{
    unsigned long tmp = 0, tx = rect.x, ty = rect.y;
    tmp = (tx << 40) | (ty << 16) | (rect.width << 8) | rect.height;
    return tmp;
}

bool Rectequal_to(const Rect &a, const Rect &b)
{
    return (a.x == b.x) && (a.y == b.y) && (a.width == b.width) && (a.height == b.height);
}

void cgefCellgem::readmask()
{
    timer ct(__FUNCTION__);
    cv::Mat img = cv::imread(cgefParam::GetInstance()->m_maskstr, -1);
    assert(!img.empty());
    m_rows = img.rows;
    m_cols = img.cols;
    printf("img row:%d col:%d\n", m_rows, m_cols);

    m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
    m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
    m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
    m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
    m_blocknum = m_block_size[2] * m_block_size[3];

    vector<Vec4i> hierarchy;
    findContours(img, m_contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
    int scnt = m_contours.size();
    std::unordered_map<Rect, int, function<size_t (const Rect &)>, function<bool (const Rect &, const Rect &)>> map_cidx(scnt, Rect_hash, Rectequal_to);
    
    for(int i=0;i<scnt;i++)
    {
        const Rect &rect = cv::boundingRect(m_contours[i]);
        map_cidx.emplace(rect, i);
    }

    Mat outimg;
    int count = connectedComponentsWithStats(img, outimg, m_stats, m_centroids);
    assert(count == scnt+1);
    m_maskcellnum = scnt;

    m_vec_veccell.reserve(m_blocknum);
    for(int i=0;i<m_blocknum;i++)
    {
        vector<celldata> vectmp;
        m_vec_veccell.emplace_back(std::move(vectmp));
    }

    int cx, cy, x, y, w, h, c_idx, blkid;
    for(int i=1;i<count;i++)
    {
        cx = m_centroids.at<double>(i,0);
        cy = m_centroids.at<double>(i,1);

        x = m_stats.at<int>(i, CC_STAT_LEFT);
        y = m_stats.at<int>(i, CC_STAT_TOP);
        w = m_stats.at<int>(i, CC_STAT_WIDTH);
        h = m_stats.at<int>(i, CC_STAT_HEIGHT);
        Rect r1(x, y, w, h);

        m_min_x = std::min(m_min_x, x);
        m_max_x = std::max(m_max_x, x+w);
        m_min_y = std::min(m_min_y, y);
        m_max_y = std::max(m_max_y, y+h);

        c_idx = map_cidx.at(r1);
        blkid = cx/m_block_size[0] + (cy/m_block_size[1])*m_block_size[2];
        m_vec_veccell[blkid].emplace_back(c_idx, i);
    }
}

void cgefCellgem::readxy()
{
    timer st(__FUNCTION__);
    cgefParam::GetInstance()->m_infile = gzopen(cgefParam::GetInstance()->m_rawgemstr.c_str(), "r");
    gzbuffer(cgefParam::GetInstance()->m_infile, READLEN);

    std::string line;
    while (readline(cgefParam::GetInstance()->m_infile, line))
    {
        if (line.substr(0, 6) == "geneID") break;
    }

    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        readCellgemTask *rtask = new readCellgemTask(2);
        m_thpoolPtr->addTask(rtask);
    }
    m_thpoolPtr->waitTaskDone();
    gzclose(cgefParam::GetInstance()->m_infile);
    printf("minx:%d miny:%d\n", cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_min_y);
}


void cgefCellgem::readcellgem()
{
    timer st(__FUNCTION__);
    cgefParam::GetInstance()->m_infile = gzopen(cgefParam::GetInstance()->m_cellgemstr.c_str(), "r");
    gzbuffer(cgefParam::GetInstance()->m_infile, READLEN);

    std::string line;
    while (readline(cgefParam::GetInstance()->m_infile, line))
    {
        if (line.substr(0, 6) == "geneID") break;
    }

    int type = 1;
    if(cgefParam::GetInstance()->m_intype == INPUTTYPE_GEM_LABEL)
    {
        type = 3;
    }

    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        readCellgemTask *rtask = new readCellgemTask(type);
        m_thpoolPtr->addTask(rtask);
    }
    m_thpoolPtr->waitTaskDone();
    gzclose(cgefParam::GetInstance()->m_infile);

    printf("cellcnt:%ld genecnt:%ld \n", cgefParam::GetInstance()->m_map_cell.size(), 
                        cgefParam::GetInstance()->m_map_gene.size());
    auto itor = cgefParam::GetInstance()->m_map_gene.begin();
    int idx = 0;
    for(;itor!=cgefParam::GetInstance()->m_map_gene.end();itor++)
    {
        m_hash_gname2gid.emplace(itor->first, idx++);//gname到geneid的映射
    }
}

void cgefCellgem::readcellgem_new()
{
    timer st(__FUNCTION__);

    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        readCellgemTask *rtask = new readCellgemTask(4);
        m_thpoolPtr->addTask(rtask);
    }

    FileReader freader(cgefParam::GetInstance()->m_cellgemstr);
    freader.readfile();

    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        cgefParam::GetInstance()->m_bpPtr->addptr(nullptr);
    }
    

    m_thpoolPtr->waitTaskDone();

    printf("cellcnt:%ld genecnt:%ld \n", cgefParam::GetInstance()->m_map_cell.size(), 
                        cgefParam::GetInstance()->m_map_gene.size());
    auto itor = cgefParam::GetInstance()->m_map_gene.begin();
    int idx = 0;
    for(;itor!=cgefParam::GetInstance()->m_map_gene.end();itor++)
    {
        m_hash_gname2gid.emplace(itor->first, idx++);//gname到geneid的映射
    }
}

void cgefCellgem::writeFile(CgefWriter *cwptr)
{
    m_cgefwPtr = cwptr;
// readmask();
//     // readxy();
//     //     printf("minx:%d maxx:%d miny:%d maxy:%d\n", cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_max_x,
//     //                                 cgefParam::GetInstance()->m_min_y, cgefParam::GetInstance()->m_max_y);

//     long total = 22747;
//     total *= 25666;
//     cgefParam::GetInstance()->m_pdata = (char*)calloc(total, 1);
//     readxy();

//     cv::Mat img(25666, 22747, CV_8UC1, cgefParam::GetInstance()->m_pdata);

//     Mat outimg;
//     int count = connectedComponentsWithStats(img, outimg, m_stats, m_centroids);
// printf("%d\n", count);
//     return;

    if(cgefParam::GetInstance()->m_intype == INPUTTYPE_GEM_ADJUST)
    {
        readmask();
        m_cgefwPtr->m_x_len = m_cols;
        m_cgefwPtr->m_y_len = m_rows;
        readxy();
        readcellgem();
        writeAttr();
        writeCell();
        writeGene();
    }
    else if(cgefParam::GetInstance()->m_intype == INPUTTYPE_GEM_LABEL)
    {
        //readcellgem_new();
        
        readcellgem();
        printf("minx:%d maxx:%d miny:%d maxy:%d\n", cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_max_x,
                                    cgefParam::GetInstance()->m_min_y, cgefParam::GetInstance()->m_max_y);
        writeAttr();
        celltype();
        writeCell_celltype();
        writeGene();
    }
}

void cgefCellgem::writeAttr()
{
    CellBinAttr cell_bin_attr = {
            .version = 2,
            .resolution = 0,
            .offsetX = cgefParam::GetInstance()->m_min_x,
            .offsetY = cgefParam::GetInstance()->m_min_y
    };
    m_cgefwPtr->storeAttr(cell_bin_attr);
}

void cgefCellgem::addCellborder(int cx, int cy, vector<char> &vec_border, celldata & cdata)
{
    vector<Point> &vborder = m_contours[cdata.c_idx];
    vector<Point> tmpborder;
    int sz = vborder.size();
    if(vborder.size() > BORDERCNT)
    {
        double epsilon = 0.01 * arcLength(vborder, true);
        approxPolyDP(vborder, tmpborder, epsilon, true);
    }
    
    sz = tmpborder.size();
    int i=0;
    for(;i<sz;i++)
    {
        vec_border.emplace_back(vborder[i].x - cx);
        vec_border.emplace_back(vborder[i].y - cy);
    }

    for(;i<BORDERCNT;i++) //不足补0
    {
        vec_border.emplace_back(0);
        vec_border.emplace_back(0);
    }
}

void cgefCellgem::writeCell()
{
    timer st(__FUNCTION__);
    unsigned int cid = 0, gid = 0, offcnt = 0; 
    unsigned short gene_count, exp_count, dnb_count, area, cell_type_id;

    int cx, cy; //细胞质点
    vector<char> vec_border;
    vec_border.reserve(m_maskcellnum*2*16);
    vector<unsigned int> vec_blkidx;
    vec_blkidx.reserve(m_blocknum+1);
    vector<unsigned int> vec_cellLabel;
    vec_cellLabel.reserve(m_maskcellnum);

    for(vector<celldata> &vec : m_vec_veccell)
    {
        int cnt = 0;
        for(celldata & cdata : vec)
        {
            cgef_cell *cellptr = cgefParam::GetInstance()->m_map_cell[cdata.l_idx];
            if(cellptr == nullptr) continue;
 
            m_hash_clabel2cid.emplace(cdata.l_idx, cid);
            vec_cellLabel.emplace_back(cdata.l_idx);
            
            // vector<Point> &vborder = m_contours[cdata.c_idx];
            cx = m_centroids.at<double>(cdata.l_idx,0);
            cy = m_centroids.at<double>(cdata.l_idx,1);
            addCellborder(cx, cy, vec_border, cdata);
            // if(vborder.size() > 16)
            // {
            //     for(int i=0;i<16;i++)
            //     {
            //         vec_border.emplace_back(vborder[i].x - cx);
            //         vec_border.emplace_back(vborder[i].y - cy);
            //     }
            // }

            gene_count = cellptr->m_map_cellexp.size();
            exp_count = cellptr->expcnt;
            dnb_count = cellptr->dnbcnt;
            auto itor = cellptr->m_map_cellexp.begin();
            for(;itor != cellptr->m_map_cellexp.end();itor++)
            {
                gid = m_hash_gname2gid[itor->first];
                m_cgefwPtr->cell_exp_list_.emplace_back(gid, itor->second);
            }

            area = m_stats.at<int>(cdata.l_idx, CC_STAT_AREA);
            cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
            CellData cell = {
                    cid,
                    cx,
                    cy,
                    m_cgefwPtr->expression_num_,
                    gene_count,
                    exp_count,
                    dnb_count,
                    area,
                    cell_type_id
            };
            cnt++;
            cid++;

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
        vec_blkidx.emplace_back(offcnt);
        offcnt += cnt;
    }
    vec_blkidx.emplace_back(cid);
    
    m_cgefwPtr->cell_num_ = cid;

    unsigned int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);

    m_cgefwPtr->storeCell(m_blocknum, vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();
    m_cgefwPtr->storeCellLabel(vec_cellLabel);
}

void cgefCellgem::writeGene()
{
    timer st(__FUNCTION__);
    m_cgefwPtr->gene_num_ = cgefParam::GetInstance()->m_map_gene.size();
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_ , sizeof(GeneData)));
    
    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count = 0;
    vector<GeneExpData> gene_exp_list;

    int cid = 0, i = 0;
    auto itor = cgefParam::GetInstance()->m_map_gene.begin();
    for(;itor != cgefParam::GetInstance()->m_map_gene.end();itor++,i++)
    {
        cgef_gene *geneptr = itor->second;
        auto itor_g = geneptr->m_map_geneexp.begin();
        for(;itor_g != geneptr->m_map_geneexp.end();itor_g++)
        {
            cid = m_hash_clabel2cid[itor_g->first];
            gene_exp_list.emplace_back(cid, itor_g->second);
            max_MID_count = std::max(max_MID_count, itor_g->second);
            m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, itor_g->second);
        }

        cell_count = geneptr->m_map_geneexp.size();
        gene_data_list[i].cell_count = cell_count;
        gene_data_list[i].exp_count = geneptr->expcnt;
        memcpy(gene_data_list[i].gene_name, itor->first.c_str(), itor->first.length());
        gene_data_list[i].max_mid_count = max_MID_count;
        gene_data_list[i].offset = offset;
        offset += cell_count;

        min_exp_count = std::min(min_exp_count, geneptr->expcnt);
        max_exp_count = std::max(max_exp_count, geneptr->expcnt);
        min_cell_count = std::min(min_cell_count, cell_count);
        max_cell_count = std::max(max_cell_count, cell_count);
    }

    m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count,
                                    gene_data_list, gene_exp_list);
    free(gene_data_list);
}

void cgefCellgem::celltype()
{
    timer st(__FUNCTION__);
    
    m_rows = cgefParam::GetInstance()->m_max_y - cgefParam::GetInstance()->m_min_y;
    m_cols = cgefParam::GetInstance()->m_max_x - cgefParam::GetInstance()->m_min_x;
    m_cgefwPtr->m_x_len = m_cols;
    m_cgefwPtr->m_y_len = m_rows;

    m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
    m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
    m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
    m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
    m_blocknum = m_block_size[2] * m_block_size[3];
    m_vec_veccell.reserve(m_blocknum);
    for(int i=0;i<m_blocknum;i++)
    {
        vector<celldata> vectmp;
        m_vec_veccell.emplace_back(std::move(vectmp));
    }

    int celltypeid = 0;
    auto itor = cgefParam::GetInstance()->m_map_cell.begin();
    for(;itor != cgefParam::GetInstance()->m_map_cell.end();itor++)
    {
        itor->second->getCenter(m_block_size, cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_min_y);
        m_vec_veccell[itor->second->m_blkid].emplace_back(0, itor->first); //没有连通域
        assert(itor->first == itor->second->m_celllabel);
        if(m_hash_celltype.find(itor->second->m_type) == m_hash_celltype.end())
        {
            m_hash_celltype.emplace(itor->second->m_type, celltypeid++);
            S32 cell_type(itor->second->m_type);
            m_cgefwPtr->cell_type_list_.emplace_back(std::move(cell_type));
        }
    }
}

void cgefCellgem::writeCell_celltype()
{
    timer st(__FUNCTION__);
    unsigned int cid = 0, gid = 0, offcnt = 0; 
    unsigned short gene_count, exp_count, dnb_count, area, cell_type_id;

    int cx, cy; //细胞质点
    vector<unsigned int> vec_blkidx;
    vec_blkidx.reserve(m_blocknum+1);
    vector<unsigned int> vec_cellLabel;
    vec_cellLabel.reserve(m_maskcellnum);
    
    for(vector<celldata> &vec : m_vec_veccell)
    {
        int cnt = 0;
        for(celldata & cdata : vec)
        {
            cgef_cell *cellptr = cgefParam::GetInstance()->m_map_cell[cdata.l_idx];
            if(cellptr == nullptr) continue;
 
            cx = cellptr->m_centerx;
            cy = cellptr->m_centery;
            m_hash_clabel2cid.emplace(cdata.l_idx, cid);
            vec_cellLabel.emplace_back(cdata.l_idx);

            gene_count = cellptr->m_map_cellexp.size();
            exp_count = cellptr->expcnt;
            dnb_count = cellptr->dnbcnt;
            auto itor = cellptr->m_map_cellexp.begin();
            for(;itor != cellptr->m_map_cellexp.end();itor++)
            {
                gid = m_hash_gname2gid[itor->first];
                m_cgefwPtr->cell_exp_list_.emplace_back(gid, itor->second);
            }

            //area = m_stats.at<int>(cdata.l_idx, CC_STAT_AREA);
            //cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
            cell_type_id = m_hash_celltype[cellptr->m_type];
            CellData cell = {
                    cid,
                    cx,
                    cy,
                    m_cgefwPtr->expression_num_,
                    gene_count,
                    exp_count,
                    dnb_count,
                    0,
                    cell_type_id
            };
            cnt++;
            cid++;

            m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
            m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

            m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
            m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

            // m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
            // m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

            m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
            m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

            m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
            m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, dnb_count);

            m_cgefwPtr->expression_num_ += gene_count;
            m_cgefwPtr->exp_count_sum_ += exp_count;
            m_cgefwPtr->dnb_count_sum_ += dnb_count;
            //m_cgefwPtr->area_sum_ += area;

            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
        }
        vec_blkidx.emplace_back(offcnt);
        offcnt += cnt;
    }
    vec_blkidx.emplace_back(cid);
    
    m_cgefwPtr->cell_num_ = cid;

    // unsigned int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
    // m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);

    m_cgefwPtr->storeCell(m_blocknum, vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList_N();
    m_cgefwPtr->storeCellLabel(vec_cellLabel);

}

