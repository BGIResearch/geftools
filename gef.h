/** @file gef.h
    @brief Declare GEFs related structs.

    Created by huangzhibo on 2021/12/23.
*/

#ifndef GEFTOOLS_GEF_H
#define GEFTOOLS_GEF_H

#include "hdf5.h"
#include <string>
#include <cstring>
#include <vector>

using namespace std;

/**
 * @brief Coordinate union
 */
union Coordinate {
    unsigned int pos[2]; ///< dnb coordinates x, y
    unsigned long long int pos_id;
};

/**
 * @brief Expression struct
 */
struct Expression {
    Expression(int x, int y, unsigned int count, unsigned int exon=0):
    x(x),y(y),count(count), exon(exon){};
    int x; ///< dnb coordinates x
    int y; ///< dnb coordinates x
    unsigned int count; ///< expression count (MIDcount)
    unsigned int exon; //exon cnt
};


struct DnbExpression {
    int x; ///< dnb coordinates x
    int y; ///< dnb coordinates x
    unsigned short count; ///< expression count (MIDcount)
    unsigned short gene_id;
};

/**
 * @brief ExpressionAttr struct, record the attributes of the dataset named expressione
 */
struct ExpressionAttr
{
    int min_x; ///< Min X of dnb coordinate
    int min_y; ///< Min Y of dnb coordinate
    int max_x; ///< Max X of dnb coordinate
    int max_y; ///< Max Y of dnb coordinate
    unsigned int max_exp;  ///< Max expression count
    unsigned int resolution; ///< The resolution of stereo chip
};

struct Gene {
    Gene(const char* g, unsigned int o, unsigned c)
    {
        int i = 0;
        while (g[i] != '\0')
        {
            gene[i] = g[i];
            ++i;
        }
        offset = o;
        count = c;
    }
    char gene[32] = {0};
    unsigned int offset;
    unsigned int count;
};

struct Cell
{
    unsigned int cellid;
    unsigned int offset;
    unsigned short count;
};

struct GeneStat
{
    GeneStat(string& g, unsigned int m, float e)
    {
        memcpy(gene, g.c_str(), g.size());
        mid_count = m;
        E10 = e;
    }
    GeneStat(const char* g, unsigned int m, float e)
    {
        int len = strlen(g);
        memcpy(gene, g, len);
        mid_count = m;
        E10 = e;
    }
    char gene[32] = {0};
    unsigned int mid_count;
    float E10;
};

//simple gene data
struct GeneS 
{
    GeneS(const char *ptr):geneid(ptr){};
    const char *geneid;
    std::vector<Expression> *vecptr;
};

struct GeneInfo
{
    GeneInfo(const char *ptr):geneid(ptr),umicnt(0){};
    const char *geneid;
    unsigned long umicnt;
    float e10;
    //float c50;
    unsigned int maxexp;
    unsigned int maxexon;
    std::vector<Expression> *vecptr;
};

// wholeExp Matrix
struct BinStat {
    unsigned int mid_count;
    unsigned short gene_count;
};

struct BinStatUS {
    unsigned short mid_count;
    unsigned short gene_count;
};

struct DnbAttr {
    int min_x;
    int len_x;
    int min_y;
    int len_y;
    unsigned int max_mid;
    unsigned int max_gene;
    unsigned int max_exon;
    unsigned long number;
};

struct DnbMatrix {
    DnbAttr dnb_attr;
    BinStatUS *pmatrix_us;
    BinStat *pmatrix;
    unsigned short *pexon16;
    unsigned int *pexon32;
};

// struct GeneErank
// {
//     GeneErank(const char *ptr):geneid(ptr){};
//     const char *geneid;
//     unsigned long umicnt;
//     float e10;
//     //float c50;
//     char attribute[10];
// };

/**
 * @brief Describe the Cell dataset in the cell bin GEF file
 */
struct CellData {
    unsigned int id;
    int x; ///< Coordinate X of center point in this cell
    int y; ///< Coordinate Y of center point in this cell
    unsigned int offset;  ///< Offset of current cell in cellExp, 0-based
    unsigned short gene_count; ///< The number of gene in this cell
    unsigned short exp_count; ///< The total expression count of all genes in this cell
    unsigned short dnb_count; ///< Dnb number in this cell
    unsigned short area; ///< The polygon area of this cell
    unsigned short cell_type_id; ///< Cell type ID to index the CellTypeList
    unsigned short cluster_id; ///< Cluster ID, should start from 1
    //unsigned short incnt;
};

struct CellAttr {
    float average_gene_count;
    float average_exp_count;
    float average_dnb_count;
    float average_area;
    float median_gene_count;
    float median_exp_count;
    float median_dnb_count;
    float median_area;
    int min_x;
    int min_y;
    unsigned short min_gene_count;
    unsigned short min_exp_count;
    unsigned short min_dnb_count;
    unsigned short min_area;
    int max_x;
    int max_y;
    unsigned short max_gene_count;
    unsigned short max_exp_count;
    unsigned short max_dnb_count;
    unsigned short max_area;
};

struct GeneData {
    GeneData(const char* g, unsigned int o, unsigned int c, unsigned int e, unsigned m)
    {
        int i = 0;
        while (g[i] != '\0')
        {
            gene_name[i] = g[i];
            ++i;
        }
        offset = o;
        cell_count = c;
        exp_count = e;
        max_mid_count = m;
    }
    char gene_name[32] = {0};
    unsigned int offset;  ///< Offset of current gene in geneExp, 0-based
    unsigned int cell_count;
    unsigned int exp_count;
    unsigned short max_mid_count;  ///< max MID count of current gene
};

struct CellExpData {
//    explicit CellExpData(unsigned int gene_exp){
//        geneID = gene_exp >> 16;
//        count = gene_exp  & 0xFFFF;
//    }
//    CellExpData(unsigned short g, unsigned short c){
//        geneID = g;
//        count = c;
//    }
    CellExpData(unsigned short id, unsigned short cnt):gene_id(id),count(cnt){}
    unsigned short gene_id;
    unsigned short count;
};

struct GeneExpData {
    GeneExpData(unsigned int id, unsigned short cnt):cell_id(id),count(cnt){}
    unsigned int cell_id;
    unsigned short count;
};

/**
 * @brief Attributes of the cell bin GEF.
 */
struct CellBinAttr
{
    unsigned int version; ///< Cell Bin GEF version
    unsigned int resolution; ///< Pitch (nm) between neighbor spots
    int offsetX; ///< Minimum value of x-axis coordinate with offset
    int offsetY; ///< Minimum value of y-axis coordinate with offset
    string omics;
};

struct Dnbs
{
    Dnbs(uint16_t gid, uint16_t cnt):geneid(gid),midcnt(cnt){};
    uint16_t geneid;
    uint16_t midcnt;
};

struct Dnbs_exon
{
    Dnbs_exon(uint16_t gid, uint16_t cnt, uint16_t exon):geneid(gid),midcnt(cnt),exon(exon){};
    uint16_t geneid;
    uint16_t midcnt;
    uint16_t exon;
};


hid_t getMemtypeOfGeneData();
hid_t getMemtypeOfGeneExpData();
hid_t getMemtypeOfCellData();
hid_t getMemtypeOfCellExpData();

#endif //GEFTOOLS_GEF_H
