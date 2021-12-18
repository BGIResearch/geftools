//
// Created by huangzhibo on 2021/12/15.
//

#ifndef GEFTOOLS__CELL_BIN_H_
#define GEFTOOLS__CELL_BIN_H_

#include <string>
#include <vector>
#include <iostream>
#include "hdf5.h"

using namespace std;

struct CellData {
  unsigned int x;
  unsigned int y;
  unsigned int offset;  //Offset of current cell in cellExp, 0-based
  unsigned short geneCount;
  unsigned short expCount;
  unsigned short dnbCount;
  unsigned short area;
  unsigned short cellTypeID;
};

struct GeneData {
    GeneData(const char* g, unsigned int m, unsigned int o, unsigned c)
    {
        int i = 0;
        while (g[i] != '\0')
        {
            geneName[i] = g[i];
            ++i;
        }
        maxMIDCount = m;
        offset = o;
        cellCount = c;
    }
    char geneName[32] = {0};
    unsigned int maxMIDCount;  //max MID count of current gene
    unsigned int offset;  //Offset of current gene in geneExp, 0-based
    unsigned int cellCount;
};

struct CellExpData {
    unsigned short geneID;
    unsigned short count;
};

struct GeneExpData {
    unsigned int cellID;
    unsigned short count;
};

struct CellBinAttr
{
    unsigned int version = 1;
    unsigned int resolution = 0;
    unsigned int offsetX = 0;
    unsigned int offsetY = 0;
};


class CellBin{
  private:
    hid_t file_id_;
    hid_t group_id_;
    void storeVersion();

  public:
    CellBin(const string& filepath,  const string& mode);
    ~CellBin();

    CellBinAttr cell_bin_attr;

    unsigned int gene_num_ = 0;
    unsigned int cell_num_ = 0;
    unsigned long long int expression_num_ = 0;

    unsigned int getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;

    void storeCell(unsigned int * x, unsigned int * y, unsigned short * area, unsigned int size);

    void storeCellExp();

    void storeCellBorder(char* borderPath, unsigned int size);

    void storeGeneList(vector<string>& geneList);

    void storeGeneList();

    void addDnbInCell(unsigned int * dnb_coordinates, unsigned int size);

    void setGeneExpMap(const string &inPath);

    void storeAttr();

};




#endif //GEFTOOLS__CELL_BIN_H_
