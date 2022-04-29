/*
 * @Author: zhaozijian
 * @Date: 2022-04-28 10:49:41
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-28 15:54:04
 * @Description: file content
 */
#ifndef GEFTOOLS_FILEREADER_H
#define GEFTOOLS_FILEREADER_H

#include <zlib.h>
#include "cgefParam.h"

class FileReader
{
public:
    FileReader(std::string &strpath);
    ~FileReader();
    void readfile();
    void readbuf();
    int cuttail(Membuf *ptr, int len);
private:
    gzFile m_infile;
    string m_leftstr;
};


#endif