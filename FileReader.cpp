/*
 * @Author: zhaozijian
 * @Date: 2022-04-28 10:50:02
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 11:03:01
 * @Description: file content
 */
#include "FileReader.h"
#include "readCellgemTask.h"

FileReader::FileReader(std::string &strpath)
{
    m_infile = gzopen(strpath.c_str(), "r");
    gzbuffer(m_infile, READLEN);
}

FileReader::~FileReader()
{
    gzclose(m_infile);
}

void FileReader::readfile()
{
    char buf[256]={0};
    while (1)
    {
        gzgets(m_infile, buf, 256);
        if(memcmp(buf, "geneID", 6) == 0) break; //跳过前面的标题
    }

    while (!gzeof(m_infile))
    {
        readbuf();
    }
}

void FileReader::readbuf()
{
    Membuf *ptr = new Membuf();
    char *pbuf = ptr->m_buf;
    int leftlen = m_leftstr.length();
    memcpy(pbuf, m_leftstr.c_str(), leftlen);
    m_leftstr.clear();
    pbuf += leftlen;
    int readlen = BUFLEN-leftlen;
    int reallen = gzread(m_infile, pbuf, readlen);
    
    cuttail(ptr, reallen+leftlen);
    cgefParam::GetInstance()->m_bpPtr->addptr(ptr);
}

int FileReader::cuttail(Membuf *ptr, int len)
{
    int i = len-1;
    for(;i>0;i--)
    {
        if(ptr->m_buf[i] == '\n')
        {
            break;
        }
    }

    ptr->m_len = i+1;
    m_leftstr.append(ptr->m_buf+ptr->m_len, len-ptr->m_len);
    return 0;
}