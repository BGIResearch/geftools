
#ifndef GENETOH5_READTASK_H
#define GENETOH5_READTASK_H

#include <unordered_map>
#include "gef.h"
#include "utils.h"
#include "thread_pool.h"
#include "bgef_options.h"

using namespace std;

typedef struct
{
    int readlen; //想要读取的长度
    int reallen; //实际读取的长度
}RLen;


class ReadTask:public ITask
{
public:
    ReadTask();
    ~ReadTask();
    void doTask();
private:
    void readbuf(RLen &rlen);
    int cuttail(char *pbuf);
    int getGeneInfo();
    int mergeGeneinfo();
private:
    int m_buflen = 0;
    unsigned int min_x = UINT_MAX, min_y = UINT_MAX, max_x = 0, max_y = 0;
    char *m_pbuf = nullptr;
    BgefOptions *m_pcmd = nullptr;
    unordered_map<string, vector<Expression>> m_map_gege;

    static string m_leftstr;
    static mutex m_readmtx; //读文件锁
    static mutex m_mergemtx; //合并锁
};


#endif //GENETOH5_READTASK_H