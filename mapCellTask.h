/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 16:56:28
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-03-25 17:25:56
 * @Description: file content
 */
#ifndef GEFTOOLS_MAPCELL_H_
#define GEFTOOLS_MAPCELL_H_

#include "thread_pool.h"
#include "gef.h"
#include "cgefParam.h"

//pointPolygonTest 
class mapCellTask:public ITask
{

public:
    mapCellTask();
    ~mapCellTask();
    void doTask()
    {
        auto itor  = vec.begin();
        for(itor;itor != vec.end();itor++)
        {
            for(GeneExp &gxp : *itor)
            {
                if(pointPolygonTest())
            }
        }
    }
private:
    /* data */
};


#endif