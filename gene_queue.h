
#ifndef GENETOH5_GENEQUEUE_H
#define GENETOH5_GENEQUEUE_H

#include <thread>
#include <mutex>
#include <unistd.h>
#include <queue>
#include <iostream>
#include <condition_variable>
#include "gef.h"

class GeneQueue
{
public:
    GeneQueue(){};
    ~GeneQueue(){};
    void addqueue(GeneInfo *ptr)
    {
        std::lock_guard<std::mutex> tlock(m_mtx_queue);
        m_qgeneptr.emplace(ptr);
        m_cv_queue.notify_one();
    }

    GeneInfo* getGeneInfo2()
    {
        GeneInfo *ptr = nullptr;
        std::unique_lock<std::mutex> tlock(m_mtx_queue);

        m_cv_queue.wait(tlock, [this] {return !m_qgeneptr.empty();});
        ptr = m_qgeneptr.front();
        m_qgeneptr.pop();
        return ptr;
    }
private:
    std::mutex m_mtx_queue; //队列锁
    std::condition_variable m_cv_queue;
    std::queue<GeneInfo *> m_qgeneptr;
};


#endif