/*
 * @Author: zhaozijian
 * @Date: 2022-04-28 10:57:26
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 11:04:25
 * @Description: file content
 */
#ifndef GEFTOOLS_BUFPOOL_H
#define GEFTOOLS_BUFPOOL_H

#include "concurrentqueue/concurrentqueue.h"
// #include <list>
// #include <mutex>
// #include <condition_variable>
const int BUFLEN = 16*1024*1024;

class Membuf
{
public:
    Membuf()
    {
        m_buf = new char[BUFLEN];
    }
    ~Membuf()
    {
        delete[] m_buf;
    }
    int m_len;
    char *m_buf;
};


class BufPool
{
public:
    BufPool(int size)
    {
        m_emptyptr = new moodycamel::ConcurrentQueue<Membuf *>(size);
        // for(int i=0;i<size;i++)
        // {
        //     Membuf *ptr = new Membuf();
        //     m_emptyptr->enqueue(ptr);
        // }
        
        // m_fullptr = new moodycamel::ConcurrentQueue<Membuf *>(size);
        // m_totalbuf = size;
        // for(int i=0;i<m_totalbuf;i++)
        // {
        //     Membuf *ptr = new Membuf();
        //     m_emptyList.push_back(ptr);
        // }
    }
    ~BufPool()
    {
        // Membuf *ptr = nullptr;
        // bool ret = false
        // do
        // {
        //     ret = m_emptyptr->try_dequeue(ptr);
        // }while (!ret)
        
    }
    void addempty(Membuf *ptr)
    {
        // bool ret = false;
        // do
        // {
        //     ret = m_emptyptr->enqueue(ptr);
        // }while (!ret);

        // std::lock_guard<std::mutex> tlock(m_mtx_empty);
        // m_emptyList.push_back(ptr);
        // m_cv_empty.notify_one();
    }
    void addfull(Membuf *ptr)
    {
        // bool ret = false;
        // do
        // {
        //     ret = m_fullptr->enqueue(ptr);
        // }while (!ret);

        // std::lock_guard<std::mutex> tlock(m_mtx_full);
        // m_fullList.push_back(ptr);
        // m_cv_full.notify_one();
    }
    Membuf * getempty()
    {
        // Membuf *ptr = nullptr;
        // bool ret = false;
        // do
        // {
        //     ret = m_emptyptr->try_dequeue(ptr);
        // }while (!ret);
        // return ptr;

        // Membuf *ptr = nullptr;
        // std::unique_lock<std::mutex> tlock(m_mtx_empty);
        // while (m_emptyList.empty())
        // {
        //     //printf("%d getEmptyJob wait\n", syscall(SYS_gettid));
        //     m_cv_empty.wait(tlock);
        // }

        // //printf("%d getEmptyJob wait a buf\n", getpid());
        // ptr = m_emptyList.front();
        // m_emptyList.pop_front();
        //return ptr;
    }
    Membuf * getfull()
    {
        // Membuf *ptr = nullptr;
        // bool ret = false;
        // do
        // {
        //     ret = m_fullptr->try_dequeue(ptr);
        // }while (!ret);
        // return ptr;

        // Membuf *ptr = nullptr;
        // std::unique_lock<std::mutex> tlock(m_mtx_full);

        // while (m_fullList.empty())
        // {
        //     //printf("%d getFullJob wait\n", syscall(SYS_gettid));
        //     m_cv_full.wait(tlock);
        // }

        // //printf("%d getFullJob front\n", getpid());
        // ptr = m_fullList.front();
        // m_fullList.pop_front();
        // return ptr;
    }

    void addptr(Membuf *ptr)
    {
        m_emptyptr->enqueue(ptr);
    }
    Membuf *getptr()
    {
        Membuf *ptr = nullptr;
        bool ret = false;
        do
        {
            ret = m_emptyptr->try_dequeue(ptr);
        }while (!ret);
        return ptr;
    }
private:
    moodycamel::ConcurrentQueue<Membuf *> *m_emptyptr = nullptr;
    // moodycamel::ConcurrentQueue<Membuf *> *m_fullptr = nullptr;

    // int m_limit = 0;
    // int m_totalbuf = 0;
    // std::list<Membuf *> m_fullList;   //保存内容为满的内存指针
    // std::list<Membuf *> m_emptyList;  //保存内容为空的内存指针
    // std::mutex m_mtx_empty;
    // std::condition_variable m_cv_empty;

    // std::mutex m_mtx_full;
    // std::condition_variable m_cv_full;
};



#endif