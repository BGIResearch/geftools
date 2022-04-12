/*
 * @Author: zhaozijian
 * @Date: 2022-04-07 15:58:56
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-07 18:11:48
 * @Description: file content
 */
#ifndef GEFTOOLS_TIMER_H
#define GEFTOOLS_TIMER_H

#include <chrono>
using namespace std::chrono;

class timer
{
public:
    timer(const char *msg):m_pmsg(msg)
    {
        m_start = system_clock::now();
        m_mid1 = m_start;
    }
    ~timer()
    {
        stop();
    }

    void start()
    {
        m_start = system_clock::now();
        m_mid1 = m_start;
    }

    void stop(const char *msg="-")
    {
        m_end = system_clock::now();
        duration<double> diff = m_end - m_start;
        printf("%s %s elapsed time: %7.5f ms\n", m_pmsg, msg, diff.count()*1000);
    }

    void showgap(const char *msg="-")
    {
        m_mid2 = system_clock::now();
        duration<double> diff = m_mid2 - m_mid1;
        m_mid1 = m_mid2;
        printf("%s %s elapsed time: %7.5f ms\n", m_pmsg, msg, diff.count()*1000);
    }
private:
    const char *m_pmsg;
    time_point<system_clock> m_start, m_mid1, m_mid2, m_end;
};
#endif