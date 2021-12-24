/** @file main_view.h
    @brief Not yet implemented.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_VIEW_H_
#define GEFTOOLS__MAIN_VIEW_H_

#include <string>

using namespace std;

/**
 * @brief Parameters parsed from command lines.
 */
struct ViewOptions {
    string input_file;
    string output_gem;
    string output_mask;
    int threads;
    bool verbose;
};

/**
 * @brief Not yet implemented
 * @param argc
 * @param argv
 * @return
 */
int view(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_VIEW_H_
