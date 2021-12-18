//
// Created by huangzhibo on 2021/12/16.
//

#ifndef GEFTOOLS__MAIN_VIEW_H_
#define GEFTOOLS__MAIN_VIEW_H_

#include <string>

using namespace std;

struct ViewOptions {
    string input_file;
    string mask_file;
    string output_file;
    int threads;
    bool verbose;
};

int view(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_VIEW_H_
