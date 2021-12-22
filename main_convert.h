//
// Created by huangzhibo on 2021/12/14.
//

#ifndef GEFTOOLS__MAIN_CONVERT_H_
#define GEFTOOLS__MAIN_CONVERT_H_

#include <string>
#include "cell_bin_writer.h"
#include "cxxopts.h"

using namespace std;

struct ConvertOptions {
  string input_file;
  string mask_file;
  string output_file;
  int threads;
  bool verbose;
};

int convert(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_CONVERT_H_
