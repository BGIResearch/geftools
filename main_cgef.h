/** @file main_cgef.h
    @brief Main entrance of the geftools cgef command.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_CGEF_H_
#define GEFTOOLS__MAIN_CGEF_H_

#include <string>
#include "cell_bin_writer.h"
#include "cxxopts.h"

using namespace std;

struct cgefOptions {
  string input_file;
  string mask_file;
  string output_file;
  int threads;
  bool verbose;
};

/**
 * @brief Main entrance of the geftools cgef command.
 *
 * Command line arguments for geftools cgef will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int cgef(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_CGEF_H_
