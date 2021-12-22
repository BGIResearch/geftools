/** @file main_convert.h
    @brief Main entrance of the geftools convert command.

    Created by huangzhibo on 2021/12/14.
*/

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

/**
 * @brief Main entrance of the geftools convert command.
 *
 * Command line arguments for geftools convert will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int convert(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_CONVERT_H_
