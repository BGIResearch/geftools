/** @file main_gen.h
    @brief Main entrance of the geftools gen command.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_GEN_H_
#define GEFTOOLS__MAIN_GEN_H_

#include <string>
#include "cell_bin_writer.h"
#include "cxxopts.h"

using namespace std;

struct genOptions {
  string input_file;
  string mask_file;
  string output_file;
  int threads;
  bool verbose;
};

/**
 * @brief Main entrance of the geftools gen command.
 *
 * Command line arguments for geftools gen will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int gen(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_GEN_H_
