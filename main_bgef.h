/** @file main_cgef.h
    @brief Main entrance of the geftools cgef command.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_BGEF_H_
#define GEFTOOLS__MAIN_BGEF_H_

#include <string>
#include <thread>
#include <zlib.h>
#include <queue>
#include "bgef_writer.h"
#include "cxxopts.h"
#include "bgef_options.h"
#include "utils.h"

using namespace std;

//struct BgefOptions {
//  string input_file;
//  string output_file;
//  vector<int> bin_sizes;
//  int offset_x;
//  int offset_y;
//  int threads;
//  bool verbose;
//};

/**
 * @brief Main entrance of the geftools cgef command.
 *
 * Command line arguments for geftools cgef will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int bgef(int argc, char *argv[]);

/**
 * @brief Function to generate cgef (cell bin GEF) file.
 * @param cgef_file  The output file path of cgef (cell bin GEF).
 * @param bgef_file  The input file path of common bin GEF, only bin1 will be used.
 * @param mask_file  The input file path of mask file.
 * @param block_size The block size of cell dataset in cgef, including two elements, 0: x_block_size, 1: y_block_size
 * @param verbose    Print the run time of this function. default : false.
 * @return
 */
int generateBgef(const string& cgef_file,
                 const string& bgef_file,
                 const string& mask_file,
                 const int* block_size,
                 bool verbose = false);

/**
 * @brief Function to generate common bin GEF file(.bgef).
 * @param input_file  The input file path of gem file or bin1 bgef.
 * @param bgef_file  The output file path of common bin GEF (.bgef).
 * @param bin_size
 * @param n_thread
 * @param verbose
 * @return
 */
int generateBgef(const string &input_file,
                 const string &bgef_file,
                 int bin_size,
                 int n_thread,
                 bool verbose);

void gem2gef(BgefOptions *opts);

int mRead(BgefOptions *opts);

unsigned int parseResolutin(string& filename);

void writednb(BgefOptions *opts, BgefWriter &bgef_writer, int bin);

#endif //GEFTOOLS__MAIN_BGEF_H_
