/** @file cell_bin_writer.h
    @brief Collection of functions for writing cell bin gef.

    Created by huangzhibo on 2021/12/16.
*/

#ifndef GEFTOOLS__CELL_BIN_WRITER_H_
#define GEFTOOLS__CELL_BIN_WRITER_H_

#include <string>
#include <vector>
#include <iostream>
#include "hdf5.h"

using namespace std;

/**
 * @brief Convert common bin gef to cell bin gef.
 * @param bin_gef common bin gef
 * @param mask_file binary image of cell shape
 * @param outfile output filepath of cell bin gef
 * @return
 */
int cell_bin_writer(const string & bin_gef, const string & mask_file, const string & outfile);

#endif //GEFTOOLS__CELL_BIN_WRITER_H_
