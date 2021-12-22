/** @file utils.h
    @brief Collection of utility functions.

    Created by huangzhibo on 2021/12/17.
*/

#ifndef GEFTOOLS__UTILS_H_
#define GEFTOOLS__UTILS_H_

#include <string>
#include <ctime>
#include <iostream>
#include <fstream>
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

/**
 * \brief Define a String type with 32 byte
 *
 * It is compatible with HDF5 H5T_C_S1 type.
 *   strtype = H5Tcopy(H5T_C_S1);
 *   H5Tset_size(strtype, 32);
 */
struct S32 {
    S32() = default;

    explicit S32(const char *c) {
        int i = 0;
        while (c[i] != '\0') {
            value[i] = c[i];
            ++i;
        }
    }

    char value[32] = {0};
};

/**
 * @brief Gets the format string of now time.
 */
S32 getStrfTime();

/*!
 * \brief Copies one file to another path
 * \param src_file The source file path
 * \param dst_file The destination file path
 * \return true if the file was copied successfully
 */
bool copyFile(const string& src_file, const string& dst_file);

/**
 *
 * @param coordinates
 * @param new_coordinates
 * @param offset_point
 */
void offsetCoordinates(vector<Point> & coordinates, vector<Point> & new_coordinates, Point & offset_point);

#endif //GEFTOOLS__UTILS_H_
