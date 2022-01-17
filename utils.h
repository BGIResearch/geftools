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
#include <array>
#include <iomanip>
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

static union
{
    char c[4];
    unsigned long l;
}endian_test = { { 'l','?','?','b' } };
#define ENDIANNESS ((char)endian_test.l)


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
 * @brief Split string based on a character delimiter.
 * @param s   The string to split.
 * @param delim  The delimiter, a character.
 * @return
 */
vector<string> split(const string &s, char delim);

/**
 * @brief Read lines from txt file.
 * @param filename
 * @return
 */
vector<string> readLines(const string &filename);

/**
 * @brief Print time.
 * @param prev   Previous time， ( the start time ).
 * @param message  Message of the step.
 * @return Current time, used as the next previous time.
 */
time_t printTime(time_t prev, string message);

/**
 * @brief Print cpu time.
 * @param prev   Previous clock time， ( the start time ).
 * @param message  Message of the step.
 * @return Current clock time, used as the next previous clock time.
 */
unsigned long printCpuTime(unsigned long prev, string message);

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
