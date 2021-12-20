/** @file utils.h
    @brief Collection of utility functions.

    Created by huangzhibo on 2021/12/17.
*/

#ifndef GEFTOOLS__UTILS_H_
#define GEFTOOLS__UTILS_H_

#include <string>
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

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
