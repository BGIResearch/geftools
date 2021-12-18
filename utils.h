/** @file utils.h
    @brief Collection of utility functions.

    Created by huangzhibo on 2021/12/17.
*/

#ifndef GEFTOOLS__UTILS_H_
#define GEFTOOLS__UTILS_H_

#include <string>

using namespace std;

/*!
 * \fn static bool copyFile(const string& inPath, const string& outPath)
 * \brief Copies one file to another path
 * \param src_file The source file path
 * \param dst_file The destination file path
 * \return true if the file was copied successfully
 */
bool copyFile(const string& src_file, const string& dst_file);


#endif //GEFTOOLS__UTILS_H_
