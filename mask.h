/** @file mask.h
    @brief Declare a mask class to describe mask file info.

    The mask file is a binary image describing the shape of cells.
    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MASK_H_
#define GEFTOOLS__MASK_H_
#include <vector>
#include <string>
#include "polygon.h"

using namespace std;

/**
 * @brief A mask class to describe mask file info.
 */
class Mask {
  private:
    vector<Polygon> polygons_;
    unsigned int cell_num_;
    int rows_, cols_;

  public:
    explicit Mask(const string& file);


    void to_mask();
    void write_polygons();

    /**
     * @brief Get the number of cells, one polygon obtained from mask is one cell
     */
    unsigned int getCellNum() const;

    /**
     * Get polygon information of all cells
     * @return A vector of polygons
     */
    const vector<Polygon> &getPolygons() const;
};

#endif //GEFTOOLS__MASK_H_
