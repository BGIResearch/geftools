/** @file mask.h
    @brief Declare a Mask class to describe mask file info.

    The mask file is a binary image describing the shape of cells.
    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MASK_H_
#define GEFTOOLS__MASK_H_
#include <string>
#include "polygon.h"
#include "utils.h"

using namespace std;

/**
 * @brief A mask class to describe mask file info.
 */
class Mask {
  private:
    unsigned int cell_num_;
    unsigned int block_num_;
    unsigned int block_size_[4] = {0};
    unsigned int * block_index_ = nullptr;
    vector<vector<Point> > contours_;
    vector<Vec4i> hierarchy_;
    vector<Polygon> polygons_;
    int min_x_{INT_MAX}, max_x_{0}, min_y_{INT_MAX}, max_y_{0}, rows_{0}, cols_{0};

  public:

    /**
     * @brief Constructor of Mask class.
     * @param file   Mask file of cell shape.
     */
    Mask(const string& file, const int block_size[]);
    /**
     * @brief Displays cell contours in a window.
     */
    void showMaskInWindow();

    /**
     * @brief Get the number of cells, one polygon obtained from mask is one cell
     */
    unsigned int getCellNum() const;

    /**
     * Get polygon information of all cells
     * @return A vector of polygons
     */
    const vector<Polygon> &getPolygons() const;

    /**
     *
     * @param border_array A pointer to a memory block, size = cell_num * 16 * 2 *sizeof(char)
     */
    void getBorders(char * border_array);

    void getEffectiveRectangle(unsigned int* effective_rect) const;

    unsigned int getBlockNum() const;

    unsigned int *getBlockIndex();

    const unsigned int *getBlockSize() const;

    void preBlockSort();

    static bool polygonComp(Polygon& p1, Polygon& p2);

    virtual ~Mask();
};

#endif //GEFTOOLS__MASK_H_
