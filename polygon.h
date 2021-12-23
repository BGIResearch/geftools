/** @file polygon.h
    @brief Declare a Polygon class to describe cell shape.

    Created by huangzhibo on 2021/12/13..
*/

#ifndef GEFTOOLS__POLYGEN_H
#define GEFTOOLS__POLYGEN_H
#include <vector>
#include "opencv2/opencv.hpp"

using namespace cv;
using namespace std;

/**
 * @brief A polygon class to describe cell shape.
 */
class Polygon {
  private:
    vector<Point> border_;
    vector<Point> relative_border_;
public:
    const vector<Point> &getRelativeBorder() const;

private:
    Point center_;
    double area_;
    short border_size_;
    short original_contour_size_;
    int min_x_{INT_MAX}, max_x_{0}, min_y_{INT_MAX}, max_y_{0}, rows_{0}, cols_{0};

    void setMinMaxXY();

  public:
    /**
     * @brief Get border points of this polygon
     * @return A vector of points, point num \<= 16
     */
    const vector<Point> &getBorder() const;

//    getBorderRelativeToCenter();

    /// Get border points number of this polygon
    short getBorderSize() const;

    /// Get original contour points number before apply approxPolyDP
    short getOriginalContourSize() const;

    /**
     * @brief Get a matrix with value 1 filled in this polygon region, outside of the region filled by 0
     * @return A matrix with type CV_8UC1, new Mat::zeros(rows_, cols_, CV_8UC1), and fill 1 into this polygon region
     */
    Mat getFillPolyMat() const;

    /// Get center point of this polygon
    const Point &getCenter() const;

    /// Get area of this polygon
    double getArea() const;

    /**
     * @brief Gets area of this polygon.
     *
     * Ceil up to the nearest (unsigned short) integer
     */

    unsigned short getAreaUshort() const;

    /// Get min X coordinate of this polygon.
    int getMinX() const;

    /// Get max X coordinate of this polygon.
    int getMaxX() const;

    /// Get min Y coordinate of this polygon.
    int getMinY() const;

    /// Get max Y coordinate of this polygon.
    int getMaxY() const;

    /// Get rows number, rows_ = max_x_ - min_x_ + 1
    int getRows() const;

    /// Get cols number, cols = max_y_ - min_y_ + 1
    int getCols() const;

    /**
     * @brief apply contour info from opencv findContours,
     * must call this function before use other functions in this class
     * @param contour A vector of points, one contour of opencv findContours contours
     * @return true if the cell polygon is proper (more then two contour points), else false
     */
    bool applyContour(const vector<Point>& contour);
};

#endif //GEFTOOLS__POLYGEN_H
