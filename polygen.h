//
// Created by huangzhibo on 2021/12/13.
//

#ifndef GEFTOOLS__POLYGEN_H
#define GEFTOOLS__POLYGEN_H
#include <vector>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

class Polygen {
  private:
    vector<Point> border_;
    Point center_;
    double area_;
    short border_size_;
    short original_contour_size_;
    int min_x_{0}, max_x_{0}, min_y_{0}, max_y_{0}, rows_{0}, cols_{0};

    void setMinMaxXY();

  public:
    const vector<Point> &getBorder() const;

    Mat getFillPolyMat() const;

    const Point &getCenter() const;

    double getArea() const;

    short getBorderSize() const;

    short getOriginalContourSize() const;

    int getMinX() const;

    int getMaxX() const;

    int getMinY() const;

    int getMaxY() const;

    int getRows() const;

    int getCols() const;

    bool applyContour(const vector<Point>& contour);
};

#endif //GEFTOOLS__POLYGEN_H
