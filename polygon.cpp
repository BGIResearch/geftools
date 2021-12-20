#include "polygon.h"

bool Polygon::applyContour(const vector<Point>& contour){
    original_contour_size_ = static_cast<short>(contour.size());
    if(original_contour_size_ <= 2 ) return false;

    if(contour.size() > 16){
        double epsilon = 0.01 * arcLength(contour, true);
        approxPolyDP(contour, border_, epsilon, true);
    }else{
        border_ = contour;
    }

    Moments mu = moments(border_, true);

    border_size_ = static_cast<short>(border_.size());

    assert(border_size_ > 2);
    assert(mu.m00 > 0);

    center_ = Point(static_cast<int>(mu.m10/mu.m00), static_cast<int>(mu.m01/mu.m00));

    area_ = mu.m00;

    for(const auto& p : border_){
        min_x_ = p.x < min_x_ ? p.x : min_x_;
        max_x_ = p.x > max_x_ ? p.x : max_x_;
        min_y_ = p.y < min_y_ ? p.y : min_y_;
        max_y_ = p.y > max_y_ ? p.y : max_y_;
    }

    for(const auto& p : border_){
        Point relative_point = Point(p.x - min_x_, p.y - min_y_);
        relative_border_.emplace_back(relative_point);
    }

    cols_ = max_x_ - min_x_ + 1;
    rows_ = max_y_ - min_y_ + 1;

    return true;
}

const vector<Point> &Polygon::getBorder() const {
    return border_;
}

const Point &Polygon::getCenter() const {
    return center_;
}

double Polygon::getArea() const {
    return area_;
}

short Polygon::getBorderSize() const {
    return border_size_;
}

short Polygon::getOriginalContourSize() const {
    return original_contour_size_;
}

void Polygon::setMinMaxXY() {
    for(const auto& p : border_){
        min_x_ = p.x < min_x_ ? p.x : min_x_;
        max_x_ = p.x > max_x_ ? p.x : max_x_;
        min_y_ = p.y < min_y_ ? p.y : min_y_;
        max_y_ = p.y > max_y_ ? p.y : max_y_;
    }
    rows_ = max_x_ - min_x_ + 1;
    cols_ = max_y_ - min_y_ + 1;
}

int Polygon::getMinX() const {
    return min_x_;
}

int Polygon::getMaxX() const {
    return max_x_;
}

int Polygon::getMinY() const {
    return min_y_;
}

int Polygon::getMaxY() const {
    return max_y_;
}

int Polygon::getRows() const {
    return rows_;
}

int Polygon::getCols() const {
    return cols_;
}

Mat Polygon::getFillPolyMat() const {
    Mat fill_points = Mat::zeros(rows_, cols_, CV_8UC1);
    fillPoly(fill_points, relative_border_, 1);
    return fill_points;
}

const vector<Point> &Polygon::getRelativeBorder() const {
    return relative_border_;
}

unsigned short Polygon::getAreaUshort() const {
    return static_cast<unsigned short>(area_);
}
