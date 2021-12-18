//
// Created by 黄志博 on 2021/12/13.
//
#include "polygen.h"

bool Polygen::applyContour(const vector<Point>& contour){
    original_contour_size_ = static_cast<short>(contour.size());
    if(original_contour_size_ <= 2 ) return false;

    if(contour.size() > 16){
        double epsilon = 0.01 * arcLength(contour, true);
        approxPolyDP(contour, border_, epsilon, true);
    }else{
        border_ = contour;
    }

    setMinMaxXY();

    Moments mu = moments(border_, true);

    border_size_ = static_cast<short>(border_.size());

    assert(border_size_ > 2);
    assert(mu.m00 > 0);

    center_ = Point(static_cast<int>(mu.m10/mu.m00), static_cast<int>(mu.m01/mu.m00));
    area_ = mu.m00;
    return true;
}

const vector<Point> &Polygen::getBorder() const {
    return border_;
}

const Point &Polygen::getCenter() const {
    return center_;
}

double Polygen::getArea() const {
    return area_;
}

short Polygen::getBorderSize() const {
    return border_size_;
}

short Polygen::getOriginalContourSize() const {
    return original_contour_size_;
}

void Polygen::setMinMaxXY() {
    for(const auto& p : border_){
        min_x_ = p.x < min_x_ ? p.x : min_x_;
        max_x_ = p.x > max_x_ ? p.x : max_x_;
        min_y_ = p.x < min_y_ ? p.x : min_y_;
        max_y_ = p.x > max_y_ ? p.x : max_y_;
    }
    rows_ = max_x_ - min_x_ + 1;
    cols_ = max_y_ - min_y_ + 1;
}

int Polygen::getMinX() const {
    return min_x_;
}

int Polygen::getMaxX() const {
    return max_x_;
}

int Polygen::getMinY() const {
    return min_y_;
}

int Polygen::getMaxY() const {
    return max_y_;
}

int Polygen::getRows() const {
    return rows_;
}

int Polygen::getCols() const {
    return cols_;
}

Mat Polygen::getFillPolyMat() const {
    Mat fill_points = Mat::zeros(rows_, cols_, CV_8UC1);
    fillPoly(border_, fill_points, 1);
    return fill_points;
}
