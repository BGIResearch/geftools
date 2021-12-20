//
// Created by huangzhibo on 2021/12/14.
//

#include "mask.h"
Mask::Mask(const string& file){
    cv::Mat img = cv::imread(file,-1);
    if( img.empty() ) { exit(-1);}

    rows_ = img.rows;
    cols_ = img.cols;

    //findContours从二值图像中检索轮廓，并返回检测到的轮廓的个数
    findContours(
      img,
      contours_,
      hierarchy_,
      RETR_EXTERNAL,
      CHAIN_APPROX_SIMPLE
    );

    for(auto & contour : contours_){
        Polygon p = Polygon();
        bool cell_polygon_is_good = p.applyContour(contour);
        if(cell_polygon_is_good){
            polygons_.emplace_back(p);
        }
    }

    cell_num_ = polygons_.size();
}

const vector<Polygon> &Mask::getPolygons() const {
    return polygons_;
}

unsigned int Mask::getCellNum() const {
    return cell_num_;
}

void Mask::showMaskInWindow() {
    Mat cnt_img = Mat::zeros(rows_, cols_, CV_8UC3);
    drawContours( cnt_img, contours_, -1, Scalar(128,255,255),
                  3, LINE_AA, hierarchy_, 3 );
    imshow("Mask Contours", cnt_img);
    waitKey(0);
}

void Mask::getBorders(char * border_array) {
    for (unsigned int i = 0; i < cell_num_; i++) {
        Polygon polygon = polygons_[i];
        vector<Point> border = polygon.getBorder();
        const Point& center = polygon.getCenter();
        unsigned int index1 = i * 32;

        auto border_size = static_cast<short>(border.size());

        for (short j = 0; j < 16; ++j) {
            unsigned int index2 = index1 + (j << 1);
            if(j >= border_size){
                border_array[index2] = 0;
                border_array[index2+1] = 0;
            }else{
                Point p = border[j];
                border_array[index2] = static_cast<char>(p.x - center.x);
                border_array[index2+1] = static_cast<char>(p.y - center.y);
            }
        }
    }
}
