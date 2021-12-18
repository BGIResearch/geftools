//
// Created by huangzhibo on 2021/12/14.
//

#include "mask.h"
Mask::Mask(const string& file){
    cv::Mat img = cv::imread(file,-1);
    if( img.empty() ) { exit(-1);}

    img = img.t();

    rows_ = img.rows;
    cols_ = img.cols;

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    //findContours从二值图像中检索轮廓，并返回检测到的轮廓的个数
    findContours(
      img,
      contours,
      hierarchy,
      RETR_EXTERNAL,
      CHAIN_APPROX_SIMPLE
    );

    for(auto & contour : contours){
        Polygon p = Polygon();
        bool cell_polygon_is_good = p.applyContour(contour);
        if(cell_polygon_is_good){
            polygons_.emplace_back(p);
        }
    }

    cell_num_ = polygons_.size();
}

void Mask::write_polygons() {

}

const vector<Polygon> &Mask::getPolygons() const {
    return polygons_;
}

unsigned int Mask::getCellNum() const {
    return cell_num_;
}
