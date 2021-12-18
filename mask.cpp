//
// Created by huangzhibo on 2021/12/14.
//

#include "mask.h"
Mask::Mask(const string& file){
    cv::Mat img = cv::imread(file,-1);
    if( img.empty() ) { exit(-1);}

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
        Polygen p = Polygen();
        bool cell_polygen_is_good = p.applyContour(contour);
        if(cell_polygen_is_good){
            polygens_.emplace_back(p);
        }
    }

    cell_num_ = polygens_.size();
}

void Mask::write_polygens() {

}

const vector<Polygen> &Mask::getPolygens() const {
    return polygens_;
}

unsigned int Mask::getCellNum() const {
    return cell_num_;
}
