//
// Created by huangzhibo on 2021/12/16.
//

#include "cell_bin_writer.h"
#include "mask.h"
#include "common_bin.h"
#include "cell_bin.h"

using namespace cv;

int cellBinWriter(const string &bin_gef, const string &mask_file, const string &outfile) {
    CommonBin common_bin_gef = CommonBin(bin_gef, 1);
//    common_bin_gef.

    CellBin cell_bin_gef = CellBin(outfile, "w");

    Mask mask = Mask(mask_file);

    const vector<Polygen>& polygens = mask.getPolygens();

    for(unsigned int i = 0; i < mask.getCellNum(); i++){
        Polygen p = polygens[i];
        Mat fill_points = p.getFillPolyMat();
//        unsigned char validDnbArray[p.getRows()][p.getCols()];
        Mat valid_dnb_mat = Mat::zeros(p.getRows(), p.getCols(), CV_8UC1);
        common_bin_gef.getDnbStatMatrix(p.getMinX(), p.getMinY(), p.getRows(), p.getCols(), valid_dnb_mat.data);
        valid_dnb_mat = valid_dnb_mat.mul(fill_points);
        vector<Point> non_zero_coordinates;
        findNonZero(valid_dnb_mat,non_zero_coordinates);
    }

    return 0;
}


bool check_shape(){
    return true;
}

