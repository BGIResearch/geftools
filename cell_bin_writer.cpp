//
// Created by huangzhibo on 2021/12/16.
//

#include "cell_bin_writer.h"

int cellBinWriter(const string &bin_gef, const string &mask_file, const string &outfile) {
    CommonBin common_bin_gef = CommonBin(bin_gef, 1);
    map<unsigned long long int, vector<unsigned int>> bin_gene_exp_map = common_bin_gef.getBinGeneExpMap();

    CellBin cell_bin_gef = CellBin(outfile, "w");

    Mask mask = Mask(mask_file);

    const vector<Polygon>& polygons = mask.getPolygons();

    for(unsigned int i = 0; i < mask.getCellNum(); i++){
        Polygon p = polygons[i];
        Rect roi = Rect(p.getMinX(),p.getMinY(),p.getCols(),p.getRows());
        Mat roi_mat = common_bin_gef.getWholeExpMatrix(roi);
        Mat fill_points = p.getFillPolyMat();
        roi_mat = roi_mat.mul(fill_points);

        vector<Point> non_zero_coordinates, non_zero_coordinates_offset;
        findNonZero(roi_mat,non_zero_coordinates);
        Point offset = Point(-p.getMinX(), -p.getMinY());
        offsetCoordinates(non_zero_coordinates, non_zero_coordinates_offset, offset);

        cell_bin_gef.addDnbInCell(
                non_zero_coordinates_offset,
                bin_gene_exp_map,
                p.getCenter(),
                p.getAreaUshort());
    }

    char* borders = static_cast<char *>(malloc(mask.getCellNum() * 16 * 2 * sizeof(char)));
    mask.getBorders(borders);

    ExpressionAttr expression_attr = common_bin_gef.getExpressionAttr();
    CellBinAttr cell_bin_attr = {
            .version = 1,
            .resolution = expression_attr.resolution,
            .offsetX = expression_attr.min_x,
            .offsetY = expression_attr.min_y
    };

    cell_bin_gef.storeAttr(cell_bin_attr);
    cell_bin_gef.storeCellBorder(borders, mask.getCellNum());
    cell_bin_gef.storeCell();
    cell_bin_gef.storeCellExp();


    return 0;
}


bool check_shape(){
    return true;
}

