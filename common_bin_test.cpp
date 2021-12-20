//
// Created by huangzhibo on 2021/12/16.
//

#include "common_bin_test.h"
#include "cell_bin_writer.h"
#include <opencv2/opencv.hpp>

using namespace cv;

void test_getDnbStatMatrix() {
    cv::Mat img = cv::imread("/Users/huangzhibo/workitems/01.github/opencv_learning/test_data/mask.tif",-1);


    img = img.t();
    cout << "rows : " << img.rows << endl;
    cout << "cols : " << img.cols << endl;

    CommonBin common_bin_gef = CommonBin("/Users/huangzhibo/workitems/01.github/gefpy/test_data/barCode_gef/stereomics.h5", 1);
    cout << "/Users/huangzhibo/workitems/01.github/gefpy/test_data/barCode_gef/stereomics.h5" << endl;
    cout << "cell num: " << common_bin_gef.getGeneNum() << endl;
    cout << "x: " << common_bin_gef.getDnbStatMatrixShape()[0] << endl;
    cout << "y: " << common_bin_gef.getDnbStatMatrixShape()[1] << endl;


    common_bin_gef.getExpression();

    unsigned int x0 = 34, x1 = 42, y0 = 56, y1 = 58;
//    unsigned int x0 = 52, x1 = 56, y0 = 39, y1 = 42;
//    unsigned int x0 = 39, x1 = 42, y0 = 52, y1 = 56;
    unsigned int rows = x1 - x0 + 1, cols = y1 - y0 + 1;
//    unsigned char m[rows][cols];
//    Mat valid_dnb_mat = Mat::zeros(232, 223, CV_8UC1);
//    unsigned char * m;
//    m = static_cast<unsigned char *>(malloc((x1 - x0 + 1) * (y1 - y0 + 1) * sizeof(unsigned char)));
//
//    common_bin_gef.getDnbStatMatrix(x0, y0, rows, cols, valid_dnb_mat.data);
//    common_bin_gef.getDnbStatMatrix(0, 0, 232, 223, valid_dnb_mat.data);



//    cv::imwrite("orig.tif",valid_dnb_mat);
//    cv::waitKey(0);




//    Mat mat = Mat(rows, cols, CV_8UC1, m[0]);
//
//    mat = mat.mul(valid_dnb_mat。);

//    cout << valid_dnb_mat << endl;

//    vector<Point> pts;
//    pts.emplace_back(Point(3, 2));
//    pts.emplace_back(Point(5, 6));
//
//    Mat test_mat = Mat(pts);
//
//    cout << test_mat << endl;

//
//    cout << "m: " << m[0][0] << endl;
//    printf("%d ", m[0][0] );

//    string mask = "/Users/huangzhibo/workitems/01.github/gefpy/test_data/FP200000617TL_B6_mask.tif";
//    string bin_gef = "/Users/huangzhibo/workitems/01.github/gefpy/test_data/barcode_gene_exp_gef/stereomics.h5";
//    string cell_gef = "/Users/huangzhibo/workitems/01.github/geftools/test_data/output_test8.gef";

//    string mask = "/Users/huangzhibo/workitems/01.github/gefpy/test_data/T324_T307/FP200000341BR_A3_T324/FP200000341BR_A3_T324_mask.tif";
//    string bin_gef = "/Users/huangzhibo/workitems/01.github/gefpy/test_data/T324_T307/FP200000341BR_A3_T324/FP200000341BR_A3.bin_gef.h5";
//    string cell_gef = "/Users/huangzhibo/workitems/01.github/geftools/test_data/FP200000341BR_A3.test.gef";

    string bin_gef = "/Users/huangzhibo/workitems/01.github/gefpy/test_data/barCode_gef/stereomics.h5";
    string mask = "/Users/huangzhibo/workitems/01.github/opencv_learning/test_data/mask.tif";
    string cell_gef = "/Users/huangzhibo/workitems/01.github/geftools/test_data/demo2.h5";

    cellBinWriter(bin_gef, mask, cell_gef);

}
