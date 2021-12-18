#include <vector>
#include <opencv2/opencv.hpp> //Include file for every supported OpenCV function
#include "polygen.h"
using namespace cv;
using namespace std;

const int w = 500;
int levels = 3;

vector<vector<Point> > contours;
vector<Vec4i> hierarchy;

static void on_trackbar(int, void*)
{
    Mat cnt_img = Mat::zeros(w, w, CV_8UC3);
    int _levels = levels - 3;
    drawContours( cnt_img, contours, _levels <= 0 ? 3 : -1, Scalar(128,255,255),
                  3, LINE_AA, hierarchy, std::abs(_levels) );
    imshow("contours", cnt_img);
}

int main() {
    cv::Mat img = cv::imread("/Users/huangzhibo/workitems/13.github/opencv_learning/test_data/mask.tif",-1);
    if( img.empty() ) return -1;



    vector<vector<Point> > contours0;

    //findContours从二值图像中检索轮廓，并返回检测到的轮廓的个数
    findContours(
            img,
            contours0,
            hierarchy,
            RETR_EXTERNAL,
            CHAIN_APPROX_SIMPLE
    );

    contours.resize(contours0.size());
//    double epsilon;

    for( size_t k = 0; k < contours0.size(); k++ ){
//        epsilon = 0.01 * arcLength(contours0[k], true);
//        approxPolyDP(contours0[k], contours[k], epsilon, true);
          Polygen p = Polygen(contours0[k]);
//        Moments mu = moments(contours[k], true);
//        cout << contours[k] << endl;
        cout << p.border << endl;
        cout << p.border.size() << endl;
        cout << contours0[k].size() << endl;
        cout << p.center << endl;
        cout << p.area << endl;
//        cout << "m00 : " << mu.m00 << endl;
//        cout << "x : " << int(mu.m10/mu.m00) << ", y : " << int(mu.m01/mu.m00) << endl;
//        double area = contourArea(contours[k]);
//        cout << "area : " << area << endl;

    }

//    namedWindow( "contours", 1 );
//    createTrackbar( "levels+3", "contours", &levels, 7, on_trackbar );
//    on_trackbar(0,0);
//    waitKey();


//    cv::namedWindow( "Example1", cv::WINDOW_AUTOSIZE );
//    cv::imshow( "Example1", img );
//    cv::waitKey( 0 );
//    cv::destroyWindow( "Example1" );
    return 0;
}
