cmake_minimum_required(VERSION 3.17)
project(geftools)

set(CMAKE_CXX_STANDARD 11)

include_directories("/Users/huangzhibo/software/app/hdf5-1.12.1/include" "/Users/huangzhibo/software/app/opencv-4.5.4/include/opencv4")
link_directories("/Users/huangzhibo/software/app/hdf5-1.12.1/lib" "/Users/huangzhibo/software/app/opencv-4.5.4/lib")

add_executable(geftools main.cpp polygon.cpp mask.cpp main_convert.cpp cell_bin.cpp common_bin.cpp common_bin.h main_view.cpp main_view.h cell_bin_writer.cpp cell_bin_writer.h common_bin_test.cpp common_bin_test.h utils.cpp utils.h)
target_link_libraries(geftools hdf5 opencv_core opencv_imgproc opencv_imgcodecs opencv_highgui opencv_video opencv_videoio)
