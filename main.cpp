#include <iostream>
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "main_convert.h"
#include "common_bin_test.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

using namespace cv;

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: geftools (Tools for manipulating GEFs)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
//  fprintf(stderr, "Contact: Huang Zhibo <huangzhibo@genomics.cn>\n\n");
  fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
  fprintf(stderr, "Command: convert       Generate cell bin GEF according to common bin GEF and mask file\n");
  fprintf(stderr, "         view          view GEF\n");
  fprintf(stderr, "         h5ls          scan h5 file\n");
  fprintf(stderr, "         mask          manipulating mask file\n");
  fprintf(stderr, "\n");
  fprintf(stderr,
          "Note: Please report issues at https://github.com/BGI-flexlab/bamqc/issues\n");
  return 1;
}

int main(int argc, const char* argv[])
{

    test_getDnbStatMatrix();
    return 0;



//  usage();
//  return convert(argc-1, const_cast<char **>(argv + 1));

}
