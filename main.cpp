#include <iostream>
#include "main_convert.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

using namespace cv;

static int usage()
{
    cerr << endl;
    cerr << "Program: geftools (Tools for manipulating GEFs)" << endl;
    cerr << "Version: " << PACKAGE_VERSION << endl;
//    cerr << "Contact: Huang Zhibo <huangzhibo@genomics.cn>\n" << endl;
    cerr << "Usage:   geftools <command> [options]\n" << endl;
    cerr << "Command: convert       Generate cell bin GEF according to common bin GEF and mask file" << endl;
    cerr << "         view          view GEF" << endl;
//    fprintf(stderr, "         h5ls          scan h5 file\n");
//    fprintf(stderr, "         mask          manipulating mask file\n");
    cerr << "\nNote: Please report issues at https://github.com/BGIResearch/geftools/issues" << endl;
    return 1;
}

int main(int argc, const char* argv[]){
    int ret;
    if (argc < 2) return usage();

    if (strcmp(argv[1], "convert") == 0) ret = convert(argc-1, const_cast<char **>(argv + 1));
    else {
        cerr << "[main] unrecognized command " << argv[1] << endl;
        return 1;
    }

    return ret;
}
