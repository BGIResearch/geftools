#include <iostream>
#include "main_bgef.h"
#include "main_cgef.h"
#include "main_view.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.5"
#endif

using namespace cv;

static int usage()
{
    cerr << endl;
    cerr << "Program: geftools (Tools for manipulating GEFs)" << endl;
    cerr << "Version: " << PACKAGE_VERSION << endl;
//    cerr << "Contact: Huang Zhibo <huangzhibo@genomics.cn>\n" << endl;
    cerr << "Usage:   geftools <command> [options]\n" << endl;
    cerr << "Command: bgef          Generate common bin GEF(.bgef) according to gem file or bin1 GEF" << endl;
    cerr << "         cgef          Generate cell bin GEF(.cgef) according to common bin GEF and mask file" << endl;
    cerr << "         view          View GEF, generate gem file" << endl;
//    fprintf(stderr, "         h5ls          scan h5 file\n");
//    fprintf(stderr, "         mask          manipulating mask file\n");
    cerr << "\nNote: Please report issues at https://github.com/BGIResearch/geftools/issues" << endl;
    return 1;
}

int main(int argc, const char* argv[]){
    time_t prev;
    time(&prev);

    int ret;
    if (argc < 2) return usage();

    if (strcmp(argv[1], "bgef") == 0) ret = bgef(argc-1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "cgef") == 0) ret = cgef(argc-1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "view") == 0) ret = view(argc-1, const_cast<char **>(argv + 1));
    else {
        cerr << "[main] unrecognized command " << argv[1] << endl;
        return 1;
    }

    printTime(prev, "Total Time");
    return ret;
}
