#include <iostream>
#include <fstream>
#include "utils.h"

bool copyFile(const string& src_file, const string& dst_file) {
    ifstream fin(src_file, ios::binary);
    ofstream fout(dst_file, ios::binary);

    bool ret = true;

    while(!fin.eof()){
        char buf;
        fin.read(&buf, sizeof(char));

        if(fin.eof()) break;

        if (fout.bad())
        {
            ret = false;
            break;
        }
        fout.write(&buf, sizeof(char));
    }

    fout.close();
    fin.close();
    return ret;
}
