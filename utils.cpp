//
// Created by huangzhibo on 2021/12/17.
//

#include <iostream>
#include <fstream>
#include "utils.h"

bool copyFile(const string &inPath, const string &outPath) {
    ifstream fin(inPath, ios::binary);
    ofstream fout(outPath, ios::binary);

    bool bRet = true;

    while(!fin.eof()){
        char szBuf;
        fin.read(&szBuf, sizeof(char));

        if(fin.eof()) break;

        if (fout.bad())
        {
            bRet = false;
            break;
        }
        fout.write(&szBuf, sizeof(char));
    }

    fout.close();
    fin.close();
    return bRet;
}
