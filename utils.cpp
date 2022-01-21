#include "utils.h"

S32 getStrfTime()
{
    time_t timep;
    time (&timep);
    S32 tmp;
    strftime(tmp.value, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    return tmp;
};

time_t printTime(time_t prev, string message){
    time_t cur;
    time(&cur);
    std::cout << std::setw (30) << message;
    double cost;
    cost = difftime(cur, prev);
    printf(" - %.f sec\n",cost);
    return cur;
}

unsigned long printCpuTime(unsigned long prev, string message){
    unsigned long cur=clock();
    std::cout << std::setw (30) << message;
    printf (" - %.6f cpu sec\n", static_cast<double>(cur - prev) / CLOCKS_PER_SEC);
    prev=cur;
    return prev;
}

vector<string> split(const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.emplace_back(item);
    }

    return result;
}

vector<string> readLines(const string &filename){
    vector<string> result;
    char data[1000] = {0};
    ifstream infile;
    infile.open(filename);

    while(infile.getline(data,1000)){
        result.emplace_back(data);
    }

    if(!infile.eof()) {
        cerr << "Error to read file : " << filename << endl;
        exit(2);
    }

    infile.close();
    return result;
}

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

void offsetCoordinates(vector<Point> & coordinates, vector<Point> & new_coordinates, Point & offset_point){
    for (auto & coordinate: coordinates) {
        Point p = Point(coordinate.x - offset_point.x, coordinate.y - offset_point.y);
        new_coordinates.emplace_back(p);
    }
}


