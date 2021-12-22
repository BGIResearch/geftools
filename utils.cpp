#include "utils.h"

S32 getStrfTime()
{
    time_t timep;
    time (&timep);
    S32 tmp;
    strftime(tmp.value, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    return tmp;
};

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
