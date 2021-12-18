//
// Created by 黄志博 on 2021/12/14.
//

#ifndef GEFTOOLS__MASK_H_
#define GEFTOOLS__MASK_H_
#include <vector>
#include <string>
#include "polygen.h"

using namespace std;

class Mask {
  private:
    vector<Polygen> polygens_;
    unsigned int cell_num_;

  public:
    explicit Mask(const string& file);

    void to_mask();
    void write_polygens();
    unsigned int getCellNum() const;
    const vector<Polygen> &getPolygens() const;
};

#endif //GEFTOOLS__MASK_H_
