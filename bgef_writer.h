/** @file bgef_writer.h
    @brief Declare a BgefWriter class for writing common bin gef.

    Created by huangzhibo on 2021/01/24.
*/

#ifndef GEFTOOLS_CGEF_WRITER_H
#define GEFTOOLS_CGEF_WRITER_H

#include "hdf5.h"
#include "utils.h"
#include "gef.h"

static constexpr unsigned int version = 2;

class BgefWriter {
  private:
    int binsize_;
    hid_t str32_type_;
    hid_t file_id_;
    hid_t gene_exp_group_id_;
    hid_t gene_exp_bin_group_id_;
    hid_t whole_exp_group_id_;

    unsigned int resolution_;
    bool verbose_ = false;

  public:
    BgefWriter(const string& output_filename, bool verbose);
    ~BgefWriter();

    bool createGMGroup(int bin);

    bool storeGene(vector<Expression> &exps, vector<Gene> &genes, DnbAttr &dnbAttr, unsigned int maxexp, int binsize);
    bool storeDnb(DnbMatrix & dnb_matrix, int binsize);
    bool storeStat(vector<GeneStat>& geneStat) const;

    unsigned int getResolution() const;

    void setResolution(unsigned int resolution);

};

#endif