#ifndef DSRPDB_ALIGN_H
#define DSRPDB_ALIGN_H
#include <dsrpdb/Protein.h>
#include <dsrpdb/Transform.h>
#include <vector>

namespace dsrpdb {
  /*!  Compute an optimal rigid alignment of the second protein to the
    first and transform it accordingly. The two proteins must have the
    same number of residues and all C_alpha atom possitions must be
    defined.

    \example pdb_align.cc
  */
  void align_second_protein_to_first(const Protein &base, Protein &o);

  /*!
    Compute the transformation matrix which maps the second point set 
    to the first. 

    \example pdb_align_points.cc
  */
  Transform compute_transform_taking_first_to_second(const std::vector<Point> &a,
						     const std::vector<Point> &b);
};
#endif
