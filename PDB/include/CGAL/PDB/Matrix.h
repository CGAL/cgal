#ifndef CGAL_DSRPDB_MATRIX_H
#define CGAL_DSRPDB_MATRIX_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d.h>

CGAL_PDB_BEGIN_NAMESPACE
  //! Use TNT::Array2D<double> as a matrix
  typedef TNT::Array2D<double> Matrix;
CGAL_PDB_END_NAMESPACE
#endif
