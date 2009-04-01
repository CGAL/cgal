#ifndef CGAL_DSRPDB_MATRIX_H
#define CGAL_DSRPDB_MATRIX_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d.h>

namespace CGAL { namespace PDB {
//! Use TNT::Array2D<double> as a matrix
typedef TNT::Array2D<double> Matrix;

double det(const Matrix& m) {
  CGAL_assertion(m.dim1() == 3);
  CGAL_assertion(m.dim2() == 3);
  
  return (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
	  m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
	  m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]));
}


Matrix transpose(const Matrix& m) {
  Matrix mt(m.dim2(), m.dim1());
  for (int i = 0; i < m.dim1(); i++) {
    for (int j = 0; j < m.dim2(); j++) {
      mt[j][i] = m[i][j];
    }
  }
  return mt;
}



}}
#endif
