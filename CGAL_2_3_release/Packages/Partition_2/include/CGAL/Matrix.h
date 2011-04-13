// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/Matrix.h
// package       : $CGAL_Package: Partition_2 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: 2D Matrix
// ============================================================================

#ifndef   CGAL_MATRIX_H
#define   CGAL_MATRIX_H

#include <vector>
#include <iostream>
#include <cstddef>

namespace CGAL {

template <class T>
class  Matrix : public std::vector< std::vector<T> > 
{
public:
   Matrix(size_t x = 0, size_t y = 0) : 
      std::vector< std::vector<T> > (x, std::vector<T>(y)), 
      _rows(x), _columns(y) 
   {}

   size_t rows() const { return _rows; }
   size_t columns() const { return _columns; }


protected:
   size_t _rows;
   size_t _columns;
};

template <class T> 
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
{
   typedef typename Matrix<T>::size_type size_type;

   for (size_type i = 0; i < m.rows(); i++) 
   {
      os << std::endl << i << " : ";
      for (size_type j = 0; j < m.columns(); j++) 
      {
         os << m[i][j] << " ";
      }
   }
   return os;
}

}

#endif // CGAL_MATRIX_H
