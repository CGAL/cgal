// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef   CGAL_MATRIX_H
#define   CGAL_MATRIX_H

#include <CGAL/license/Partition_2.h>


#include <vector>
#include <iostream>
#include <cstddef>

namespace CGAL {

template <class T>
class  Matrix : public std::vector< std::vector<T> > 
{
public:
   Matrix(std::size_t x = 0, std::size_t y = 0) : 
      std::vector< std::vector<T> > (x, std::vector<T>(y)), 
      _rows(x), _columns(y) 
   {}

   std::size_t rows() const { return _rows; }
   std::size_t columns() const { return _columns; }


protected:
   std::size_t _rows;
   std::size_t _columns;
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
