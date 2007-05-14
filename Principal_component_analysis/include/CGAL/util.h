// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://gankit@scm.gforge.inria.fr/svn/cgal/trunk/Principal_component_analysis/include/CGAL/linear_least_squares_fitting_triangles.h $
// $Id: linear_least_squares_fitting_2.h 37882 2007-04-03 15:15:30Z spion $
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/Linear_algebraCd.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// Initialize a matrix in n dimension by an array or numbers
template<typename K>
typename CGAL::Linear_algebraCd<typename K::FT>::Matrix
init_Matrix(int n, typename K::FT entries[]) {
  CGAL_assertion(n>1); // the dimension > 1
  typedef typename CGAL::Linear_algebraCd<typename K::FT>::Matrix Matrix;
  Matrix ret(n);
  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      ret[i][j]=entries[i*n+j];
    }
  }
  return ret;
} // end initialization of matrix

} // end namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H
