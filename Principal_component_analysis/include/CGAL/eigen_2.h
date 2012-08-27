// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s) : Pierre Alliez

#ifndef CGAL_EIGEN_2_H
#define CGAL_EIGEN_2_H

namespace CGAL {

namespace internal {

  // extract eigenvalues and eigenvectors from a 2x2 symmetric
  // positive definite matrix. 
  // Note: computations involve a square root.
  // Matrix numbering:
  // a b
  // b c 
  // Eigen values and vectors are sorted in descendent order.
  template <typename K>
  void eigen_symmetric_2(const typename K::FT *matrix, // a b c
                         std::pair<typename K::Vector_2,
                                   typename K::Vector_2>& eigen_vectors,
                         std::pair<typename K::FT,
                                   typename K::FT>& eigen_values) 
  {
    // types
    typedef typename K::FT FT;
    typedef typename K::Vector_2 Vector;

    // for better reading
    FT a = matrix[0];
    FT b = matrix[1];
    FT c = matrix[2];
    FT p = c*c - 2*a*c + 4*b*b + a*a;
    CGAL_assertion(a >= 0.0 && c >= 0.0);

    // degenerate or isotropic case
    if(p == 0.0) 
    {
      // unit eigen values by default
      eigen_values.first = eigen_values.second = (FT)1.0;

      // any vector is eigen vector
      // the 2D canonical frame is output by default
      eigen_vectors.first  = Vector((FT)1.0,(FT)0.0);
      eigen_vectors.second = Vector((FT)0.0,(FT)1.0);
    }
    else 
    {
      if(b == 0.0) 
      {
	if(a>=c) 
	{
	  eigen_values.first  = a;
	  eigen_values.second = c;
	  eigen_vectors.first  = Vector((FT)1.0, (FT)0.0);
	  eigen_vectors.second  = Vector((FT)0.0, (FT)1.0);
	}
	else
	{
	  eigen_values.first  = c;
	  eigen_values.second = a;
	  eigen_vectors.first  = Vector((FT)0.0, (FT)1.0);
	  eigen_vectors.second  = Vector((FT)1.0, (FT)0.0);
	}
      }
      else // generic case
      {
        FT l1 = (FT)(0.5 * ( -1*std::sqrt(p) + c + a));
        FT l2 = (FT)(0.5 * (    std::sqrt(p) + c + a));

        // all eigen values of a symmetric positive
        // definite matrix must be real and positive
        // we saturate the values if this is not the 
        // case for floating point computations.
        l1 = (l1 < (FT)0.0) ? (FT)0.0 : l1;
        l2 = (l2 < (FT)0.0) ? (FT)0.0 : l2;

        // sort eigen values and vectors in descendent order.
        if(l1 >= l2)
        {
          eigen_values.first  = l1;
          eigen_values.second = l2;
          eigen_vectors.first  = Vector((FT)1.0, (FT)(-(std::sqrt(p)-c+a) / (2*b)));
          eigen_vectors.second = Vector((FT)1.0, (FT)( (std::sqrt(p)+c-a) / (2*b)));
        }
        else
        {
          eigen_values.first  = l2;
          eigen_values.second = l1;
          eigen_vectors.first  = Vector((FT)1.0, (FT)( (std::sqrt(p)+c-a) / (2*b)));
          eigen_vectors.second = Vector((FT)1.0, (FT)(-(std::sqrt(p)-c+a) / (2*b)));
        }
      } // end generic case
    } // end non-degenerate case
  } // end eigen_symmetric_2
} // end namespace internal

} //namespace CGAL

#endif // CGAL_EIGEN_2_H
