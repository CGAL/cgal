// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_TC_UTILITIES_H
#define CGAL_TC_UTILITIES_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>

#include <vector>

namespace CGAL {
namespace Tangential_complex_ {

  template <typename K>
  std::vector<typename K::Vector_d>
  compute_gram_schmidt_basis(
    std::vector<typename K::Vector_d> const& input_basis,
    K const& kernel)
  {
    typedef typename K::Vector_d Vector;
    typedef std::vector<Vector> Basis;
    const int D = Ambient_dimension<Vector>::value;
    
    // Kernel functors
    K::Squared_length_d        sqlen      = kernel.squared_length_d_object();
    K::Scaled_vector_d         scaled_vec = kernel.scaled_vector_d_object();
    //K::Scalar_product_d        inner_pdct = kernel.scalar_product_d_object();
    //K::Difference_of_vectors_d diff_vec   = kernel.difference_of_vectors_d_object();
    Get_functor<K, Scalar_product_tag>::type inner_pdct(kernel); // CJTODO TEMP
    Get_functor<K, Difference_of_vectors_tag>::type diff_vec(kernel);

    Basis output_basis;

    Basis::const_iterator inb_it = input_basis.begin();
    Basis::const_iterator inb_it_end = input_basis.end();
    for (int i = 0 ; inb_it != inb_it_end ; ++inb_it, ++i)
    {
      Vector u = *inb_it;

      Basis::iterator outb_it = output_basis.begin();
      for (int j = 0 ; j < i ; ++j)
      {
        Vector const& ej = *outb_it;
        Vector u_proj = scaled_vec(ej, inner_pdct(u, ej));
        u = diff_vec(u, u_proj);
      }

      output_basis.push_back(scaled_vec(u, 1./CGAL::sqrt(sqlen(u))));
    }

    return output_basis;
  }

} // namespace Tangential_complex_
} //namespace CGAL

#endif // CGAL_TC_UTILITIES_H
