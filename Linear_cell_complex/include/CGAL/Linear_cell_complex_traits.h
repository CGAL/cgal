// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_TRAITS_H
#define CGAL_LINEAR_CELL_COMPLEX_TRAITS_H 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/Origin.h>

namespace CGAL {

  template <unsigned int d>
  struct LCC_default_kernel
  { typedef Cartesian_d<double> type; };
  template <>
  struct LCC_default_kernel<2>
  { typedef Exact_predicates_inexact_constructions_kernel type; };
  template <>
  struct LCC_default_kernel<3>
  { typedef Exact_predicates_inexact_constructions_kernel type; };  

  /** Trait class for Linear_cell_complex class.
   *  dD version (for the moment there is only one dD kernel in CGAL).
   */
  template <unsigned int d_,
            class Kernel=typename LCC_default_kernel<d_>::type >
  struct Linear_cell_complex_traits : public Kernel
  {
    static const unsigned int ambient_dimension = d_;
    
    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_d     Point;
    typedef typename Kernel::Vector_d    Vector;
    
    // Constructions
    struct Construct_translated_point
    {
      Point operator() (const Point&p, const Vector& v)
      { return p+v; }
      Point operator() (const CGAL::Origin&, const Vector& v)
      { return operator() (Point(ambient_dimension, CGAL::Origin()), v); }
    };

    struct Construct_vector
    {
      Vector operator() (const Point& p1, const Point& p2)
      { return p2-p1; }
      Vector operator() (const CGAL::Origin&, const Point& p)
      { return operator() (Point(ambient_dimension, CGAL::Origin()), p); }
    };
    
    struct Construct_sum_of_vectors
    {
      Vector operator() (const Vector&v1, const Vector& v2)
      { return v1+v2; }
    };

    struct Construct_scaled_vector
    {
      Vector operator() (const Vector& v, typename Kernel::FT scale)
      { return scale*v; }
    };

    struct Construct_midpoint
    {
      Point operator() (const Point&p1, const Point& p2)
      { return typename Kernel::Midpoint_d()(p1, p2); }
    };
  };

  /** Trait class for Linear_cell_complex class.
   *  2D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<2,Kernel> : public Kernel
  {
    static const unsigned int ambient_dimension = 2;
    
    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_2     Point;
    typedef typename Kernel::Vector_2    Vector;

    // Constructions
    typedef typename Kernel::Construct_translated_point_2
    Construct_translated_point;

    typedef typename Kernel::Construct_vector_2 Construct_vector;

    typedef typename Kernel::Construct_sum_of_vectors_2
    Construct_sum_of_vectors;
    
    typedef typename Kernel::Construct_scaled_vector_2
    Construct_scaled_vector;

    typedef typename Kernel::Construct_midpoint_2
    Construct_midpoint;
  };

  /** Trait class for Linear_cell_complex class.
   *  3D version specialization.
   */
  template <class Kernel>
  struct Linear_cell_complex_traits<3,Kernel> : public Kernel
  {
    static const unsigned int ambient_dimension = 3;

    typedef typename Kernel::FT          FT;
    typedef typename Kernel::Point_3     Point;
    typedef typename Kernel::Vector_3    Vector;
    
    // Constructions
    typedef typename Kernel::Construct_translated_point_3 
    Construct_translated_point;

    typedef typename Kernel::Construct_vector_3 Construct_vector;

    typedef typename Kernel::Construct_sum_of_vectors_3
    Construct_sum_of_vectors;

    typedef typename Kernel::Construct_scaled_vector_3 
    Construct_scaled_vector;
    
    typedef typename Kernel::Construct_midpoint_3
    Construct_midpoint;
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_TRAITS_H //
// EOF //
