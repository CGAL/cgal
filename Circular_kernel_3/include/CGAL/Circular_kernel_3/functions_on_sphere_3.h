// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_FUNCTIONS_ON_SPHERE_3_H
#define CGAL_SPHERICAL_KERNEL_FUNCTIONS_ON_SPHERE_3_H

namespace CGAL {

  // At the moment we dont need those functions
  // But in the future maybe (some make_x_monotone? etc..)
  template <class SK>
  typename SK::Circular_arc_point_3
  x_extremal_point(const Sphere_3<SK> & c, bool i)
  {
    return SphericalFunctors::x_extremal_point<SK>(c,i);
  }

  template <class SK, class OutputIterator>
  OutputIterator
  x_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
  {
    return SphericalFunctors::x_extremal_point<SK>(c,res);
  }

  template <class SK>
  typename SK::Circular_arc_point_3
  y_extremal_point(const Sphere_3<SK> & c, bool i)
  {
    return SphericalFunctors::y_extremal_point<SK>(c,i);
  }

  template <class SK, class OutputIterator>
  OutputIterator
  y_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
  {
    return SphericalFunctors::y_extremal_point<SK>(c,res);
  }

  template <class SK>
  typename SK::Circular_arc_point_3
  z_extremal_point(const Sphere_3<SK> & c, bool i)
  {
    return SphericalFunctors::z_extremal_point<SK>(c,i);
  }

  template <class SK, class OutputIterator>
  OutputIterator
  z_extremal_points(const Sphere_3<SK> & c, OutputIterator res)
  {
    return SphericalFunctors::z_extremal_point<SK>(c,res);
  }

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_FUNCTIONS_ON_SPHERE_3_H
