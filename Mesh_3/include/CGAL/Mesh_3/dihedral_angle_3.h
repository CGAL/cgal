// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DIHEDRAL_ANGLE_3_H
#define CGAL_DIHEDRAL_ANGLE_3_H

#include <CGAL/Kernel_traits.h>
#include <cmath>


namespace CGAL {
  
template <typename Tetrahedron_3>
double
min_dihedral_angle(const Tetrahedron_3& t)
{  
  return min_dihedral_angle(t, typename Kernel_traits<Tetrahedron_3>::Kernel());
}

template <typename Tetrahedron_3, typename K>
double
min_dihedral_angle(const Tetrahedron_3& t, K k=K())
{  
  typedef typename K::Point_3 Point_3;
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  
  const Point_3& a = vertex(t, 0);
  const Point_3& b = vertex(t, 1);
  const Point_3& c = vertex(t, 2);
  const Point_3& d = vertex(t, 3);
  
  double min_angle = dihedral_angle(a,b,c,d,k);
  min_angle = (std::min)(dihedral_angle(b,a,c,d,k), min_angle);
  min_angle = (std::min)(dihedral_angle(c,a,b,d,k), min_angle);
  min_angle = (std::min)(dihedral_angle(d,a,b,c,k), min_angle);
  
  return min_angle;
}  

template <typename Tetrahedron_3>
double
dihedral_angle(const Tetrahedron_3& t)
{
  return dihedral_angle(t, typename Kernel_traits<Tetrahedron_3>::Kernel());
}
  

template <typename Tetrahedron_3>
double
dihedral_angle(const Tetrahedron_3& t)
{
  return dihedral_angle(t, typename Kernel_traits<Tetrahedron_3>::Kernel());
}

template <typename Tetrahedron_3,
          typename K>
double
dihedral_angle(const Tetrahedron_3& t, K k = K())
{
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();

  return dihedral_angle(vertex(t, 0),
                        vertex(t, 1),
                        vertex(t, 2),
                        vertex(t, 3),
                        k);
}

template <typename Point_3>
double
dihedral_angle(const Point_3& a, const Point_3& b,
               const Point_3& c, const Point_3& d)
{
  return dihedral_angle(a, b, c, d, typename Kernel_traits<Point_3>::Kernel());
}

template <typename Point_3,
          typename K>
double
dihedral_angle(const Point_3& a, const Point_3& b,
               const Point_3& c, const Point_3& d, K k = K())
{
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();
  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();

  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  const Vector_3 ab = vector(a,b);
  const Vector_3 ac = vector(a,c);
  const Vector_3 ad = vector(a,d);

  const Vector_3 acad = cross_product(ac, ad);
  const double x = scalar_product(cross_product(ac, ab), acad);
  const FT l_ac = CGAL::sqrt(sq_distance(a, c));
  const double y = l_ac * scalar_product(ab, acad);

  return std::atan2(y, x);
}

} // end namespace CGAL

#endif // CGAL_DIHEDRAL_ANGLE_3_H
