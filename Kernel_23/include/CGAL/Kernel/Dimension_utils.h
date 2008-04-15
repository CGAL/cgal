// Copyright (c) 2005  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_DIMENSION_UTILS_H
#define CGAL_KERNEL_DIMENSION_UTILS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>

CGAL_BEGIN_NAMESPACE

// Some tools to find the type of a kernel object given its dimension.
// E.g. : Vector<2, K>::type  is  K::Vector_2.
//
// Currently undocumented => for CGAL internal use only.

// TODO :
// - What about functors ?
//   At least those with a dimensional-independent interface.
// - Another thing which would be nice would be to make d-dimensional
//   algorithms use the 2d-3d kernel interface with a smaller wrapper.
//   (again, this forces a full kernel, not a sub-set traits, but at least...)
//   Then, based on Dimension<>, it's possibly to automatically use it.
//   cf Dimension_mapper<> at the end of the file.

template < typename Dim_tag, typename K >  struct Point;
template < typename Dim_tag, typename K >  struct Vector;
template < typename Dim_tag, typename K >  struct Iso_box;
template < typename Dim_tag, typename K >  struct Direction;
template < typename Dim_tag, typename K >  struct Line;
template < typename Dim_tag, typename K >  struct Ray;
template < typename Dim_tag, typename K >  struct Segment;
template < typename Dim_tag, typename K >  struct Triangle;
template < typename Dim_tag, typename K >  struct Hypersphere;
template < typename Dim_tag, typename K >  struct Hyperplane;
template < typename Dim_tag, typename K >  struct Aff_transformation;

// Not-so generalizable stuff : Conic_2, Tetrahedron_3.
template < typename Dim_tag, typename K >  struct Tetrahedron;


template < typename K >
struct Point <Dimension_tag<2>, K> { typedef typename K::Point_2 type; };

template < typename K >
struct Point <Dimension_tag<3>, K> { typedef typename K::Point_3 type; };

template < typename K >
struct Point <Dynamic_dimension_tag, K> { typedef typename K::Point_d type; };


template < typename K >
struct Vector <Dimension_tag<2>, K> { typedef typename K::Vector_2 type; };

template < typename K >
struct Vector <Dimension_tag<3>, K> { typedef typename K::Vector_3 type; };

template < typename K >
struct Vector <Dynamic_dimension_tag, K> { typedef typename K::Vector_d type; };


template < typename K >
struct Iso_box <Dimension_tag<2>, K> { typedef typename K::Iso_rectangle_2 type; };

template < typename K >
struct Iso_box <Dimension_tag<3>, K> { typedef typename K::Iso_cuboid_3 type; };

template < typename K >
struct Iso_box <Dynamic_dimension_tag, K> { typedef typename K::Iso_box_d type; };


template < typename K >
struct Direction <Dimension_tag<2>, K> { typedef typename K::Direction_2 type; };

template < typename K >
struct Direction <Dimension_tag<3>, K> { typedef typename K::Direction_3 type; };

template < typename K >
struct Direction <Dynamic_dimension_tag, K> { typedef typename K::Direction_d type; };


template < typename K >
struct Line <Dimension_tag<2>, K> { typedef typename K::Line_2 type; };

template < typename K >
struct Line <Dimension_tag<3>, K> { typedef typename K::Line_3 type; };

template < typename K >
struct Line <Dynamic_dimension_tag, K> { typedef typename K::Line_d type; };


template < typename K >
struct Ray <Dimension_tag<2>, K> { typedef typename K::Ray_2 type; };

template < typename K >
struct Ray <Dimension_tag<3>, K> { typedef typename K::Ray_3 type; };

template < typename K >
struct Ray <Dynamic_dimension_tag, K> { typedef typename K::Ray_d type; };


template < typename K >
struct Segment <Dimension_tag<2>, K> { typedef typename K::Segment_2 type; };

template < typename K >
struct Segment <Dimension_tag<3>, K> { typedef typename K::Segment_3 type; };

template < typename K >
struct Segment <Dynamic_dimension_tag, K> { typedef typename K::Segment_d type; };


template < typename K >
struct Triangle <Dimension_tag<2>, K> { typedef typename K::Triangle_2 type; };

template < typename K >
struct Triangle <Dimension_tag<3>, K> { typedef typename K::Triangle_3 type; };

template < typename K >
struct Triangle <Dynamic_dimension_tag, K> { typedef typename K::Triangle_d type; };


template < typename K >
struct Tetrahedron <Dimension_tag<3>, K> { typedef typename K::Tetrahedron_3 type; };

template < typename K >
struct Tetrahedron <Dynamic_dimension_tag, K> { typedef typename K::Tetrahedron_d type; };


template < typename K >
struct Hypersphere <Dimension_tag<2>, K> { typedef typename K::Circle_2 type; };

template < typename K >
struct Hypersphere <Dimension_tag<3>, K> { typedef typename K::Sphere_3 type; };

template < typename K >
struct Hypersphere <Dynamic_dimension_tag, K> { typedef typename K::Sphere_d type; };


template < typename K >
struct Hyperplane <Dimension_tag<2>, K> { typedef typename K::Line_2 type; };

template < typename K >
struct Hyperplane <Dimension_tag<3>, K> { typedef typename K::Plane_3 type; };

template < typename K >
struct Hyperplane <Dynamic_dimension_tag, K> { typedef typename K::Hyperplane_d type; };


template < typename K >
struct Aff_transformation <Dimension_tag<2>, K>
{ typedef typename K::Aff_transformation_2 type; };

template < typename K >
struct Aff_transformation <Dimension_tag<3>, K>
{ typedef typename K::Aff_transformation_3 type; };

template < typename K >
struct Aff_transformation <Dynamic_dimension_tag, K>
{ typedef typename K::Aff_transformation_d type; };


#if 0
// K is a kernel, T one of its object.
// Dimension_mapper provides a Kernel_d like interface to 2d/3d/dd kernel API.

// NB : one issue is that it _requires_ a full kernel interface as parameter,
// so you can't provide a simple traits (sub-set of the kernel).
// Might be problematic...
// But it can also be useful.
// Should this be part of a grand scheme for better handling of d-dim in CGAL ?
template < typename T,
           typename K = typename Kernel_traits<T>::Kernel,
           typename Dim_tag = Dimension_of<T, K>::value >
class Dimension_mapper;

template < typename T, typename K >
struct Dimension_mapper <T, K, 2>
  : public K
{
  Dimension_mapper() {}
  Dimension_mapper(const K& k) : K(k) {}

  typedef typename K::Point_2    Point_d;
  typedef typename K::Vector_2   Vector_d;
  // ...

  typedef typename K::Less_x_2   Less_x_d;
  // ...

  Less_x_d less_x_d_object() const { return this->less_x_2_object(); }
  // ...
};

template < typename T, typename K >
struct Dimension_mapper <T, K, 3>
  : public K
{
  Dimension_mapper() {}
  Dimension_mapper(const K& k) : K(k) {}

  typedef typename K::Point_3    Point_d;
  typedef typename K::Vector_3   Vector_d;
  // ...
};

template < typename T, typename K >
struct Dimension_mapper <T, K, Dynamic_dimension_tag>
  : public K
{
  Dimension_mapper() {}
  Dimension_mapper(const K& k) : K(k) {}
};
#endif

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_DIMENSION_UTILS_H
