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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_DIMENSION_UTILS_H
#define CGAL_KERNEL_DIMENSION_UTILS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/Dimension.h>

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

template < int dim, typename K >  struct Point;
template < int dim, typename K >  struct Vector;
template < int dim, typename K >  struct Iso_box;
template < int dim, typename K >  struct Direction;
template < int dim, typename K >  struct Line;
template < int dim, typename K >  struct Ray;
template < int dim, typename K >  struct Segment;
template < int dim, typename K >  struct Triangle;
template < int dim, typename K >  struct Hypersphere;
template < int dim, typename K >  struct Hyperplane;
template < int dim, typename K >  struct Aff_transformation;

// Not-so generalizable stuff : Conic_2, Tetrahedron_3.
template < int dim, typename K >  struct Tetrahedron;


template < typename K >
struct Point <2, K> { typedef typename K::Point_2 type; };

template < typename K >
struct Point <3, K> { typedef typename K::Point_3 type; };

template < typename K >
struct Point <0, K> { typedef typename K::Point_d type; };


template < typename K >
struct Vector <2, K> { typedef typename K::Vector_2 type; };

template < typename K >
struct Vector <3, K> { typedef typename K::Vector_3 type; };

template < typename K >
struct Vector <0, K> { typedef typename K::Vector_d type; };


template < typename K >
struct Iso_box <2, K> { typedef typename K::Iso_rectangle_2 type; };

template < typename K >
struct Iso_box <3, K> { typedef typename K::Iso_cuboid_3 type; };

template < typename K >
struct Iso_box <0, K> { typedef typename K::Iso_box_d type; };


template < typename K >
struct Direction <2, K> { typedef typename K::Direction_2 type; };

template < typename K >
struct Direction <3, K> { typedef typename K::Direction_3 type; };

template < typename K >
struct Direction <0, K> { typedef typename K::Direction_d type; };


template < typename K >
struct Line <2, K> { typedef typename K::Line_2 type; };

template < typename K >
struct Line <3, K> { typedef typename K::Line_3 type; };

template < typename K >
struct Line <0, K> { typedef typename K::Line_d type; };


template < typename K >
struct Ray <2, K> { typedef typename K::Ray_2 type; };

template < typename K >
struct Ray <3, K> { typedef typename K::Ray_3 type; };

template < typename K >
struct Ray <0, K> { typedef typename K::Ray_d type; };


template < typename K >
struct Segment <2, K> { typedef typename K::Segment_2 type; };

template < typename K >
struct Segment <3, K> { typedef typename K::Segment_3 type; };

template < typename K >
struct Segment <0, K> { typedef typename K::Segment_d type; };


template < typename K >
struct Triangle <2, K> { typedef typename K::Triangle_2 type; };

template < typename K >
struct Triangle <3, K> { typedef typename K::Triangle_3 type; };

template < typename K >
struct Triangle <0, K> { typedef typename K::Triangle_d type; };


template < typename K >
struct Tetrahedron <3, K> { typedef typename K::Tetrahedron_3 type; };

template < typename K >
struct Tetrahedron <0, K> { typedef typename K::Tetrahedron_d type; };


template < typename K >
struct Hypersphere <2, K> { typedef typename K::Circle_2 type; };

template < typename K >
struct Hypersphere <3, K> { typedef typename K::Sphere_3 type; };

template < typename K >
struct Hypersphere <0, K> { typedef typename K::Sphere_d type; };


template < typename K >
struct Hyperplane <2, K> { typedef typename K::Line_2 type; };

template < typename K >
struct Hyperplane <3, K> { typedef typename K::Plane_3 type; };

template < typename K >
struct Hyperplane <0, K> { typedef typename K::Hyperplane_d type; };


template < typename K >
struct Aff_transformation <2, K>
{ typedef typename K::Aff_transformation_2 type; };

template < typename K >
struct Aff_transformation <3, K>
{ typedef typename K::Aff_transformation_3 type; };

template < typename K >
struct Aff_transformation <0, K>
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
           int dim = Dimension_of<T, K>::value >
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
struct Dimension_mapper <T, K, 0>
  : public K
{
  Dimension_mapper() {}
  Dimension_mapper(const K& k) : K(k) {}
};
#endif

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_DIMENSION_UTILS_H
