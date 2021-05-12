// Copyright (c) 2005,2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_DIMENSION_UTILS_H
#define CGAL_KERNEL_DIMENSION_UTILS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>

namespace CGAL {

namespace Access {

// Some tools to find the type of a kernel object given its dimension.
// E.g. : Vector<K, 2>::type  is  K::Vector_2.
//
// Currently undocumented => for CGAL internal use only.

// TODO :
// - What about functors ?
//   At least those with a dimensional-independent interface.
// - Another thing which would be nice would be to make d-dimensional
//   algorithms use the 2d-3d kernel interface with a smaller wrapper.
//   (again, this forces a full kernel, not a sub-set traits, but at least...)
//   Then, based on Dimension<>, it's possible to automatically use it.

template < typename K, typename Dim_tag = typename K::Dimension>  struct Point;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Vector;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Iso_box;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Direction;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Line;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Ray;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Segment;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Triangle;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Hypersphere;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Hyperplane;
template < typename K, typename Dim_tag = typename K::Dimension>  struct Aff_transformation;

// Not-so generalizable stuff : Conic_2, Tetrahedron_3.
template < typename K, typename Dim_tag = typename K::Dimension>  struct Tetrahedron;


template < typename K >
struct Point <K, Dimension_tag<2> > { typedef typename K::Point_2 type; };

template < typename K >
struct Point <K, Dimension_tag<3> > { typedef typename K::Point_3 type; };

template < typename K >
struct Point <K, Dynamic_dimension_tag> { typedef typename K::Point_d type; };


template < typename K >
struct Vector <K, Dimension_tag<2> > { typedef typename K::Vector_2 type; };

template < typename K >
struct Vector <K, Dimension_tag<3> > { typedef typename K::Vector_3 type; };

template < typename K >
struct Vector <K, Dynamic_dimension_tag> { typedef typename K::Vector_d type; };


template < typename K >
struct Iso_box <K, Dimension_tag<2> > { typedef typename K::Iso_rectangle_2 type; };

template < typename K >
struct Iso_box <K, Dimension_tag<3> > { typedef typename K::Iso_cuboid_3 type; };

template < typename K >
struct Iso_box <K, Dynamic_dimension_tag> { typedef typename K::Iso_box_d type; };


template < typename K >
struct Direction <K, Dimension_tag<2> > { typedef typename K::Direction_2 type; };

template < typename K >
struct Direction <K, Dimension_tag<3> > { typedef typename K::Direction_3 type; };

template < typename K >
struct Direction <K, Dynamic_dimension_tag> { typedef typename K::Direction_d type; };


template < typename K >
struct Line <K, Dimension_tag<2> > { typedef typename K::Line_2 type; };

template < typename K >
struct Line <K, Dimension_tag<3> > { typedef typename K::Line_3 type; };

template < typename K >
struct Line <K, Dynamic_dimension_tag> { typedef typename K::Line_d type; };


template < typename K >
struct Ray <K, Dimension_tag<2> > { typedef typename K::Ray_2 type; };

template < typename K >
struct Ray <K, Dimension_tag<3> > { typedef typename K::Ray_3 type; };

template < typename K >
struct Ray <K, Dynamic_dimension_tag> { typedef typename K::Ray_d type; };


template < typename K >
struct Segment <K, Dimension_tag<2> > { typedef typename K::Segment_2 type; };

template < typename K >
struct Segment <K, Dimension_tag<3> > { typedef typename K::Segment_3 type; };

template < typename K >
struct Segment <K, Dynamic_dimension_tag> { typedef typename K::Segment_d type; };


template < typename K >
struct Triangle <K, Dimension_tag<2> > { typedef typename K::Triangle_2 type; };

template < typename K >
struct Triangle <K, Dimension_tag<3> > { typedef typename K::Triangle_3 type; };

template < typename K >
struct Triangle <K, Dynamic_dimension_tag> { typedef typename K::Triangle_d type; };


template < typename K >
struct Tetrahedron <K, Dimension_tag<3> > { typedef typename K::Tetrahedron_3 type; };

template < typename K >
struct Tetrahedron <K, Dynamic_dimension_tag> { typedef typename K::Tetrahedron_d type; };


template < typename K >
struct Hypersphere <K, Dimension_tag<2> > { typedef typename K::Circle_2 type; };

template < typename K >
struct Hypersphere <K, Dimension_tag<3> > { typedef typename K::Sphere_3 type; };

template < typename K >
struct Hypersphere <K, Dynamic_dimension_tag> { typedef typename K::Sphere_d type; };


template < typename K >
struct Hyperplane <K, Dimension_tag<2> > { typedef typename K::Line_2 type; };

template < typename K >
struct Hyperplane <K, Dimension_tag<3> > { typedef typename K::Plane_3 type; };

template < typename K >
struct Hyperplane <K, Dynamic_dimension_tag> { typedef typename K::Hyperplane_d type; };


template < typename K >
struct Aff_transformation <K, Dimension_tag<2> >
{ typedef typename K::Aff_transformation_2 type; };

template < typename K >
struct Aff_transformation <K, Dimension_tag<3> >
{ typedef typename K::Aff_transformation_3 type; };

template < typename K >
struct Aff_transformation <K, Dynamic_dimension_tag>
{ typedef typename K::Aff_transformation_d type; };

} // namespace Access

} //namespace CGAL

#endif // CGAL_KERNEL_DIMENSION_UTILS_H
