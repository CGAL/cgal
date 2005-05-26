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

#ifndef CGAL_KERNEL_DIMENSION_H
#define CGAL_KERNEL_DIMENSION_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>

CGAL_BEGIN_NAMESPACE

// This is a tool to obtain the static dimension of a kernel object.
// Currently undocumented => for CGAL internal use only.

// TODO :
// - Should the default be that kernel objects export this value themselves ?
//   This way, user defined objects will work as well ?
//   But an object could have a dimension depending on the traits...
//   => Well, let's make it work first.
// - Should it be named Dimension_traits<> ?

template < typename T, typename K = typename Kernel_traits<T>::Kernel >
struct Dimension;

// 2D objects

template < typename K >
struct Dimension < typename K::Point_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Vector_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Direction_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Line_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Ray_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Segment_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Triangle_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Iso_rectangle_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Circle_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Conic_2, K >
{ static const int value = 2; };

template < typename K >
struct Dimension < typename K::Aff_transformation_2, K >
{ static const int value = 2; };


// 3D objects

template < typename K >
struct Dimension < typename K::Point_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Plane_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Vector_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Direction_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Line_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Ray_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Segment_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Triangle_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Tetrahedron_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Iso_cuboid_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Sphere_3, K >
{ static const int value = 3; };

template < typename K >
struct Dimension < typename K::Aff_transformation_3, K >
{ static const int value = 3; };


// dD objects

template < typename K >
struct Dimension < typename K::Point_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Vector_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Direction_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Hyperplane_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Aff_transformation_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Sphere_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Iso_box_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Segment_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Ray_d, K >
{ static const int value = 0; };

template < typename K >
struct Dimension < typename K::Line_d, K >
{ static const int value = 0; };

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_DIMENSION_H
