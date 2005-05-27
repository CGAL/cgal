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

// Some tools to find the type of a kernel object with the same dimension
// as another given one.
// E.g. : Same_dimension_vector<K::Point_2>::type  is  K::Vector_2.
//
// Currently undocumented => for CGAL internal use only.

// TODO :
// - Do the other missing types.
// - What about functors ?  At least with a dimensional-insensitive interface.
// - Another thing which would be nice would be to make d-dimensional
//   algorithms use the 2d-3d kernel interface with a smaller wrapper.
//   (again, this forces a full kernel, not a sub-set traits, but at least...)
//   Then, based on Dimension<>, it's possibly to automatically use it.
//   cf Dimension_mapper<> at the end of the file.

template < typename T,
           typename K = typename Kernel_traits<T>::Kernel,
           int dim = Dimension<T, K>::value >
struct Same_dimension_point;

template < typename T, typename K >
struct Same_dimension_point <T, K, 2>
{ typedef typename K::Point_2  type; };

template < typename T, typename K >
struct Same_dimension_point <T, K, 3>
{ typedef typename K::Point_3  type; };

template < typename T, typename K >
struct Same_dimension_point <T, K, 0>
{ typedef typename K::Point_d  type; };


template < typename T,
           typename K = typename Kernel_traits<T>::Kernel,
           int dim = Dimension<T, K>::value >
struct Same_dimension_vector;

template < typename T, typename K >
struct Same_dimension_vector <T, K, 2>
{ typedef typename K::Vector_2  type; };

template < typename T, typename K >
struct Same_dimension_vector <T, K, 3>
{ typedef typename K::Vector_3  type; };

template < typename T, typename K >
struct Same_dimension_vector <T, K, 0>
{ typedef typename K::Vector_d  type; };


template < typename T,
           typename K = typename Kernel_traits<T>::Kernel,
           int dim = Dimension<T, K>::value >
struct Same_dimension_iso_box;

template < typename T, typename K >
struct Same_dimension_iso_box <T, K, 2>
{ typedef typename K::Iso_rectangle_2  type; };

template < typename T, typename K >
struct Same_dimension_iso_box <T, K, 3>
{ typedef typename K::Iso_cuboid_3  type; };

template < typename T, typename K >
struct Same_dimension_iso_box <T, K, 0>
{ typedef typename K::Iso_box_d  type; };


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
