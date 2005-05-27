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

#ifndef CGAL_BOUNDING_BOX_H
#define CGAL_BOUNDING_BOX_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/Dimension.h>
#include <CGAL/Kernel/Dimension_utils.h>

CGAL_BEGIN_NAMESPACE

// TODO :
// - Generalize to d-dimension.
// - Maybe the internal code itself could be factorized against dimensions...

namespace CGALi {

// internal functions specialized by dimension.

template < class ForwardIterator, class Traits >
typename Traits::Iso_rectangle_2
bounding_box_2(ForwardIterator f, ForwardIterator l, const Traits& t)
{
  CGAL_precondition(f != l);
  typedef typename Traits::Less_x_2                  Less_x_2;
  typedef typename Traits::Less_y_2                  Less_y_2;
  typedef typename Traits::Construct_iso_rectangle_2 Rect;

  Less_x_2 lessx = t.less_x_2_object();
  Less_y_2 lessy = t.less_y_2_object();
  Rect     rect  = t.construct_iso_rectangle_2_object();

  ForwardIterator xmin = f;
  ForwardIterator xmax = f;
  ForwardIterator ymin = f;
  ForwardIterator ymax = f;

  while (++f != l) {
    if (lessx(*f, *xmin))
      xmin = f;
    else if (lessx(*xmax, *f))
      xmax = f;

    if (lessy(*f, *ymin))
      ymin = f;
    else if (lessy(*ymax, *f))
      ymax = f;
  }

  return rect(*xmin, *xmax, *ymin, *ymax);
}

template < class ForwardIterator, class Traits >
typename Traits::Iso_cuboid_3
bounding_box_3(ForwardIterator f, ForwardIterator l, const Traits& t)
{
  CGAL_precondition(f != l);
  typedef typename Traits::Less_x_3                  Less_x_3;
  typedef typename Traits::Less_y_3                  Less_y_3;
  typedef typename Traits::Less_z_3                  Less_z_3;
  typedef typename Traits::Construct_iso_cuboid_3    Cub;

  Less_x_3 lessx = t.less_x_3_object();
  Less_y_3 lessy = t.less_y_3_object();
  Less_z_3 lessy = t.less_z_3_object();
  Cub      cub  = t.construct_iso_cuboid_3_object();

  ForwardIterator xmin = f;
  ForwardIterator xmax = f;
  ForwardIterator ymin = f;
  ForwardIterator ymax = f;
  ForwardIterator zmin = f;
  ForwardIterator zmax = f;

  while (++f != l) {
    if (lessx(*f, *xmin))
      xmin = f;
    else if (lessx(*xmax, *f))
      xmax = f;

    if (lessy(*f, *ymin))
      ymin = f;
    else if (lessy(*ymax, *f))
      ymax = f;

    if (lessz(*f, *zmin))
      zmin = f;
    else if (lessz(*zmax, *f))
      zmax = f;
  }

  return cub(*xmin, *xmax, *ymin, *ymax, *zmin, *zmax);
}

#if 0
template < class ForwardIterator, class Traits >
typename Traits::Iso_box_d
bounding_box_d(ForwardIterator f, ForwardIterator l, const Traits& t);

// To be written...

#endif


template < int dim >
struct bbox;

template <>
struct bbox<2>
{
  template < class ForwardIterator, class Traits >
  typename Traits::Iso_rectangle_2
  operator()(ForwardIterator f, ForwardIterator l, const Traits& t) const
  {
    return bounding_box_2(f, l, t);
  }
};

template <>
struct bbox<3>
{
  template < class ForwardIterator, class Traits >
  typename Traits::Iso_cuboid_3
  operator()(ForwardIterator f, ForwardIterator l, const Traits& t) const
  {
    return bounding_box_3(f, l, t);
  }
};

template <>
struct bbox<0>
{
  template < class ForwardIterator, class Traits >
  typename Traits::Iso_box_d
  operator()(ForwardIterator f, ForwardIterator l, const Traits& t) const
  {
    return bounding_box_d(f, l, t);
  }
};

}

template < class ForwardIterator >
inline
typename Same_dimension_iso_box<typename std::iterator_traits< ForwardIterator >
                                ::value_type>::type
bounding_box(ForwardIterator f, ForwardIterator l)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type Pt;
  typedef typename Kernel_traits<Pt>::Kernel Kernel;
  return bounding_box(f, l, Kernel());
}

template < class ForwardIterator, class Traits >
inline
typename Same_dimension_iso_box<typename std::iterator_traits< ForwardIterator >
                                ::value_type, Traits>::type
bounding_box(ForwardIterator f, ForwardIterator l, const Traits& t)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type Pt;
  typedef typename Kernel_traits<Pt>::Kernel Kernel;
  return bounding_box(f, l, Kernel());
}

CGAL_END_NAMESPACE

#endif // CGAL_BOUNDING_BOX_H
