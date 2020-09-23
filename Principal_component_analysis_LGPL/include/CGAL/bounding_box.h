// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_BOUNDING_BOX_H
#define CGAL_BOUNDING_BOX_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel/Dimension_utils.h>

namespace CGAL {

// TODO :
// - Add support for more types of objects.
// - Generalize to d-dimension.
// - Maybe the internal code itself could be factorized across dimensions...

namespace internal {

// internal functions specialized by dimension.

template < class ForwardIterator, class Traits >
typename Traits::Iso_rectangle_2
bounding_box(ForwardIterator f, ForwardIterator l, const Traits& t, Dimension_tag<2>)
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
bounding_box(ForwardIterator f, ForwardIterator l, const Traits& t, Dimension_tag<3>)
{
  CGAL_precondition(f != l);
  typedef typename Traits::Less_x_3                  Less_x_3;
  typedef typename Traits::Less_y_3                  Less_y_3;
  typedef typename Traits::Less_z_3                  Less_z_3;
  typedef typename Traits::Construct_iso_cuboid_3    Cub;

  Less_x_3 lessx = t.less_x_3_object();
  Less_y_3 lessy = t.less_y_3_object();
  Less_z_3 lessz = t.less_z_3_object();
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
bounding_box(ForwardIterator f, ForwardIterator l, const Traits& t, Dynamic_dimension_tag);

// To be written...

#endif

}

template < class ForwardIterator, class K >
inline
typename Access::Iso_box<K, typename Ambient_dimension<typename std::iterator_traits<ForwardIterator>
                                              ::value_type, K>::type>::type
bounding_box(ForwardIterator f, ForwardIterator l, const K& k)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type Pt;
  return internal::bounding_box(f, l, k, typename Ambient_dimension<Pt>::type() );
}

template < class ForwardIterator >
inline
typename Access::Iso_box<typename Kernel_traits<typename std::iterator_traits<ForwardIterator>::value_type>::Kernel,
                   typename Ambient_dimension<typename std::iterator_traits<ForwardIterator>::value_type,
                   typename Kernel_traits<typename std::iterator_traits<ForwardIterator>::value_type>::Kernel>::type >::type
bounding_box(ForwardIterator f, ForwardIterator l)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type Pt;
  typedef typename Kernel_traits<Pt>::Kernel Kernel;
  return bounding_box(f, l, Kernel());
}

} //namespace CGAL

#endif // CGAL_BOUNDING_BOX_H
