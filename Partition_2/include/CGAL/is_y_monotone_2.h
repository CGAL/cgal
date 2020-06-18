// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_IS_Y_MONOTONE_2_H
#define CGAL_IS_Y_MONOTONE_2_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/circulator.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <iterator>
#include <list>

namespace CGAL {

template<class InputIterator, class Traits>
bool
is_y_monotone_2(InputIterator first, InputIterator last,
                const Traits& traits)
{
   typedef typename Traits::Point_2                  Point_2;
   typedef typename Traits::Less_yx_2                Less_yx_2;
   typedef std::list<Point_2>                        Vertex_list;
   typedef typename Vertex_list::iterator            I;
   typedef Circulator_from_iterator<I>               Circulator;

   Vertex_list  polygon(first, last);
   I n_point_ref;
   I s_point_ref;

   // find the extreme north and south points
   ch_n_point(polygon.begin(), polygon.end(), n_point_ref, traits);
   ch_s_point(polygon.begin(), polygon.end(), s_point_ref, traits);

   Circulator north_c(polygon.begin(), polygon.end(), n_point_ref);
   Circulator south_c(polygon.begin(), polygon.end(), s_point_ref);

   Less_yx_2 less_yx = traits.less_yx_2_object();

   Circulator prev_c = south_c;
   Circulator c = south_c;

   for (c++; c != north_c; c++)
   {
      if (!less_yx(*prev_c, *c)) return false;
      prev_c = c;
   }

   prev_c = north_c;
   c = north_c;
   for (c++; c != south_c; c++)
   {
      if (!less_yx(*c, *prev_c)) return false;
      prev_c = c;
   }
   return true;
}

template<class InputIterator>
bool
is_y_monotone_2(InputIterator first, InputIterator last)
{
   typedef typename std::iterator_traits<InputIterator>::value_type Point_2;
   typedef typename Kernel_traits<Point_2>::Kernel  K;
   return is_y_monotone_2(first, last, K());
}

}

#endif // CGAL_IS_Y_MONOTONE_2_H
