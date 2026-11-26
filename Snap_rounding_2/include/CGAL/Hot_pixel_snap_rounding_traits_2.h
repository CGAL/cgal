// Copyright (c) 2001,2009,2014  Tel-Aviv University (Israel), Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>,
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_HOT_PIXEL_SNAP_ROUNDING_2_TRAITS_H
#define CGAL_HOT_PIXEL_SNAP_ROUNDING_2_TRAITS_H

#include <CGAL/license/Snap_rounding_2.h>

#include <cmath>

#include <CGAL/Arr_segment_traits_2.h>

namespace CGAL {

/*! \ingroup Snap_rounding_hot_pixel_grp
 *
 *  The class `Hot_pixel_snap_rounding_traits_2<Kernel>` is a model of the
 * `HotPixelSnapRoundingTraits_2` concept.  This class should be instantiated
 * with an exact geometric kernel that conforms to the \cgal kernel-concept that
 * is model of the concept `Kernel`.

 *
 * This geometric kernel must provide an (arbitrary-precision) rational number
 * type (`FT`), `Point_2`, and `Segment_2`. `FT` must be model of
 * `RealEmbeddableTraits`.
 *
 * The `CGAL::to_double()` function is used to implement the operation that
 * rounds the coordinates of a point to a center of a pixel. This operation is
 * one of the traits-concept requirements. The coordinates are converted to
 * double, rounded down to the nearest grid point, and finally adjusted to lie
 * on a center of a pixel. Notice that if `CGAL::to_double()` returns the
 * closest value, then when it rounds up a given coordinate, the resulting ISR,
 * may be imprecise, and the distance between some vertex and a non-incident
 * edge can be slightly less than the guaranteed half-the width-of-a-pixel.
 *
 * \cgalModels{HotPixelSnapRoundingTraits_2}
 *
 * \sa `CGAL::snap_rounding_2()`
 * \sa `CGAL::hot_pixel_snap_rounding_2()`
 */
template<typename K>
class Hot_pixel_snap_rounding_traits_2 :
  public CGAL::Arr_segment_traits_2<K> {
public:
  using NT = typename K::FT;
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Segment_2 = typename K::Segment_2;

  using Base_traits = CGAL::Arr_segment_traits_2<K>;
  using Compare_x_2 = typename Base_traits::Compare_x_2;

public:
  class Snap_2 {
  public:
    void operator()(const Point_2& p, NT pixel_size, NT& x, NT& y) {
      NT x_tmp = p.x() / pixel_size;
      NT y_tmp = p.y() / pixel_size;

      NT x_floor = std::floor(CGAL::to_double(x_tmp));
      NT y_floor = std::floor(CGAL::to_double(y_tmp));

      // fix the usage of to_double and floor:
      // there is no guarantee that the double approximation of *_tmp
      // has the same integer part as the NT version. The while loops
      // is a simple way to ensure that *_tmp are the true floor.
      while (x_floor > x_tmp) x_floor -= NT(1);
      while (y_floor > y_tmp) y_floor -= NT(1);
      while (x_floor + 1 <= x_tmp) x_floor += NT(1);
      while (y_floor + 1 <= y_tmp) y_floor += NT(1);

      CGAL_assertion(x_floor <= x_tmp);
      CGAL_assertion(y_floor <= y_tmp);
      CGAL_assertion(x_floor + 1 > x_tmp);
      CGAL_assertion(y_floor + 1 > y_tmp);

      x = x_floor * pixel_size + pixel_size / NT(2.0);
      y = y_floor * pixel_size + pixel_size / NT(2.0);
    }
  };

  Snap_2 snap_2_object() const { return Snap_2(); }

  /*! Functor */
  class Integer_grid_point_2 {
  public:
    Point_2 operator()(const Point_2& p, NT pixel_size) {
      NT x = (p.x() - pixel_size / NT(2.0)) / pixel_size;
      NT y = (p.y() - pixel_size / NT(2.0)) / pixel_size;
      Point_2 out_p(x,y);
      return out_p;
    }
  };

  /*! */
  Integer_grid_point_2 integer_grid_point_2_object() const {
    return Integer_grid_point_2();
  }

};

} // namespace CGAL

#endif
