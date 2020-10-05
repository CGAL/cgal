// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SEGMENT_ASSERTIONS_H
#define CGAL_SEGMENT_ASSERTIONS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


namespace CGAL {

template <class Traits_>
class Segment_assertions
{
  typedef Traits_                                  Traits_2;

  typedef typename Traits_2::Point_2               Point_2;
  typedef typename Traits_2::Kernel                Kernel;
  typedef typename Kernel::Line_2                  Line_2;
  typedef typename Traits_2::X_monotone_curve_2    X_monotone_curve_2;

public:

  static bool _assert_is_point_on (const Point_2& pt,
                                   const X_monotone_curve_2& cv,
                                   Tag_true /* tag */)
  {
    Traits_2    traits;
    return (traits.compare_y_at_x_2_object() (pt, cv) == EQUAL);
  }

  static bool _assert_is_point_on (const Point_2& /* pt */,
                                   const X_monotone_curve_2& /* cv */,
                                   Tag_false /* tag */)
  {
    return (true);
  }

  static bool _assert_is_point_on (const Point_2& pt,
                                   const Line_2& l,
                                   Tag_true /* tag */)
  {
    Kernel      kernel;
    return (kernel.has_on_2_object() (l, pt));
  }

  static bool _assert_is_point_on (const Point_2& /* pt */,
                                   const Line_2& /* l */,
                                   Tag_false /* tag */)
  {
    return (true);
  }
};

} //namespace CGAL

#endif
