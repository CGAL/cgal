// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert

#ifndef CGAL_CH_FUNCTION_OBJECTS_2_H
#define CGAL_CH_FUNCTION_OBJECTS_2_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/enum.h>

namespace CGAL {

namespace internal {

template <class R>
class r_Less_dist_to_line
{
public:
  typedef bool    result_type;

  typedef typename R::Point_2  Point;
  typedef typename R::Line_2   Line;

        r_Less_dist_to_line() : line_constructed( false )
        { }

  bool  operator()(const Point& a, const Point& b,
                   const Point& c, const Point& d) const
        {
          if (!line_constructed)
          {
             line_constructed = true;
             l_ab = Line(a,b);
          }
          Comparison_result res = compare_signed_distance_to_line(l_ab, c, d);
          if ( res == SMALLER )
          {
              return true;
          }
          else if ( res == EQUAL )
          {
              return lexicographically_xy_smaller( c, d );
          }
          else
          {
              return false;
          }
        }

private:
  mutable bool line_constructed;
  mutable Line l_ab;
};

} // namespace internal

} //namespace CGAL

#endif // CGAL_CH_FUNCTION_OBJECTS_2_H
