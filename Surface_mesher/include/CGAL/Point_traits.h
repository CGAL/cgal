// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_TRAITS_H
#define CGAL_POINT_TRAITS_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Weighted_point_3.h>
#include <CGAL/assertions.h>

namespace CGAL {

  template <class P>
  struct Is_weighted : public Tag_false {} ;
  
  template <typename K>
  struct Is_weighted< ::CGAL::Weighted_point_3<K> > :
    public Tag_true {} ;
  
  namespace details {

    template <class P, typename FT, bool>
    struct Point_traits_aux 
    {
      // should give errors
    };

    template <class P, typename FT>
    struct Point_traits_aux<P, FT, false>
    {
      typedef P Point;
      typedef P Bare_point;
      typedef typename Kernel_traits<P>::type K;
      typedef typename ::CGAL::Weighted_point_3<K> Weighted_point;
      typedef Tag_false Is_weighted;
     
      const Bare_point& bare_point(const Point& bp)
      {
        return bp;
      }

      Weighted_point weighted_point(const Point& bp)
      {
        return Weighted_point(bp);
      }

      const Point& point(const Bare_point& bp)
      {
        return bp;
      }

      const Point& point(const Weighted_point& wp)
      {
        return wp.point();
      }
    }; // end class Point_traits_aux<P>

    template <class P, typename FT>
    struct Point_traits_aux<P, FT, true>
    {
      typedef P Point;
      typedef P Weighted_point;
      typedef typename Point::Point Bare_point;
      typedef Tag_true Is_weighted;
     
      const Bare_point& bare_point(const Point& wp)
      {
        return wp.point();
      }

      const Weighted_point& weighted_point(const Point& wp)
      {
        return wp;
      }

      Point point(const Bare_point& bp)
      {
        return Weighted_point(bp);
      }

      const Point& point(const Weighted_point& wp)
      {
        return wp;
      }
    }; // end class Point_traits_aux<P, FT, true>

    template <class Point>
    struct FT_of_point 
    {
      typedef typename CGAL::Kernel_traits<Point>::Kernel::FT type;
    };

  } // end namespace details

template <class Point>
class Point_traits
  : public details::Point_traits_aux<
      Point,
      typename details::FT_of_point<Point>::type,
      Is_weighted<Point>::value
    >
{
}; // end class Point_traits<T>

} // end namespace CGAL

#endif // CGAL_POINT_TRAITS_H
