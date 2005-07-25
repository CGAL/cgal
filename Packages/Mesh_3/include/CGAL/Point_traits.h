// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_TRAITS_H_H
#define CGAL_POINT_TRAITS_H_H

#include <CGAL/Weighted_point.h>

namespace CGAL {

  namespace details {

    template <class P, typename FT>
    struct Point_traits_aux
    {
      typedef P Point;
      typedef P Bare_point;
      typedef typename ::CGAL::Weighted_point<Bare_point, FT> Weighted_point;
      typedef Tag_false Has_weithed_point;
     
      const Bare_point& bare_point(const Point& bp)
      {
        return bp;
      }

      Weighted_point weighted_point(const Point& bp)
      {
        return weighted_point(bp);
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
    struct Point_traits_aux< ::CGAL::Weighted_point<P, FT>, FT >
    {
      typedef ::CGAL::Weighted_point<P, FT> Point;
      typedef Point Weighted_point;
      typedef typename Point::Point Bare_point;
      typedef Tag_true Has_weithed_point;
     
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
    }; // end class Point_traits_aux<P>

  } // end namespace details

template <class Point>
class Point_traits : 
    public details::Point_traits_aux<Point,
       typename CGAL::Kernel_traits<Point>::Kernel::FT >
{
}; // end class Point_traits<T>

} // end namespace CGAL

#endif // CGAL_POINT_TRAITS_H_H
