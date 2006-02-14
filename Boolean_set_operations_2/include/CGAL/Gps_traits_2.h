// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Revision$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_TRAITS_2_H
#define GPS_TRAITS_2_H

#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>

CGAL_BEGIN_NAMESPACE

template <class Arr_traits,
          class Gen_polygon_ = General_polygon_2<Arr_traits> >
class Gps_traits_2 : public Arr_traits
{
  typedef Arr_traits                                    Base;
  typedef Gps_traits_2<Arr_traits>                      Self;
  
  public:
    typedef typename Base::Point_2                       Point_2;
    typedef typename Base::X_monotone_curve_2            X_monotone_curve_2;
    typedef Gen_polygon_                                 Polygon_2;
    typedef General_polygon_with_holes_2<Polygon_2>      Polygon_with_holes_2;
    typedef typename Polygon_2::Curve_const_iterator     Curve_const_iterator;
    typedef typename Base::Equal_2                       Equal_2;
    typedef typename Base::Compare_endpoints_xy_2        Compare_endpoints_xy_2;
    typedef typename Base::Construct_min_vertex_2        Construct_min_vertex_2;
    typedef typename Base::Construct_max_vertex_2        Construct_max_vertex_2;


  class Construct_polygon_2
  {
  public:

    template<class XCurveIterator>
    void operator()(XCurveIterator begin,
                    XCurveIterator end,
                    Polygon_2& pgn)
    {
      pgn.init(begin, end);
    }
  };

  Construct_polygon_2 construct_polygon_2_object() const
  {
    return Construct_polygon_2();
  }

  class Construct_curves_2
  {
  public:

    std::pair<Curve_const_iterator,
              Curve_const_iterator> operator()(const Polygon_2& pgn)
    {
      return std::make_pair(pgn.curves_begin(), pgn.curves_end());
    }
  };

  Construct_curves_2 construct_curves_2_object()
  {
    return Construct_curves_2();
  }

  typedef Is_valid_2<Self>    Is_valid_2;
 
  Is_valid_2 is_valid_2_object()
  {
    return Is_valid_2();
  }
};

CGAL_END_NAMESPACE

#endif
