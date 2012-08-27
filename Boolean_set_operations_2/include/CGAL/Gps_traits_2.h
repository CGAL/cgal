// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_TRAITS_2_H
#define CGAL_GPS_TRAITS_2_H

#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

namespace CGAL {

template <typename Arr_traits,
          typename General_polygon_t = General_polygon_2<Arr_traits> >
class Gps_traits_2 : public Arr_traits
{
  typedef Arr_traits                                    Base;
  typedef Gps_traits_2<Arr_traits,General_polygon_t>    Self;
  
public:

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  //Polygon_2 type is required by GeneralPolygonSetTraits Concept    
  typedef General_polygon_t                             Polygon_2;
  //Polygon_2 is a model of the GeneralPolygon2 concept
  typedef Polygon_2                                     General_polygon_2;
  
  //Polygon_with_holes_2 type required by GeneralPolygonSetTraits Concept.
  typedef CGAL::General_polygon_with_holes_2<General_polygon_2>
                                                        Polygon_with_holes_2;
  //Polygon_with_Holes_2 is a model of the GeneralPolygonWithHoles2 concept.
  typedef Polygon_with_holes_2
    General_polygon_with_holes_2;
  
  typedef typename General_polygon_2::Curve_const_iterator
                                                        Curve_const_iterator;
  
  typedef typename General_polygon_with_holes_2::Hole_const_iterator
                                                        Hole_const_iterator;
                                           
  typedef typename Base::Equal_2                        Equal_2;
  typedef typename Base::Compare_endpoints_xy_2         Compare_endpoints_xy_2;
  typedef typename Base::Construct_min_vertex_2         Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2         Construct_max_vertex_2;


  /*!
   * A functor for constructing a polygon from a range of x-monotone curves.
   */
  class Construct_polygon_2 {
  public:
    template<class XCurveIterator>
    void operator()(XCurveIterator begin, XCurveIterator end,
                    General_polygon_2& pgn)
    { pgn.init(begin, end); }
  };

  Construct_polygon_2 construct_polygon_2_object() const
  { return Construct_polygon_2(); }

  /*!
   * A functor for scanning all x-monotone curves that form a polygon boundary.
   */
  class Construct_curves_2
  {
  public:

    std::pair<Curve_const_iterator, Curve_const_iterator>
    operator()(const General_polygon_2& pgn)
    { return std::make_pair(pgn.curves_begin(), pgn.curves_end()); }
  };

  Construct_curves_2 construct_curves_2_object()
  { return Construct_curves_2(); }

  /*!
   * An auxiliary functor used for validity checks.
   */
  typedef Gps_traits_adaptor<Base>                      Traits_adaptor;
 
  /*typedef CGAL::Is_valid_2<Self, Traits_adaptor>           Is_valid_2;
    Is_valid_2 is_valid_2_object()
    {
    Traits_adaptor   tr_adp;
 
    return (Is_valid_2 (*this, tr_adp));
    }*/  
  
  //Added Functionality from GeneralPolygonWithHoles Concept to the traits. 
  
  /*A functor for constructing the outer boundary of a polygon with holes*/   
  class Construct_outer_boundary {
  public:
    General_polygon_2 operator()(const  General_polygon_with_holes_2& pol_wh) 
    { return pol_wh.outer_boundary(); }
  };
  
  Construct_outer_boundary construct_outer_boundary_object() const
  { return Construct_outer_boundary(); }
  
  /* typedef from General_polygon_with_holes_2. Hole_const_iterator nested type
   * is required by GeneralPolygonWithHoles2 concept
   */
  /*A functor for constructing the container of holes of a polygon with holes*/  
  class Construct_holes {
  public:
    std::pair<Hole_const_iterator, Hole_const_iterator>
    operator()(const General_polygon_with_holes_2& pol_wh) 
    { return std::make_pair(pol_wh.holes_begin(), pol_wh.holes_end()); }
  };
  
  Construct_holes construct_holes_object() const
  { return Construct_holes(); }
    
  /* A functor for constructing a General_polygon_with_holes from a
   * General_Polygon (and possibly a range of holes).
   *
   * constructs a general polygon with holes using a given general polygon
   * outer as the outer boundary and a given range of holes. If outer is an
   * empty general polygon, then an unbounded polygon with holes will be
   * created. The holes must be contained inside the outer boundary, and the
   * polygons representing the holes must be strictly simple and pairwise
   * disjoint, except perhaps at the vertices.
   */
  class Construct_general_polygon_with_holes_2 {
  public:     
    General_polygon_with_holes_2
    operator()(const General_polygon_2& pgn_boundary) 
    { return General_polygon_with_holes_2(pgn_boundary); }

    template <class HolesInputIterator>
    General_polygon_with_holes_2
    operator()(const General_polygon_2& pgn_boundary,
               HolesInputIterator h_begin,
               HolesInputIterator h_end)
    { return General_polygon_with_holes_2(pgn_boundary, h_begin,h_end); }
  };

  Construct_general_polygon_with_holes_2 construct_polygon_with_holes_2_object()
    const
  { return Construct_general_polygon_with_holes_2(); }
  
  // Return true if the outer boundary is empty, and false otherwise.
  class Is_unbounded {
  public:
    bool operator()(const  General_polygon_with_holes_2& pol_wh) 
    { return pol_wh.is_unbounded(); }  
  };
  
  Is_unbounded construct_is_unbounded_object() { return Is_unbounded(); }
};

} //namespace CGAL

#endif
