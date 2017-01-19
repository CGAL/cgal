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

#ifndef CGAL_GPS_SEGMENT_TRAITS_2_H
#define CGAL_GPS_SEGMENT_TRAITS_2_H

#include <CGAL/license/Boolean_set_operations_2.h>



#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Polygon_2_curve_iterator.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

namespace CGAL {

template < class Kernel_, 
           class Container_ = std::vector<typename Kernel_::Point_2>,
           class Arr_seg_traits_ = Arr_segment_traits_2<Kernel_> >
class Gps_segment_traits_2 : public Arr_seg_traits_
{
  typedef Arr_seg_traits_                                               Base;
  typedef Gps_segment_traits_2<Kernel_, Container_, Arr_seg_traits_>    Self;

public:

  // Polygon_2 type is required by GeneralPolygonSetTraits Concept
  typedef CGAL::Polygon_2<Kernel_, Container_>          Polygon_2;
  // Polygon_2 is a model of the GeneralPolygon2 concept. 
  typedef  Polygon_2                                    General_polygon_2;

  // Polygon_with_holes_2 can be a simple polygon , with holes that are 
  // entirely inside him , or some vertices of the polygon and its holes
  // may overlap.
  
  // Polygon_with_holes_2 type required by GeneralPolygonSetTraits Concept.
  typedef CGAL::Polygon_with_holes_2<Kernel_, Container_>    
                                                Polygon_with_holes_2;
  // Polygon_with_Holes_2 is a model of the GeneralPolygonWithHoles2 concept. 
  typedef  Polygon_with_holes_2                 General_polygon_with_holes_2;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;

  typedef Polygon_2_curve_iterator<X_monotone_curve_2, Polygon_2>
                                                Curve_const_iterator;

  typedef typename Polygon_with_holes_2::Hole_const_iterator
                                                Hole_const_iterator;
  typedef typename Base::Point_2                Point_2;

 
  /*!
   * A functor for constructing a polygon from a range of segments.
   */
  class Construct_polygon_2 {
    typedef Gps_segment_traits_2<Kernel_, Container_, Arr_seg_traits_> Self;
    typedef Gps_traits_adaptor<Self>            Traits_adaptor;

    /*! The traits (in case it has state) */
    const Traits_adaptor* m_traits;
    
  public:
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Construct_polygon_2(const Self* traits) :
      m_traits(static_cast<const Traits_adaptor*>(traits))
    {}
    
    template <typename XCurveIterator>
    void operator()(XCurveIterator begin, XCurveIterator end, Polygon_2& pgn)
      const
    {
      typename Traits_adaptor::Construct_vertex_2 ctr_v =
        m_traits->construct_vertex_2_object();

      for (XCurveIterator itr = begin; itr != end; ++itr)
        pgn.push_back(ctr_v(*itr, 1));
    }
  };

  Construct_polygon_2 construct_polygon_2_object() const
  {
    return Construct_polygon_2(this);
  }

  /*!
   * A functor for scanning all segments that form a polygon boundary.
   */
  class Construct_curves_2 {
  public:
    std::pair<Curve_const_iterator, Curve_const_iterator>
    operator()(const General_polygon_2& pgn) const
    {
      Curve_const_iterator c_begin(&pgn, pgn.edges_begin());
      Curve_const_iterator c_end(&pgn, pgn.edges_end());

      return (std::make_pair(c_begin, c_end));
    }
  };

  Construct_curves_2 construct_curves_2_object() const
  {
    return Construct_curves_2();
  }

  // Added Functionality from GeneralPolygonWithHoles Concept to the traits.

  /* A functor for constructing the outer boundary of a polygon with holes
   */
  class Construct_outer_boundary {
  public:
    General_polygon_2 operator()(const  General_polygon_with_holes_2& pol_wh)
      const
    {
      return pol_wh.outer_boundary();
    }
  };

  Construct_outer_boundary construct_outer_boundary_object() const
  {
    return Construct_outer_boundary();
  }

  /* typedef from General_polygon_with_holes_2.
   * Hole_const_iterator nested type is required by
   * GeneralPolygonWithHoles2 concept
   */

  /*A functor for constructing the container of holes of a polygon with holes.
   * It returns ths begin/end iterators for the holes
   */
  class Construct_holes {
  public:
    std::pair<Hole_const_iterator, Hole_const_iterator>
    operator()(const General_polygon_with_holes_2& pol_wh) const
    {
      return std::make_pair(pol_wh.holes_begin(), pol_wh.holes_end());
    }
  };

  Construct_holes construct_holes_object() const
  {
    return Construct_holes();
  }
    
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
    General_polygon_with_holes_2 operator()(const General_polygon_2&
                                            pgn_boundary) const
    {
      return General_polygon_with_holes_2(pgn_boundary);
    }
    template <class HolesInputIterator>
    General_polygon_with_holes_2 operator()(const General_polygon_2&
                                            pgn_boundary,
                                            HolesInputIterator h_begin,
                                            HolesInputIterator h_end) const
    {
      return General_polygon_with_holes_2(pgn_boundary, h_begin,h_end);
    }
  };

  Construct_general_polygon_with_holes_2 construct_polygon_with_holes_2_object()
    const
  {
    return Construct_general_polygon_with_holes_2();
  }
  
  //functor returns true if the outer boundary is unbounded, and false otherwise.
  class Is_unbounded {
  public:
    bool operator()(const General_polygon_with_holes_2& pol_wh) const
    {
      return pol_wh.is_unbounded();
    }
  };

  Is_unbounded construct_is_unbounded_object() const
  {
    return Is_unbounded();
  }
};

} //namespace CGAL

#endif
