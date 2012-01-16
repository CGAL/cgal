// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_TYPE_NON_INTERSECTING_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_TYPE_NON_INTERSECTING_C2_H

#include <CGAL/enum.h>
#include <CGAL/determinant.h>
#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h>
#include <CGAL/Segment_Delaunay_graph_2/Arrangement_enum.h>

#include <iostream>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class K>
class Arrangement_type_non_intersecting_C2
  : public Basic_predicates_C2<K>, public Internal::Arrangement_enum
{
private:
  typedef Internal::Arrangement_enum         Enum;

  typedef Basic_predicates_C2<K>             Base;

  typedef typename Base::Point_2             Point_2;
  typedef typename Base::Segment_2           Segment_2;
  typedef typename Base::Site_2              Site_2;
  typedef typename Base::Line_2              Line_2;

  typedef typename Base::FT                  FT;
  typedef typename Base::RT                  RT;

  typedef typename Base::Orientation         Orientation;
  typedef typename Base::Comparison_result   Comparison_result;
  typedef typename Base::Oriented_side       Oriented_side;
  typedef typename Base::Sign                Sign;
  typedef typename Base::Boolean             Boolean;

  typedef typename K::Orientation_2          Orientation_2;

  typedef Are_same_points_C2<K>              Are_same_points_2;

private:
  Are_same_points_2     same_points;

public:
  typedef typename Enum::Arrangement_type    result_type;

private:

  //--------------------------------------------------------------------

  result_type
  arrangement_type_ss(const Site_2& p, const Site_2& q) const
  {
    bool same_p1q1 = same_points(p.source_site(), q.source_site());
    bool same_p1q2 = same_points(p.source_site(), q.target_site());
    bool same_p2q1 = same_points(p.target_site(), q.source_site());
    bool same_p2q2 = same_points(p.target_site(), q.target_site());

    if ( (same_p1q1 && same_p2q2) || (same_p1q2 && same_p2q1) ) {
      return IDENTICAL;
    }

    if ( same_p1q1 ) {
      return TOUCH_11;
    } if ( same_p1q2 ) {
      return TOUCH_12;
    } else if ( same_p2q1 ) {
      return TOUCH_21;
    } else if ( same_p2q2 ) {
      return TOUCH_22;
    }

    return DISJOINT;
  }

  //--------------------------------------------------------------------

  result_type
  arrangement_type_ps(const Site_2& p, const Site_2& q) const
  {
    if ( same_points(p, q.source_site()) ) {
      return TOUCH_1;
    } else if ( same_points(p, q.target_site()) ) {
      return TOUCH_2;
    } else {
      return DISJOINT;
    }
  }

  //--------------------------------------------------------------------

  result_type
  arrangement_type_pp(const Site_2& p, const Site_2& q) const
  {
    if ( same_points(p, q) ) {
      return IDENTICAL;
    } else {
      return DISJOINT;
    }
  }

  //--------------------------------------------------------------------

public:
  typedef Site_2                   argument_type;


  result_type
  operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_defined() && q.is_defined() );

    if ( p.is_point() ) {
      if ( q.is_point() ) {
	return arrangement_type_pp(p, q);
      } else {
	return arrangement_type_ps(p, q);
      }
    } else {
      if ( q.is_point() ) {
	return opposite( arrangement_type_ps(q, p) );
      } else {
	return arrangement_type_ss(p, q);
      }
    }
  }
};

//---------------------------------------------------------------------

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_TYPE_NON_INTERSECTING_C2_H
