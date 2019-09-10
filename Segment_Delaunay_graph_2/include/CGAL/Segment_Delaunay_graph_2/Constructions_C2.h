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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCTIONS_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCTIONS_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/enum.h>

#include <CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h>

#include <CGAL/Parabola_2.h>
#include <CGAL/Parabola_segment_2.h>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {


//***********************************************************************
//***********************************************************************
//                            CONSTRUCTIONS
//***********************************************************************
//***********************************************************************


//-----------------------------------------------------------------------
//                  Segment Delaunay graph site
//-----------------------------------------------------------------------
template<class Site,class ITag> class Construct_sdg_site_2;

template<class Site>
class Construct_sdg_site_2<Site,Tag_true>
{
public:
  typedef Site                             Site_2;
  typedef typename Site_2::Point_2         Point_2;
  typedef Site_2                           result_type;

public:
  result_type operator()(const Point_2& p) const {
    return Site_2(p);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1) const {
    return Site_2(p0, p1);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
			 const Point_2& q0, const Point_2& q1) const {
    return Site_2(p0, p1, q0, q1);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
			 const Point_2& q0, const Point_2& q1,
			 bool b) const {
    return Site_2(p0, p1, q0, q1, b);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1,
			 const Point_2& q0, const Point_2& q1,
			 const Point_2& r0, const Point_2& r1) const {
    return Site_2(p0, p1, q0, q1, r0, r1);
  }
};


template<class Site>
class Construct_sdg_site_2<Site,Tag_false>
{
public:
  typedef Site                             Site_2;
  typedef typename Site_2::Point_2         Point_2;
  typedef Site_2                           result_type;

public:
  result_type operator()(const Point_2& p) const {
    return Site_2(p);
  }

  result_type operator()(const Point_2& p0, const Point_2& p1) const {
    return Site_2(p0, p1);
  }
};




//-----------------------------------------------------------------------
//                  Segment Delaunay graph Voronoi vertex
//-----------------------------------------------------------------------

template<class K, class M>
class Construct_svd_vertex_2
{
public:
  typedef typename K::Site_2                Site_2;
  typedef Voronoi_vertex_C2<K,M>            Voronoi_vertex_2;
  typedef typename K::Point_2               Point_2;
  typedef Point_2                           result_type;

public:
  Point_2 operator()(const Site_2& s1, const Site_2& s2,
		     const Site_2& s3) const
  {
    Voronoi_vertex_2 v(s1, s2, s3);
    return v.point();
  }
};


//-----------------------------------------------------------------------
//                  Segment Delaunay graph Voronoi circle
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_circle_2
{
public:
  typedef typename Gt::Site_2                 Site_2;
  typedef Voronoi_vertex_C2<Gt,M>             Voronoi_vertex_2;
  typedef typename Gt::Circle_2               Circle_2;
  typedef Circle_2                            result_type;

public:
  Circle_2 operator() (const Site_2& s1, const Site_2& s2,
		       const Site_2& s3) const
  {
    Voronoi_vertex_2 v(s1, s2, s3);
    return v.circle();
  }
};



//-----------------------------------------------------------------------
//                    Segment Delaunay graph Voronoi bisector
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_bisector_2
{
public:
  typedef typename Gt::Site_2        Site_2;
  typedef typename Gt::Point_2       Point_2;
  typedef typename Gt::Line_2        Line_2;
  typedef Line_2                     result_type;

private:
  static
  Point_2 midpoint(const Point_2& p, const Point_2& q, Integral_domain_without_division_tag) {
    typedef typename Gt::FT  FT;
    FT half(0.5);
    return Point_2((p.x() + q.x()) * half,(p.y() + q.y()) * half);
  }

  static
  Point_2 midpoint(const Point_2& p, const Point_2& q, Field_tag) {
    return CGAL::midpoint(p, q);
  }

  static Point_2 midpoint(const Point_2& p, const Point_2& q) {
    typedef typename Gt::FT  FT;
    typedef Algebraic_structure_traits<FT> AST;
    return midpoint(p, q, typename AST::Algebraic_category());
  }

public:
  Line_2 operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    if ( p.is_point() && q.is_point() ) {
      Point_2 mid = midpoint(p.point(), q.point());
      Line_2 l(p.point(), q.point());
      return l.perpendicular(mid);
    }
    if ( p.is_segment() && q.is_point() ) {
      // in this case q has to be one of the two endpoints of the
      // segment p...
      Line_2 l = p.segment().supporting_line();
      return l.perpendicular(q.point());
    }
    // in this case p has to be one of the two endpoints of the
    // segment q...
    Line_2 l = q.segment().supporting_line();
    return l.perpendicular(p.point());
  }
};

//-----------------------------------------------------------------------
//                 Segment Delaunay graph Voronoi bisecting ray
//-----------------------------------------------------------------------

template<class Gt, class M>
class Construct_sdg_bisector_ray_2
{
public:
  typedef typename Gt::Site_2                   Site_2;
  typedef typename Gt::Point_2                  Point_2;
  typedef typename Gt::Line_2                   Line_2;
  typedef typename Gt::Ray_2                    Ray_2;
  typedef typename Gt::Construct_svd_vertex_2   Construct_svd_vertex_2;
  typedef typename Gt::Equal_2                  Equal_2;
  typedef Ray_2                                 result_type;

  Ray_2 operator()(const Site_2& p, const Site_2& q,
		   const Site_2& r) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    Equal_2 are_same_points;

    Point_2 v = Construct_svd_vertex_2()(p, q, r);
    Point_2 p1, p2;
    if ( p.is_point() && q.is_point() ) {
      p1 = q.point();
      p2 = p.point();
    } else if ( p.is_point() && q.is_segment() ) {
      CGAL_assertion( are_same_points(p, q.source_site()) ||
		      are_same_points(p, q.target_site()) );
      p1 = are_same_points(p, q.source_site()) ? q.target() : q.source();
      p2 = p.point();
    } else {
      // p is a segment and q a point
      p1 = q.point();
      p2 = are_same_points(q, p.source_site()) ? p.target() : p.source();
    }
    Line_2 l(p1, p2);
    Line_2 lperp = l.perpendicular( v );
    return Ray_2(v, lperp.direction());
  }
};


//-----------------------------------------------------------------------
//              Segment Delaunay graph Voronoi bisecting segment
//-----------------------------------------------------------------------


template<class Gt, class M>
class Construct_sdg_bisector_segment_2
{
public:
  typedef typename Gt::Site_2                  Site_2;
  typedef typename Gt::Point_2                 Point_2;
  typedef typename Gt::Line_2                  Line_2;
  typedef typename Gt::Ray_2                   Ray_2;
  typedef typename Gt::Segment_2               Segment_2;
  typedef CGAL::Parabola_segment_2<Gt>         Parabola_segment_2;

  typedef typename Gt::Construct_svd_vertex_2  Construct_svd_vertex_2;
  typedef typename Gt::Equal_2                 Equal_2;

  typedef CGAL::Object                         Object_2;
  typedef Object_2                             result_type;

  result_type operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r, const Site_2& s) const
  {
    Construct_svd_vertex_2 circumcenter;
    Point_2 vpqr = circumcenter(p, q, r);
    Point_2 vqps = circumcenter(q, p, s);

    Equal_2 same_points;

    if ( (p.is_point() && q.is_point()) ||
	 (p.is_segment() && q.is_segment()) ) {
      Segment_2 vorseg(vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    if ( p.is_point() ) {
      // check is p is an endpoint of q
      if (  same_points( p, q.source_site() ) ||
	    same_points( p, q.target_site() )  ) {
	Segment_2 vorseg(vpqr, vqps);
	return CGAL::make_object(vorseg);
      }
      Line_2 l = q.segment().supporting_line();
      Parabola_segment_2 vorseg(p.point(), l, vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    // check is q is an endpoint of p
    if ( same_points(q, p.source_site()) ||
	 same_points(q, p.target_site()) ) {
      Segment_2 vorseg(vpqr, vqps);
      return CGAL::make_object(vorseg);
    }
    Line_2 l = p.segment().supporting_line();
    Parabola_segment_2 vorseg(q.point(), l, vpqr, vqps);
    return CGAL::make_object(vorseg);
  }
};

//-----------------------------------------------------------------------


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCTIONS_C2_H
