// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_INTERSECTION_H
#define CGAL_SNC_INTERSECTION_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Circulator_project.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 37
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template < class Node, class Object>
struct Project_shalfedge_point {
  typedef Node         argument_type;
  typedef Object       result_type;
  Object& operator()( Node& x) const   {
    return x.source()->source()->point();
    /* a Point_3& reference must be returned by D.point() */
  }
  const Object& operator()( const Node& x) const   {
    return x.source()->source()->point();
    /* a Point_3& reference must be returned by D.point() */
  }
};

template<typename SNC_structure_>
class SNC_intersection {

  typedef SNC_structure_                     SNC_structure;
  typedef SNC_intersection<SNC_structure>    Self;

  typedef typename SNC_structure::SHalfedge               SHalfedge;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_const_handle
                                  Halffacet_const_handle;
  typedef typename SNC_structure::SHalfedge_const_handle  SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfloop_const_handle  SHalfloop_const_handle;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
                                  SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator
                                  Halffacet_cycle_const_iterator;

  typedef typename SNC_structure::Point_3        Point_3;
  typedef typename SNC_structure::Vector_3       Vector_3;
  typedef typename SNC_structure::Segment_3      Segment_3;
  typedef typename SNC_structure::Line_3         Line_3;
  typedef typename SNC_structure::Ray_3          Ray_3;
  typedef typename SNC_structure::Plane_3        Plane_3;

 public:

  static bool does_contain_internally(const Point_3& s,
                                      const Point_3& t,
                                      const Point_3& p) {
    return are_strictly_ordered_along_line (s, p, t);
  }

  static bool does_contain_internally(Halffacet_const_handle f,
                                      const Point_3& p) {
    if(!f->plane().has_on(p))
      return false;
    return point_in_facet_interior( p, f);
  }

  static bool does_intersect_internally(const Segment_3& s1,
                                        const Segment_3& s2,
                                        Point_3& p) {
    if(s2.has_on(s1.target()))
      return false;
    Ray_3 r(s1.source(), s1.target());
    if(!does_intersect_internally(r, s2, p))
      return false;
    Plane_3 pl(s1.target(), r.to_vector());
    return (pl.oriented_side(p) == CGAL::NEGATIVE);
  }

  static bool does_intersect_internally(const Ray_3& s1,
                                        const Segment_3& s2,
                                        Point_3& p) {
    if (!coplanar( s1.source(), s1.point(1), s2.source(), s2.target()))
      // the segments doesn't define a plane
      return false;
    if ( s1.has_on(s2.source()) || s1.has_on(s2.target()) ||
         s2.has_on(s1.source()))
      // the segments does intersect at one endpoint
      return false;
    Line_3 ls1(s1), ls2(s2);
    if ( ls1.direction() ==  ls2.direction() ||
         ls1.direction() == -ls2.direction() )
      // the segments are parallel
      return false;
    Vector_3 vs1(s1.to_vector()), vs2(s2.to_vector()),
      vt(cross_product( vs1, vs2)),
      ws1(cross_product( vt, vs1));
    Plane_3 hs1( s1.source(), ws1);
    Object o = intersection(hs1, ls2);
    CGAL_assertion(CGAL::assign( p, o));
    // since line(s1) and line(s2) are not parallel they intersects in only
    //   one point
    CGAL::assign( p ,o);
    Plane_3 pl(s1.source(), vs1);
    if(pl.oriented_side(p) != CGAL::POSITIVE)
      return false;
    pl = Plane_3(s2.source(), vs2);
    if(pl.oriented_side(p) != CGAL::POSITIVE)
      return false;
    pl = Plane_3(s2.target(), vs2);
    return (pl.oriented_side(p) == CGAL::NEGATIVE);
  }

  static bool does_intersect_internally(const Ray_3& ray,
                                        Halffacet_const_handle f,
                                        Point_3& p) {
    CGAL_NEF_TRACEN("-> Intersection facet - ray");
    Plane_3 h( f->plane());
    CGAL_NEF_TRACEN("-> facet's plane: " << h);
    CGAL_NEF_TRACEN("-> a point on the plane: " << h.point());
    CGAL_NEF_TRACEN("-> ray: " << ray);
    CGAL_assertion(!ray.is_degenerate());
    if(h.has_on(ray.source()))
      return false;
    Object o = intersection( h, ray);
    if( !CGAL::assign( p, o))
      return false;
    CGAL_NEF_TRACEN( "-> intersection point: " << p );
    CGAL_NEF_TRACEN( "-> point in facet interior? "<<point_in_facet_interior( f, p));
    return point_in_facet_interior( p, f);
  }

  static bool does_intersect_internally(const Segment_3& seg,
                                        Halffacet_const_handle f,
                                        Point_3& p) {
    CGAL_NEF_TRACEN("-> Intersection facet - segment");
    Plane_3 h( f->plane());
    CGAL_NEF_TRACEN("-> facet's plane: " << h);
    CGAL_NEF_TRACEN("-> a point on the plane: " << h.point());
    CGAL_NEF_TRACEN("-> segment: " << seg);
    CGAL_assertion(!seg.is_degenerate());
    if( h.has_on( seg.source()) || h.has_on(seg.target()))
      /* no possible internal intersection */
      return false;
    Object o = intersection( h, seg);
    if( !CGAL::assign( p, o))
      return false;
    CGAL_NEF_TRACEN( "-> intersection point: " << p );
    CGAL_NEF_TRACEN( "-> point in facet interior? "<<point_in_facet_interior( f, p));
    return point_in_facet_interior( p, f);
  }

 private:

  static bool point_in_facet_interior(const Point_3& p,
                                      Halffacet_const_handle f) {
    return (locate_point_in_halffacet( p, f) == CGAL::ON_BOUNDED_SIDE);
  }

  static Bounded_side locate_point_in_halffacet(const Point_3& p,
                                                Halffacet_const_handle f) {
    CGAL_NEF_TRACEN("locate point in halffacet " << p << ", " << f->plane());
    typedef Project_shalfedge_point
      < SHalfedge, const Point_3> Project;
    typedef Circulator_project
      < SHalfedge_around_facet_const_circulator, Project,
      const Point_3&, const Point_3*> Circulator;
    typedef Container_from_circulator<Circulator> Container;

    Plane_3 h(f->plane());
    CGAL_assertion(h.has_on(p));
    Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
    Bounded_side outer_bound_pos(CGAL::ON_BOUNDARY);
    if (fc.is_shalfedge() ) {
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hfc(se);
      Circulator c(hfc);
      Container ct(c);
      CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
      outer_bound_pos = bounded_side_3(ct.begin(), ct.end(), p, h);
    }
    else
      CGAL_error_msg( "is facet first cycle a SHalfloop?");
    if( outer_bound_pos != CGAL::ON_BOUNDED_SIDE )
      return outer_bound_pos;
    /* The point p is not in the relative interior of the outer face cycle
       so it is not necessary to know the position of p with respect to the
       inner face cycles */
    Halffacet_cycle_const_iterator fe = f->facet_cycles_end();
    ++fc;
    if( fc == fe )
      return outer_bound_pos;
    Bounded_side inner_bound_pos(CGAL::ON_BOUNDARY);
    CGAL_For_all(fc, fe) {
      if (fc.is_shalfloop() ) {
        SHalfloop_const_handle l(fc);
        if(l->incident_sface()->center_vertex()->point() == p )
          inner_bound_pos = CGAL::ON_BOUNDARY;
        else
          inner_bound_pos = CGAL::ON_UNBOUNDED_SIDE;
      }
      else if (fc.is_shalfedge() ) {
        SHalfedge_const_handle se(fc);
        SHalfedge_around_facet_const_circulator hfc(se);
        Circulator c(hfc);
        Container ct(c);
        CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
        inner_bound_pos = bounded_side_3( ct.begin(), ct.end(),
                                          p, h.opposite());
      }
      else
        CGAL_error_msg( "Damn wrong handle.");
      if( inner_bound_pos != CGAL::ON_UNBOUNDED_SIDE )
        return opposite(inner_bound_pos);
      /* At this point the point p belongs to relative interior of the facet's
         outer cycle, and its position is completely known when it belongs
         to the clousure of any inner cycle */
    }
    return CGAL::ON_BOUNDED_SIDE;
  }

}; // SNC_intersection

} //namespace CGAL

#endif //CGAL_SNC_INTERSECTION_H
