// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de> 
#ifndef CGAL_SNC_INTERSECTION_H
#define CGAL_SNC_INTERSECTION_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>

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
class SNC_intersection : public SNC_const_decorator<SNC_structure_> {
  // TODO: granados: is it really necessary to inherit from the decorator?

  typedef SNC_structure_                     SNC_structure;
  typedef SNC_intersection<SNC_structure>    Self;
  typedef SNC_const_decorator<SNC_structure> Base;
  //  typedef SNC_const_decorator<SNC_structure> SNC_const_decorator;

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

#ifdef CGAL_NEF3_FACET_WITH_BOX
  typedef typename SNC_structure::Partial_facet Partial_facet;
#endif

  typedef typename SNC_structure::Point_3        Point_3;
  typedef typename SNC_structure::Vector_3       Vector_3;
  typedef typename SNC_structure::Segment_3      Segment_3;
  typedef typename SNC_structure::Line_3         Line_3;
  typedef typename SNC_structure::Ray_3          Ray_3;
  typedef typename SNC_structure::Plane_3        Plane_3;
  typedef typename SNC_structure::Triangle_3     Triangle_3;

 public:

  SNC_intersection() : Base() {}
  SNC_intersection(const SNC_structure& W) : Base(W) {}

  bool does_contain_internally(const Segment_3& s, const Point_3& p) const {
    if(!are_strictly_ordered_along_line (s.source(), p, s.target()))
      return false;
    if(!s.supporting_line().has_on(p))
      return false;
    return true;
  }

  bool does_contain_internally( Halffacet_const_handle f, 
				const Point_3& p,
				bool check_has_on = true) const {
    if(check_has_on && !f->plane().has_on(p))
      return false;
    return (locate_point_in_halffacet( p, f) == CGAL::ON_BOUNDED_SIDE); 
  }

#ifdef CGAL_NEF3_FACET_WITH_BOX
  bool does_contain_internally( Partial_facet& pf, 
				const Point_3& p) const {
    CGAL_NEF_TRACEN("does point lie in partial facet" << p);
    //    pf.debug();
    if( !pf.f->plane().has_on(p))
      return false;
    return (locate_point_in_halffacet( p, pf) == CGAL::ON_BOUNDED_SIDE); 
  }
#endif

  bool does_contain_on_boundary( Halffacet_const_handle f, const Point_3& p) const {
    typedef Project_shalfedge_point
      < SHalfedge, const Point_3> Project;
    typedef Circulator_project
      < SHalfedge_around_facet_const_circulator, Project, 
      const Point_3&, const Point_3*> Circulator;
    Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
    CGAL_assertion(fc.is_shalfedge());
    if (fc.is_shalfedge() ) {
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hfc(se);
      Circulator c(hfc), cp(c), cend(c);
      do {
	c++;
	CGAL_NEF_TRACEN("contained on edge "<<Segment_3( *c, *cp)<<"? "<<
	       Segment_3( *c, *cp).has_on(p));
	if( Segment_3( *c, *cp).has_on(p))
	  return true;
	cp++;
      }
      while( c != cend); 
    } 
    Halffacet_cycle_const_iterator fe = f->facet_cycles_end();
    ++fc;
    CGAL_For_all(fc, fe) {
      if (fc.is_shalfloop() ) { 
	SHalfloop_const_handle l(fc);
	CGAL_NEF_TRACEN("isolated point on "<<l->incident_sface()->center_vertex()->point()<<"? ");
	if( l->incident_sface()->center_vertex()->point() == p)
	  return true;
      } 
      else if (fc.is_shalfedge() ) {
	SHalfedge_const_handle se(fc);
	SHalfedge_around_facet_const_circulator hfc(se);
	Circulator c(hfc), cp(c), cend(c);
	do {
	  c++;
	  CGAL_NEF_TRACEN("contained on edge "<<Segment_3( *c, *cp)<<"? "<<
		 Segment_3( *c, *cp).has_on(p));
	  if( Segment_3( *c, *cp).has_on(p))
	    return true;
	  cp++;
	} 
	while( c != cend);
      }
      else 
	CGAL_error_msg( "Damn wrong handle.");
    }
    return false;
  }
  
#ifdef LINE3_LINE3_INTERSECTION
  
  bool does_intersect_internally( const Segment_3& s1, 
				  const Segment_3& s2, 
				  Point_3& p) const  {
    CGAL_NEF_TRACEN("does intersect internally with  LINE3_LINE3_INTERSECTION");
    if ( s1.is_degenerate() || s2.is_degenerate())
      /* the segment is degenerate so there is not internal intersection */
      return false;
    if ( s1.has_on(s2.source()) || s1.has_on(s2.target()) ||
	 s2.has_on(s1.source()) || s2.has_on(s1.target()))
      /* the segments does intersect at one endpoint */
      return false;
    Object o = intersection(Line_3(ray), Line_3(s)); 
    if ( !CGAL::assign(p, o))
      return false;
    return( does_contain_internally( s, p));
  }

#else // LINE3_LINE3_INTERSECTION

  bool does_intersect_internally( const Segment_3& s1, 
				  const Segment_3& s2, 
				  Point_3& p) const {
    if(s2.has_on(s1.target()))
      return false;
    Ray_3 r(s1.source(), s1.target());
    if(!does_intersect_internally(r, s2, p))
      return false;
    Plane_3 pl(s1.target(), r.to_vector());
    return (pl.oriented_side(p) == CGAL::NEGATIVE);
  }

  bool does_intersect_internally( const Ray_3& s1, 
				  const Segment_3& s2, 
				  Point_3& p) const {
    CGAL_NEF_TRACEN("does intersect internally without  LINE3_LINE3_INTERSECTION");    
    CGAL_assertion(!s1.is_degenerate());
    CGAL_assertion(!s2.is_degenerate());
    if ( orientation( s1.source(), s1.point(1), s2.source(), s2.target()) 
	 != COPLANAR)
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
      ws1(cross_product( vt, vs1)); // , ws2(cross_product( vt, vs2));
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
    
#endif // LINE3_LINE3_INTERSECTION

  bool does_intersect( const Ray_3& r, const Triangle_3& tr,
		       Point_3& ip) const {
    // Intersection between an open ray and
    // a closed 2d-triangular region in the space
    CGAL_NEF_TRACEN("-> Intersection triangle - ray");
    CGAL_NEF_TRACEN(" -> Ray: "<<r);
    CGAL_NEF_TRACEN(" -> Triangle: "<<tr);
    CGAL_assertion( !r.is_degenerate());
    Plane_3 h( tr.supporting_plane());
    CGAL_assertion( !h.is_degenerate());
    if( h.has_on( r.source()))
      return false;
    Object o = intersection( h, r);
    if( !CGAL::assign( ip, o))
      return false;
    CGAL_NEF_TRACEN(" -> intersection point: "<<ip);
    return tr.has_on(ip);
  }

  bool does_intersect( const Segment_3& s, const Triangle_3& tr,
		       Point_3& ip) const {
    // Intersection between a open segment and
    // a closed 2d-triangular region in the space
    CGAL_NEF_TRACEN("-> Intersection triangle - segment");
    CGAL_NEF_TRACEN(" -> Segment: "<<s);
    CGAL_NEF_TRACEN(" -> Triangle: "<<tr);
    CGAL_assertion( !s.is_degenerate());
    Plane_3 h( tr.supporting_plane());
    CGAL_assertion( !h.is_degenerate());
    if( h.has_on( s.source()) || h.has_on( s.target()))
      return false;
    Object o = intersection( h, s);
    if( !CGAL::assign( ip, o))
      return false;
    CGAL_NEF_TRACEN(" -> intersection point: "<<ip);
    return tr.has_on(ip);
  }

  bool does_intersect_internally( const Ray_3& ray,
				  Halffacet_const_handle f,
				  Point_3& p,
				  bool checkHasOn = true) const { 
    CGAL_NEF_TRACEN("-> Intersection facet - ray");
    Plane_3 h( f->plane());
    CGAL_NEF_TRACEN("-> facet's plane: " << h);
    CGAL_NEF_TRACEN("-> a point on the plane: " << h.point());
    CGAL_NEF_TRACEN("-> ray: " << ray);
    CGAL_assertion(!ray.is_degenerate());
    if(checkHasOn) {
      if(h.has_on(ray.source()))
	return false;
    } else
      CGAL_assertion(!h.has_on(ray.source()));
    Object o = intersection( h, ray);
    if( !CGAL::assign( p, o))
      return false;
    CGAL_NEF_TRACEN( "-> intersection point: " << p );
    CGAL_NEF_TRACEN( "-> point in facet interior? "<<does_contain_internally( f, p));
    return does_contain_internally( f, p, false);
  }

#ifdef CGAL_NEF3_FACET_WITH_BOX
  bool does_intersect_internally( const Ray_3& ray,
				  Partial_facet pf,
				  Point_3& p) const { 
    CGAL_NEF_TRACEN("-> Intersection facet - ray");
    Plane_3 h( pf.f->plane());
    CGAL_NEF_TRACEN("-> facet's plane: " << h);
    CGAL_NEF_TRACEN("-> a point on the plane: " << h.point());
    CGAL_NEF_TRACEN("-> ray: " << ray);
    CGAL_assertion(!ray.is_degenerate());
    if( h.has_on( ray.source()))
      /* no possible internal intersection */
	return false;
    Object o = intersection( h, ray);
    if( !CGAL::assign( p, o))
      return false;
    CGAL_NEF_TRACEN( "-> intersection point: " << p );
    //    CGAL_NEF_TRACEN( "-> point in facet interior? "<<does_contain_internally( f, p));
    return does_contain_internally( pf, p, false);
  }
#endif

  bool does_intersect_internally( const Segment_3& seg,
				  Halffacet_const_handle f,
				  Point_3& p) const { 
    CGAL_NEF_TRACEN("-> Intersection facet - segment");
    Plane_3 h( f->plane());
    CGAL_NEF_TRACEN("-> facet's plane: " << h);
    CGAL_NEF_TRACEN("-> a point on the plane: " << h.point());
    CGAL_NEF_TRACEN("-> segment: " << seg);
    CGAL_assertion(!seg.is_degenerate());
    if( h.has_on( seg.source()) || h.has_on(seg.target()))
      /* no possible internal intersection */
      return false;
    return does_intersect(seg, f, p);
  }

  bool does_intersect(const Segment_3& seg,
		      Halffacet_const_handle f,
		      Point_3& p) const {
    Plane_3 h( f->plane());
    Object o = intersection( h, seg);
    if( !CGAL::assign( p, o))
      return false;
    CGAL_NEF_TRACEN( "-> intersection point: " << p );
    CGAL_NEF_TRACEN( "-> point in facet interior? "<<does_contain_internally( f, p));
    return( does_contain_internally( f, p, false));
  }

#ifdef CGAL_NEF3_FACET_WITH_BOX
  bool does_intersect_internally( const Segment_3& seg,
				  Partial_facet pf,
				  Point_3& p) const { 
    CGAL_NEF_TRACEN("-> Intersection partial facet - segment");
    Plane_3 h( pf.f->plane());
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
    //    CGAL_NEF_TRACEN( "-> point in facet interior? "<<does_contain_internally( f, p));
    return( does_contain_internally( pf, p, false));
  }
#endif

  Bounded_side locate_point_in_halffacet( const Point_3& p, 
					  Halffacet_const_handle f) const {
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
       so it is not necesary to know the possition of p with respect to the 
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
	 outer cycle, and its possition is completely known when it belongs
	 to the clousure of any inner cycle */
    }
    return CGAL::ON_BOUNDED_SIDE;
  }

#ifdef CGAL_NEF3_FACET_WITH_BOX
  Bounded_side locate_point_in_halffacet( const Point_3& p, 
					 Partial_facet pf) const {

    if(p.x() < pf.f->b.min_coord(0) || p.x() > pf.f->b.max_coord(0) ||
       p.y() < pf.f->b.min_coord(1) || p.y() > pf.f->b.max_coord(1) ||
       p.z() < pf.f->b.min_coord(2) || p.z() > pf.f->b.max_coord(2))
      return CGAL::ON_UNBOUNDED_SIDE;
    
    typedef Project_shalfedge_point
      < SHalfedge, Point_3> Project;
    typedef Circulator_project
      < SHalfedge_around_facet_const_circulator, Project, 
      const Point_3&, const Point_3*> Circulator;
    typedef Container_from_circulator<Circulator> Container;

    typedef typename Partial_facet::Outer_cycle_iterator  Outer_cycle_iterator;
    typedef typename Partial_facet::Inner_cycle_iterator  Inner_cycle_iterator;
    typedef typename Partial_facet::Isolated_vertex_iterator Isolated_vertex_iterator;

    Plane_3 h(pf.f->plane());
    CGAL_assertion(h.has_on(p));
    
    Bounded_side outer_bound_pos(CGAL::ON_BOUNDED_SIDE);

    Outer_cycle_iterator oc = pf.outer_cycles_begin();
    while(oc != pf.outer_cycles_end() && 
	  outer_bound_pos == CGAL::ON_BOUNDED_SIDE) {
      if(oc->first == oc->second) {
	SHalfedge_around_facet_const_circulator hfc(oc->first);
	Circulator c(hfc);
	Container ct(c);
	CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
	outer_bound_pos = bounded_side_3(ct.begin(), ct.end(), p, h);	
      } else {
	outer_bound_pos = bounded_side_3(Circulator(SHalfedge_around_facet_const_circulator(oc->first)), 
					 Circulator(SHalfedge_around_facet_const_circulator(oc->second)), p, h);
      }
      ++oc;
    }
    if(outer_bound_pos != CGAL::ON_BOUNDED_SIDE )
      return outer_bound_pos;    

    Bounded_side inner_bound_pos(CGAL::ON_UNBOUNDED_SIDE);

    Inner_cycle_iterator ic = pf.inner_cycles_begin();
    while(ic != pf.inner_cycles_end() && 
	  inner_bound_pos == CGAL::ON_UNBOUNDED_SIDE) {
      if(ic->first == ic->second) {
	SHalfedge_around_facet_const_circulator hfc(ic->first);
	Circulator c(hfc);
	Container ct(c);
	CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
	inner_bound_pos = bounded_side_3(ct.begin(), ct.end(), p, h);	
      } else {
	inner_bound_pos = bounded_side_3(Circulator(SHalfedge_around_facet_const_circulator(ic->first)), 
					 Circulator(SHalfedge_around_facet_const_circulator(ic->second)), p, h);
      }
      ++ic;
    }
    if(inner_bound_pos != CGAL::ON_UNBOUNDED_SIDE )
      return opposite(inner_bound_pos);
    
    Isolated_vertex_iterator iv = pf.isolated_vertices_begin();
    while(iv != pf.isolated_vertices_end()) {
      if(*iv == p)
	return CGAL::ON_BOUNDARY;
      ++iv;
    }
   
    return CGAL::ON_BOUNDED_SIDE;
  }
#endif
}; // SNC_intersection

} //namespace CGAL

#endif //CGAL_SNC_INTERSECTION_H
