// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_ray_shoter.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_ray_shoter.h       Ray shoting and  point location on SNC structures
// ============================================================================
#ifndef CGAL_SNC_RAY_SHOTER_H
#define CGAL_SNC_RAY_SHOTER_H

#include <CGAL/basic.h>
#include <CGAL/functional.h> 
#include <CGAL/function_objects.h> 
#include <CGAL/Circulator_project.h>
#include <CGAL/Nef_3/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#include <map>
#include <list>
#undef _DEBUG
#define _DEBUG 0 
//37
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template < class Node, class Object, class DClass>
struct Project_halfedge_point {
  typedef Node         argument_type;
  typedef Object       result_type;
  Object& operator()( Node& x) const   { 
    DClass D;
    return D.point(D.source(x));
    /* a Point_3& reference must be returned by D.point() */
  }
  const Object& operator()( const Node& x) const   { 
    DClass D;
    return D.point(D.source(x)); 
    /* a Point_3& reference must be returned by D.point() */
  }
};

// ----------------------------------------------------------------------------
// SNC_ray_shoting
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_ray_shoting}{SNC}{ray shot functionality}{O}}*/

template <typename SNC_structure_>
class SNC_ray_shoter : public SNC_decorator<SNC_structure_>
{ 
  typedef SNC_structure_ SNC_structure;

protected:
  typedef SNC_ray_shoter<SNC_structure>           Self;
  typedef SNC_decorator<SNC_structure>            Base;

public:
  typedef typename SNC_structure_::Kernel         Kernel;
  typedef SNC_decorator<SNC_structure>            SNC_decorator;
  typedef SNC_SM_decorator<SNC_structure>         SM_decorator;
  typedef SNC_SM_point_locator<SNC_structure>     SM_point_locator;
  typedef SNC_SM_const_decorator<SNC_structure>   SM_const_decorator;

  #define USING(t) typedef typename SNC_structure::t t
  USING(Vertex);
  USING(Halfedge);
  USING(Halffacet);
  USING(Volume);
  
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Volume_iterator);

  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Volume_handle);

  USING(Vertex_const_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_handle);

  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SFace_iterator);
  USING(SHalfloop_iterator);

  USING(SVertex);
  USING(SHalfedge);
  USING(SFace);
  USING(SHalfloop);

  USING(SVertex_handle);
  USING(SHalfedge_handle);
  USING(SFace_handle);
  USING(SHalfloop_handle);

  USING(SVertex_const_handle); 
  USING(SHalfedge_const_handle); 
  USING(SHalfloop_const_handle); 
  USING(SFace_const_handle); 

  USING(Object_handle);
  USING(SObject_handle);

  USING(SHalfedge_around_facet_const_circulator);
  USING(SHalfedge_around_facet_circulator);
  USING(SFace_cycle_iterator);
  USING(SFace_cycle_const_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Halffacet_cycle_const_iterator);
  USING(Shell_entry_iterator);
  USING(Shell_entry_const_iterator);

  USING(Point_3);
  USING(Vector_3);
  USING(Direction_3);
  USING(Segment_3);
  USING(Line_3);
  USING(Plane_3);

  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Sphere_direction);

  USING(Mark);
  #undef USING

  #define DECUSING(t) typedef typename SM_decorator::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
  DECUSING(SHalfedge_around_svertex_circulator);
  #undef DECUSING

  typedef void* GenPtr;

  SNC_ray_shoter(SNC_structure& W) : Base(W) {}
  /*{\Mcreate makes |\Mvar| a ray shoter on |W|.}*/

  Object_handle shot( Segment_3& ray) const 
     /*{\Mop returns the nearest object hit by a ray |ray|. }*/ {
    TRACEN( "Shoting ray " << ray);
    Object_handle o;
    Vertex_handle v;
    CGAL_nef3_forall_vertices( v, *sncp()) {
      if ( ray.source() != point(v) && ray.has_on(point(v)) ) {
	TRACEN("ray hit vertex case "<<point(v));
	shorten( ray, point(v));
	o = Object_handle(v);
      }
    }
    Halfedge_handle e;
    CGAL_nef3_forall_edges( e, *sncp()) {
      Point_3 q;
      if ( does_contain_internally( segment(e), ray.target())) {
	TRACEN("ray target on edge "<<segment(e));
	o = Object_handle(e);
      }
      else if( does_intersect_internally( ray, e, q) ) { 
	shorten( ray, q); 
	o = Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_nef3_forall_halffacets( f, *sncp()) {
      Point_3 q;
      if ( does_contain_internally( f, ray.target())) {
	TRACEN("ray target on facet "<<plane(f));
	o = Object_handle(f);
      }
      else if( does_intersect_internally( ray, f, q) ) {
	TRACEN("ray hit facet "<<plane(f)<<" on "<<q);
	shorten( ray, q); 
	o = Object_handle(f);
      }
    }
    return o;
  }

  Object_handle locate( const Point_3& p) const
    /*{\Mop returns the lowest dimension object on an SNC structure
      which contais |p| in its interior. }*/{
    TRACEN( "Point locator for " << p);
    Vertex_handle v;
    CGAL_nef3_forall_vertices( v, *sncp()) {
      if ( p == point(v)) {
	TRACEN("on vertex.");
	return Object_handle(v);
      }
    }
    Halfedge_handle e;
    CGAL_nef3_forall_edges( e, *sncp()) {
      if ( does_contain_internally( segment(e), p) ) {
	TRACEN("on edge.");
	return Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_nef3_forall_halffacets( f, *sncp()) {
      if ( does_contain_internally( f, p) ) {
	TRACEN("on facet.");
	return Object_handle(f);
      }
    }
    Volume_handle c;
    /* lets be |s| be the segment that connects |p| to any fixed vertex |va| */
    Vertex_handle va = --(sncp()->vertices_end()); 
    Segment_3 s( p, point(va));
    /* prune |s| by |o| if |o| intersects |s| in its relative interior */
    Object_handle o = shot(s);
    /* determine the volume that contains |s| from the last pruning object */
    if( assign( v, o))
      c = volume(get_visible_facet( v, s));
    else if( assign( e, o))
      c = volume(get_visible_facet( e, s));
    else if( assign( f, o))
      c = volume(get_visible_facet( f, s));
    else CGAL_nef3_assertion_msg(0, "where is our point, eh?");
    return Object_handle(c);
  }

  void shorten(Segment_3& s, const Point_3& p) const { 
    s = Segment_3( s.source(), p); 
    TRACEN("shoted ray "<<s);
  }
#undef _DEBUG
#define _DEBUG 2
  bool does_contain_internally(const Segment_3& s, const Point_3& p) const {
    if(!s.has_on(p))
      return false;
    Comparison_result r1 = compare_xyz(s.source(),p); 
    Comparison_result r2 = compare_xyz(s.target(),p); 
    return (r1 == opposite(r2));
  }

  bool does_intersect_internally( const Segment_3& s,
				  const Halfedge_handle e,
				  Point_3& p) const {
    return does_intersect_internally( s, segment(e), p);
  }

#ifdef LINE3_LINE3_INTERSECTION

  bool does_intersect_internally( const Segment_3& s1, 
				  const Segment_3& s2, 
				  Point_3& p) const  {
    if ( s1.is_degenerate() || s2.is_degenerate())
      /* the segment is degenerate so there is not internal intersection */
      return false;
    if ( s1.has_on(s2.source()) || s1.has_on(s2.target()) ||
	 s2.has_on(s1.source()) || s2.has_on(s1.target()))
      /* the segments does intersect at one endpoint */
      return false;
    Object o = intersection(Line_3(ray), Line_3(s)); 
    if ( !assign(p, o))
      return false;
    return( does_contain_internally( s, p));
  }

#else // LINE3_LINE3_INTERSECTION

  bool does_intersect_internally( const Segment_3& s1, 
				  const Segment_3& s2, 
				  Point_3& p) const {
    if ( s1.is_degenerate() || s2.is_degenerate())
      /* the segment is degenerate so there is not internal intersection */
      return false;
    if ( s1.has_on(s2.source()) || s1.has_on(s2.target()) ||
	 s2.has_on(s1.source()) || s2.has_on(s1.target()))
      /* the segments does intersect at one endpoint */
      return false;
    if ( orientation( s1.source(), s1.target(), s2.source(), s2.target()) 
	 != COPLANAR)
      /* the segments doesn't define a plane */
      return false;
    if ( collinear( s1.source(), s1.target(), s2.source()) &&
	 collinear( s1.source(), s1.target(), s2.target()) )
      /* the segments are collinear */
      return false;
    Line_3 ls1(s1), ls2(s2);
    if ( ls1.direction() ==  ls2.direction() ||
	 ls1.direction() == -ls2.direction() )
      /* the segments are parallel */
      return false;
    Oriented_side os1, os2;
    Vector_3 vs1(s1.direction()), vs2(s2.direction()), 
      vt(cross_product( vs1, vs2)), 
      ws1(cross_product( vt, vs1)), ws2(cross_product( vt, vs2));
    Plane_3 hs1( s1.source(), ws1);
    /* hs1 is a plane which contains line(s1) and is perpendicular to the
       plane defined by the s1 and s2 */
    os1 = hs1.oriented_side(s2.source());
    os2 = hs1.oriented_side(s2.target());
    if(os1 != opposite(os2))
      return false;
    Plane_3 hs2( s2.source(), ws2);
    /* hs is a plane which contains line(s2) and is perpendicular to the
       plane defined by the s1 and s */
    os1 = hs2.oriented_side(s1.source());
    os2 = hs2.oriented_side(s1.target());
    if(os1 != opposite(os2))
      return false;
    Object o = intersection(hs1, ls2);
    CGAL_nef3_assertion(assign( p, o));
    /* since line(s1) and line(s2) are not parallel they intersects in only
       one point */
    assign( p ,o);
    return( does_contain_internally( s2, p));
  }
#endif // LINE3_LINE3_INTERSECTION

  bool does_contain_internally( const Halffacet_handle f, 
				const Point_3& p) const {
    if( !plane(f).has_on(p))
      return false;
    return (locate_point_in_halffacet( p, f) == CGAL::ON_BOUNDED_SIDE);
  }

  bool does_intersect_internally( const Segment_3& seg,
				  const Halffacet_handle f,
				  Point_3& p) const { 
    TRACEN("-> Intersection face - ray");
    Plane_3 h( plane(f));
    TRACEN("-> facet plane " << h);
    TRACEN("-> a point on " << h.point());
    TRACEN("-> seg segment " << seg);
    CGAL_nef3_assertion( !h.is_degenerate());
    if( seg.is_degenerate())
      /* no possible internal intersection */
      return false;
    if( h.has_on( seg.source()) || h.has_on( seg.target()))
      /* no possible internal intersection */
	return false;
#ifdef REDUNDANT_CODE
    /* This optimization might be inside of |intersection()| code */
    Oriented_side os1 = h.oriented_side(seg.source());
    Oriented_side os2 = h.oriented_side(seg.target());
    TRACEN( "-> endpoint plane side " << os1 << " " << os2);
    CGAL_nef3_assertion( h.has_on(p));
    CGAL_nef3_assertion( seg.has_on(p));
    if (os1 == os2)
      return false;
#endif //REDUNDANT_CODE
    Object o = intersection( h, seg);
    Segment_3 s;
    if ( assign( s, o) ) {
      CGAL_nef3_assertion( s == seg );
      TRACEN( "-> seg belongs to facet's plane." << p );
      return false;
    }
    else if( !assign( p, o))
      return false;
    TRACEN( "-> intersection point " << p );
    TRACEN( "-> point in facet? "<<locate_point_in_halffacet(p, f));
    return does_contain_internally( f, p);
  }

  Bounded_side locate_point_in_halffacet( const Point_3& p, 
					  const Halffacet_handle f) const {
    typedef Project_halfedge_point
      < SHalfedge, const Point_3, SNC_decorator> Project;
    typedef Circulator_project
      < SHalfedge_around_facet_circulator, Project, 
      const Point_3&, const Point_3*> Circulator;
    typedef Container_from_circulator<Circulator> Container;

    Plane_3 h(plane(f));
    CGAL_nef3_assertion(h.has_on(p));
    Halffacet_cycle_iterator fc = f->facet_cycles_begin();
    SHalfedge_handle se;
    Bounded_side outer_bound_pos;
    if ( assign(se,fc) ) {
      SHalfedge_around_facet_circulator hfc(se);
      Circulator c(hfc);
      Container ct(c);
      CGAL_nef3_assertion( !is_empty_range(ct.begin(), ct.end()));
      outer_bound_pos = bounded_side_3(ct.begin(), ct.end(), p, h);
    } 
    else 
      CGAL_nef3_assertion_msg(0, "is facet first cycle a SHalfloop?");
    if( outer_bound_pos != CGAL::ON_BOUNDED_SIDE )
      return outer_bound_pos;
    /* The point p is not in the relative interior of the outer face cycle
       so it is not necesary to know the possition of p with respect to the 
       inner face cycles */
    Halffacet_cycle_iterator fe = f->facet_cycles_end();
    ++fc;
    if( fc == fe )
      return outer_bound_pos;
    Bounded_side inner_bound_pos;
    CGAL_For_all(fc, fe) {
      SHalfloop_handle l;
      if ( assign(l,fc) ) { 
        if( point(vertex(sface(l))) == p )
	  inner_bound_pos = CGAL::ON_BOUNDARY;
	else
	  inner_bound_pos = CGAL::ON_UNBOUNDED_SIDE;
      } 
      else if ( assign(se,fc) ) {
	SHalfedge_around_facet_circulator hfc(se);
	Circulator c(hfc);
	Container ct(c);
	CGAL_nef3_assertion( !is_empty_range(ct.begin(), ct.end()));
        inner_bound_pos = bounded_side_3( ct.begin(), ct.end(), 
					  p, h.opposite());
      } 
      else 
	CGAL_nef3_assertion_msg(0, "Damn wrong handle.");
      if( inner_bound_pos != CGAL::ON_UNBOUNDED_SIDE )
	return opposite(inner_bound_pos);
      /* At this point the point p belongs to relative interior of the facet's
	 outer cycle, and its possition is completely known when it belongs
	 to the clousure of any inner cycle */
    }
    return CGAL::ON_BOUNDED_SIDE;
  }
}; // SNC_ray_shot

CGAL_END_NAMESPACE

#endif //CGAL_SNC_RAY_SHOTER_H

