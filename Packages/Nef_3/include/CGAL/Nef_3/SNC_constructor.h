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
// file          : include/CGAL/Nef_3/SNC_constructor.h
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
// SNC_constructor.h       construction of basic SNCs and global construction
// ============================================================================
#ifndef CGAL_SNC_CONSTRUCTOR_H
#define CGAL_SNC_CONSTRUCTOR_H

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
#define _DEBUG 43
#include <CGAL/Nef_3/debug.h>

#define IMMN 12345

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
  
template <typename Point, typename Edge, class Decorator>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge,Decorator> Self;
  Point p; int i; Edge e;
  Decorator& D;
  Halfedge_key(Point pi, int ii, Edge ei, Decorator& Di ) : 
    p(pi), i(ii), e(ei), D(Di) {}
  Halfedge_key(const Self& k) : p(k.p), i(k.i), e(k.e), D(k.D) {}
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge,Decorator> Key;
  typedef typename Point::R R;
  typedef typename R::Vector_3 Vector;
  typedef typename R::Direction_3 Direction;
  bool operator()( const Key& k1, const Key& k2) const { 
    if ( k1.p == k2.p ) 
      return (k1.i < k2.i);
    /* previous code: 
       else return CGAL::lexicographically_xyz_smaller(k1.p,k2.p); */
    Direction l(Vector(CGAL::ORIGIN, k1.D.tmp_point(k1.e)));
    if( k1.i < 0) l = -l;
    return (Direction( k2.p - k1.p) == l); 
  }
};

template <typename Point, typename Edge, class Decorator>
std::ostream& operator<<(std::ostream& os, 
                         const Halfedge_key<Point,Edge,Decorator>& k )
{ os << k.p << " " << k.i; return os; }

template <typename R>
int sign_of(const CGAL::Plane_3<R>& h)
{ if ( h.a() != 0 ) return CGAL_NTS sign(h.a());
  if ( h.b() != 0 ) return CGAL_NTS sign(h.b());
  return CGAL_NTS sign(h.c());
}

template <typename R>
CGAL::Plane_3<R> normalized(CGAL::Plane_3<R>& h)
{ typedef typename R::RT RT;
  RT a(h.a()),b(h.b()),c(h.c()),d(h.d());
  RT x = ( a != 0 ? a : 1);
  x = ( b != 0 ? gcd(x,b) : x );
  x = ( c != 0 ? gcd(x,c) : x );
  x = ( d != 0 ? gcd(x,d) : x );
  return CGAL::Plane_3<R>(a/x,b/x,c/x,d/x); 
}

struct Plane_lt {
  template <typename R>
  bool operator()(const CGAL::Plane_3<R>& h1,
                  const CGAL::Plane_3<R>& h2) const
  { typedef typename R::RT RT;
    RT diff = h1.a()-h2.a();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.b()-h2.b();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.c()-h2.c();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.d()-h2.d(); return CGAL_NTS sign(diff) < 0;
  }
};

// ----------------------------------------------------------------------------
// SNC_constructor 
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_constructor}{SNC}{overlay functionality}{O}}*/

template <typename SNC_structure_>
class SNC_constructor : public SNC_decorator<SNC_structure_>
{ 
public:
  typedef SNC_structure_ SNC_structure;
  typedef typename SNC_structure_::Sphere_kernel  Sphere_kernel;
  typedef typename SNC_structure_::Kernel         Kernel;
  typedef SNC_constructor<SNC_structure>          Self;
  typedef SNC_decorator<SNC_structure>            Base;
  typedef SNC_decorator<SNC_structure>            SNC_decorator;
  typedef SNC_FM_decorator<SNC_structure>         FM_decorator;
  typedef SNC_SM_decorator<SNC_structure>         SM_decorator;
  typedef SNC_SM_overlayer<SNC_structure>         SM_overlayer;
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

  typedef CGAL::Unique_hash_map<SFace_const_handle,unsigned int>  
                                                         Shell_number_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> Shell_closed_hash;

  struct Shell_explorer {
    const SNC_decorator& D;
    Shell_number_hash&  Shell;
    Shell_closed_hash& Closed;
    SFace_visited_hash& Done;
    Vertex_handle v_min;
    int n;

    Shell_explorer(const SNC_decorator& Di, Shell_number_hash& Si, 
                   Shell_closed_hash& Sc, SFace_visited_hash& Vi) 
      : D(Di), Shell(Si), Closed(Sc), Done(Vi), n(0) {}

    void visit(SFace_handle h) { 
      TRACEN("visit sf "<<D.point(D.vertex(h)));
      Shell[h]=n;
      Done[h]=true;
    }
    void visit(Vertex_handle h) { 
      TRACEN("visit v  "<<D.point(h));
      if ( CGAL::lexicographically_xyz_smaller(
           D.point(h),D.point(v_min)) ) 
	v_min = h; 
    }
    void visit(Halfedge_handle h) { 
      TRACEN("visit he "<<D.point(D.source(h)));
      SFace_handle sf = D.sface(h);
      if( Closed[sf] ) {
	SM_decorator SD(D.vertex(h));
	if( SD.first_out_edge(h) == SD.last_out_edge(h) )
	  Closed[sf] = false;
      }
    }
    void visit(Halffacet_handle h) { /* do nothing */ }

    Vertex_handle& minimal_vertex() { return v_min; }

    void increment_shell_number() { 
      TRACEN("leaving shell "<<n);
      ++n; 
    }
  };

  SNC_constructor(SNC_structure& W) : Base(W) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/


  Vertex_handle create_box_corner(int x, int y, int z,
                                  bool boundary=true) const; 
  /*{\Mop produces the sphere map representing the box corner in
          direction $(x,y,z)$.}*/

  Vertex_handle create_from_facet(Halffacet_handle f,
				  const Point_3& p) const; 
  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $f$. \precond $p$ is part of $f$.}*/

  Vertex_handle create_from_edge(Halfedge_handle e,
				 const Point_3& p) const; 
  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $e$. \precond $p$ is part of $e$.}*/

  void pair_up_halfedges() const;
  /*{\Mop pairs all halfedge stubs to create the edges in 3-space.}*/

  void link_shalfedges_to_facet_cycles() const;
  /*{\Mop creates all non-trivial facet cycles from sedges. 
  \precond |pair_up_halfedges()| was called before.}*/

  void categorize_facet_cycles_and_create_facets() const;
  /*{\Mop collects all facet cycles incident to a facet and creates
  the facets. \precond |link_shalfedges_to_facet_cycles()| was called
  before.}*/

  void create_volumes() const;
  /*{\Mop collects all shells incident to a volume and creates the
  volumes.  \precond |categorize_facet_cycles_and_creating_facets()| was
  called before.}*/

  Volume_handle determine_volume( SFace_handle sf, 
                const std::vector< Vertex_handle>& MinimalVertex, 
                const Shell_number_hash&  Shell ) const {
    Vertex_handle v_min = MinimalVertex[Shell[sf]];
    Halffacet_handle f_below = get_facet_below(v_min);
    if ( f_below == Halffacet_handle())
      return Base(*this).volumes_begin();
    Volume_handle c = volume(f_below);
    if( c != Volume_handle()) {
      TRACE( "Volume " << &*c << " hit ");
      TRACEN("(Shell #" << Shell[adjacent_sface(f_below)] << ")");
      return c;
    }
    SFace_handle sf_below = adjacent_sface(f_below);
    TRACE( "Shell not assigned to a volume hit ");
    TRACEN( "(Inner shell #" << Shell[sf_below] << ")");
    c = determine_volume( sf_below, MinimalVertex, Shell);
    link_as_inner_shell( sf_below, c);
    return c;
  }

  Halffacet_handle get_visible_facet( const Vertex_handle v, 
				      const Segment_3& ray) const {
    Halffacet_handle f_visible;
    CGAL_assertion( ray.target() == point(v));
    Sphere_point sp(ray.source() - point(v));
    TRACEN( "Locating " << sp <<" in " << point(v));
    SM_point_locator L(v);
    SObject_handle o = L.locate(sp);
    SFace_const_handle sf;
    CGAL_assertion( assign( sf, o));
    /* the ray must not belong to the 2-skeleton incident to v */
    assign( sf, o);
    SFace_cycle_const_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) { // TO TEST: is v an isolated vertex?
	TRACEN( "no adjacent facet found.");
	f_visible =  Halffacet_handle();
    }
    else {
      SHalfedge_handle se; 
      SHalfloop_handle sl;
      if ( assign( se, fc)) {
	TRACEN( "adjacent facet found (SEdges cycle).");
	f_visible = facet(twin(se));
      }
      else if ( assign( sl, fc)) {
	TRACEN( "adjacent facet found (SHalfloop cycle).");
	f_visible = facet(twin(sl));
      }
      else 
	CGAL_assertion_msg(0, "Damn, wrong handle.");
    }
    return f_visible;
  }

  typedef typename Kernel::Direction_2 Direction_2;
  typedef typename Kernel::Point_2 Point_2;

  bool strictly_ordered_ccw(const Direction_2& d1, 
			    const Direction_2& d2, 
			    const Direction_2& d3) const
  /*{\Mop returns |true| iff |d2| is in the interior of the
  counterclockwise angular sector between |d1| and |d3|.}*/
  { 
    if ( d1 < d2 )  return ( d2 < d3 )||( d3 <= d1 );
    if ( d1 > d2 )  return ( d2 < d3 )&&( d3 <= d1 );
    return false;
  }

  Halffacet_handle get_visible_facet( const Halfedge_handle e,
				      const Segment_3& ray) const {
    SM_decorator SD;
    if( SD.is_isolated(e))
      return Halffacet_handle();

    Direction_3 ed(segment(e).direction()), rd(-ray.direction());
    Vector_3 ev(ed), rv(ed);
    SHalfedge_around_svertex_circulator sh(SD.first_out_edge(e)), sg(sh);
    Vector_3 h(plane(facet(twin(sh))).orthogonal_vector());
    TRACEN("initial face candidate "<<&*facet(twin(sh)));
    sg++;
    while ( true ) {
      Vector_3 g(plane(facet(twin(sg))).orthogonal_vector());
      if( CGAL_NTS is_positive( cross_product(g, ev) * h)) {
	if( CGAL_NTS is_negative( rv * g))
	  return facet(twin(sh));
	else {
	  sh = sg;
	  h = g;
	  TRACEN("new candidate "<<&*facet(twin(sh)));
	}
      }
      else
	return facet(twin(sh));
    }
    return Halffacet_handle(); // never reached
  }
  
  Halffacet_handle get_visible_facet( const Halffacet_handle f,
				      const Segment_3& ray) const {
    Halffacet_handle f_visible;
    Plane_3 h = plane(f);
    CGAL_assertion( h.has_on(ray.target()));
    CGAL_assertion( !CGAL_NTS is_zero(Vector_3(ray.direction()) * 
				      Vector_3(h.orthogonal_vector())));
    /* is imposible to reach the interior or f using a ray coplanar with f */
    if( CGAL_NTS is_negative(Vector_3(ray.direction()) * 
			     Vector_3(h.orthogonal_vector())))
      f_visible = f;
    else
      f_visible = twin(f);
    CGAL_assertion( CGAL_NTS is_negative( Vector_3(ray.direction()) *
					  Vector_3(plane(f_visible).
						   orthogonal_vector())));
    return f_visible;
  }

  Halffacet_handle get_facet_below( Vertex_handle vi) const {
    Halffacet_handle f_below;
    Segment_3 ray(point(vi), Point_3( 0, 0, -IMMN));
    /* TODO: replace IMMN constant for a real infimaximal number */
    Object_handle o = ray_shot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    if( assign(v, o)) {
	f_below = get_visible_facet(v, ray);
	if( f_below == Halffacet_handle())
	  f_below = get_facet_below(v);
    }
    else if( assign(e, o)) {
	f_below = get_visible_facet(e, ray);
	if( f_below == Halffacet_handle())
	  f_below = get_facet_below(vertex(e));
    }
    else if( assign(f, o)) {
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
    }
    return f_below;
  }
  
  Object_handle ray_shot( Segment_3& ray) const {
    TRACEN( "Shoting ray " << ray);
    Object_handle o;
    Vertex_handle v;
    CGAL_nef3_forall_vertices( v, *sncp()) {
      if ( ray.source() != point(v) && ray.has_on(point(v)) ) {
	TRACEN("ray hit vertex case");
	shorten( ray, point(v));
	o = Object_handle(v);
      }
    }
    Halfedge_handle e;
    CGAL_nef3_forall_edges( e, *sncp()) {
      Point_3 q;
      if ( does_ray_intersect_internally( ray, e, q) ) {
	TRACEN("ray hit edge case");
	shorten( ray, q); 
	o = Object_handle(e);
      }
    }
    Halffacet_handle f;
    CGAL_nef3_forall_halffacets( f, *sncp()) {
      Point_3 q;
      if ( does_ray_intersect_internally( ray, f, q) ) {
	TRACEN("ray hit facet case");
	shorten( ray, q); 
	o = Object_handle(f);
      }
    }
    return o;
  }

  void shorten(Segment_3& s, const Point_3& p) const { 
    s = Segment_3( s.source(), p); 
  }

  bool contains_internally(const Segment_3& s, const Point_3& p) const {
    if(!s.has_on(p))
      return false;
    Comparison_result r1 = compare_xyz(s.source(),p); 
    Comparison_result r2 = compare_xyz(s.target(),p); 
    return (r1 == opposite(r2));
  }

  bool does_ray_intersect_internally( const Segment_3& ray,
				      const Halfedge_handle e,
				      Point_3& p) const {
    return does_ray_intersect_internally( ray, segment(e), p);
  }

#ifdef LINE3_LINE3_INTERSECTION

  bool does_ray_intersect_internally( const Segment_3& ray, 
				      const Segment_3& s, 
				  Point_3& p) const  {
    CGAL_assertion( !ray.is_degenerate());
    if ( s.is_degenerate())
      return false;
    if ( s.has_on(ray.source()) )
      return false;
    Object o = intersection(Line_3(ray), Line_3(s)); 
    if ( !assign(p, o))
      return false;
    return( contains_internally( s, p));
  }

#else // LINE3_LINE3_INTERSECTION

  bool does_ray_intersect_internally( const Segment_3& ray, 
				      const Segment_3& s, 
				      Point_3& p) const {
    CGAL_assertion( !ray.is_degenerate());
    if ( s.is_degenerate())
      return false;
    /* the segment is degenerate so there is not internal intersection */
    if ( s.has_on(ray.source()) )
      return false;
    /* the segment contains the ray source */
    if ( orientation( ray.source(), ray.target(), s.source(), s.target()) 
	 != COPLANAR)
      return false;
    /* the segments doesn't define a plane */
    if ( collinear( ray.source(), ray.target(), s.source()) &&
	 collinear( ray.source(), ray.target(), s.target()) )
      return false;
    /* the segments are collinear */
    Line_3 lray(ray), ls(s);
    if ( lray.direction() ==  ls.direction() ||
	 lray.direction() == -ls.direction() )
      return false;
    /* the segments are parallel */
    Oriented_side os1, os2;
    Vector_3 vray(ray.direction()), vs(s.direction()), 
      vt(cross_product( vray, vs)), 
      wray(cross_product( vt, vray)), ws(cross_product( vt, vs));
    Plane_3 hray( ray.source(), wray);
    /* hray is a plane which contains line(ray) and is perpendicular to the
       plane defined by the ray and s */
    os1 = hray.oriented_side(s.source());
    os2 = hray.oriented_side(s.target());
    if(os1 != opposite(os2))
      return false;
    Plane_3 hs( s.source(), ws);
    /* hs is a plane which contains line(s) and is perpendicular to the
       plane defined by the ray and s */
    os1 = hs.oriented_side(ray.source());
    os2 = hs.oriented_side(ray.target());
    if(os1 != opposite(os2))
      return false;
    Object o = intersection(hray, ls);
    CGAL_assertion(assign( p, o));
    /* since line(ray) and line(s) are not parallel they intersects in only
       one point */
    assign( p ,o);
    return( contains_internally( s, p));
  }
#endif // LINE3_LINE3_INTERSECTION

  bool does_ray_intersect_internally( const Segment_3& ray,
				      const Halffacet_handle f,
				      Point_3& p) const { 
    // TRACEN("-> Intersection face - ray");
    Plane_3 h( plane(f));
    // TRACEN("-> facet plane " << h);
    // TRACEN("-> a point on " << h.point());
    // TRACEN("-> ray segment " << ray);
    CGAL_assertion( !h.is_degenerate());
    CGAL_assertion( !ray.is_degenerate());
    if( h.has_on( ray.source()))
	return false;
    Object o = intersection( h, ray);
    Segment_3 s;
    if ( assign( s, o) ) {
      CGAL_assertion( s == ray );
      // TRACEN( "-> ray belongs to facet's plane." << p );
      return false;
    }
    else if( !assign( p, o))
      return false;
    // TRACEN( "-> intersection point " << p );
    Oriented_side os1 = h.oriented_side(ray.source());
    Oriented_side os2 = h.oriented_side(ray.target());
    // TRACEN( "-> endpoint plane side " << os1 << " " << os2);
    CGAL_assertion( h.has_on(p));
    CGAL_assertion( ray.has_on(p));
    if (os1 == os2)
      return false;
    TRACEN( "-> point in facet? "<<locate_point_in_halffacet(p, f));
    return (locate_point_in_halffacet( p, f) == CGAL::ON_BOUNDED_SIDE);
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
    CGAL_assertion(h.has_on(p));
    Halffacet_cycle_iterator fc = f->facet_cycles_begin();
    SHalfedge_handle se;
    Bounded_side outer_bound_pos;
    if ( assign(se,fc) ) {
      SHalfedge_around_facet_circulator hfc(se);
      Circulator c(hfc);
      Container ct(c);
      CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
      outer_bound_pos = bounded_side_3(ct.begin(), ct.end(), p, h);
    } 
    else 
      CGAL_assertion_msg(0, "is facet first cycle a SHalfloop?");
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
	CGAL_assertion( !is_empty_range(ct.begin(), ct.end()));
        inner_bound_pos = bounded_side_3( ct.begin(), ct.end(), 
					  p, h.opposite());
      } 
      else 
	CGAL_assertion_msg(0, "Damn wrong handle.");
      if( inner_bound_pos != CGAL::ON_UNBOUNDED_SIDE )
	return opposite(inner_bound_pos);
      /* At this point the point p belongs to relative interior of the facet's
	 outer cycle, and its possition is completely known when it belongs
	 to the clousure of any inner cycle */
    }
    return CGAL::ON_BOUNDED_SIDE;
  }

}; // SNC_constructor<SNC>



// ----------------------------------------------------------------------------
// create_box_corner()
// creates the local graph at the corner of a cube in direction (x,y,z)
// boundary specifies if the bounding planes should be included

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_box_corner(int x, int y, int z, bool boundary) const { 
  CGAL_nef3_assertion(x*y*z != 0);
  CGAL_nef3_assertion(x==y==z==1);
  Vertex_handle v = sncp()->new_vertex();
  int R=IMMN; point(v) = Point_3(x*R,y*R,z*R);
  Sphere_point px(-x,0,0), py(0,-y,0), pz(0,0,-z);
  std::list<Sphere_segment> L;
  L.push_back(Sphere_segment(px,py));
  L.push_back(Sphere_segment(py,pz));
  L.push_back(Sphere_segment(px,pz));
  // construction of sphere map from sphere segments:
  SM_overlayer D(v); 
  D.create_from_segments(L.begin(),L.end()); 
  D.simplify();
  SHalfedge_iterator e = D.shalfedges_begin();
  SFace_handle f;
  Sphere_point p1 = D.point(D.source(e));
  Sphere_point p2 = D.point(D.target(e));
  Sphere_point p3 = D.point(D.target(D.next(e)));
  if ( spherical_orientation(p1,p2,p3) > 0 ) f = D.face(e);
  else f = D.face(D.twin(e));
  D.mark(f) = true;
  CGAL_nef3_forall_sedges_of(e,v)
    D.mark(e) = D.mark(D.source(e)) = boundary;
  D.mark_of_halfsphere(-1) = (x<0 && y>0 && z>0);
  D.mark_of_halfsphere(+1) = (x>0 && y>0 && z<0);
  return v;
} 


// ----------------------------------------------------------------------------
// create_from_facet() 
// Creates the local graph of a facet f at point p.
// Precondition is that p ist part of f.

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_from_facet(Halffacet_handle f, const Point_3& p) const
{ 
  CGAL_nef3_assertion(FM_decorator(f).contains(p));
  Vertex_handle v = sncp()->new_vertex();
  point(v) = q;
  Sphere_circle c(plane(f)); // circle through origin parallel to h
  SM_decorator D(v);
  SHalfloop_handle l = D.new_loop_pair();
  SFace_handle f1 = D.new_face(), f2 = D.new_face();
  D.link_as_loop(l,f1);
  D.link_as_loop(twin(l),f2);

  D.circle(l) = c; D.circle(twin(l)) = c.opposite();
  D.mark(f1) = mark(volume(f));
  D.mark(l) = mark(f);
  Sphere_point q(0,-1,0);
  CGAL::Oriented_side os = c.oriented_side(q);
  switch ( os ) {
    case ON_POSITIVE_SIDE: 
      D.mark_of_halfsphere(-1) = D.mark_of_halfsphere(+1) = true;
      break;
    case ON_NEGATIVE_SIDE:
      D.mark_of_halfsphere(-1) = D.mark_of_halfsphere(+1) = false;
      break;
    case ON_ORIENTED_BOUNDARY:
      if ( c.a()<=0 && c.c()>=0 ) // normal(c) dx<=0&&dz>=0
        D.mark_of_halfsphere(+1) = true;
      if ( c.a()>=0 && c.c()<=0 ) // normal(c) dx<=0&&dz>=0
        D.mark_of_halfsphere(-1) = true;
  }
  return v;
}


// ----------------------------------------------------------------------------
// create_from_edge()
// Creates the local graph of an edge e at point p.
// Precondition is that p ist part of segment(e).

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_from_edge(Halfedge_handle e,
		 const Point_3& p) const
{ CGAL_nef3_assertion(segment(e).has_on(p));
  Vertex_handle v = sncp()->new_vertex();
  point(v) = p;
  Sphere_circle c(h); // circle through origin parallel to h
  SM_decorator D(v);
  SM_const_decorator E(source(e));
  Sphere_point ps = calc_point(e);
  SVertex_handle v1 = D.new_vertex(ps);
  SVertex_handle v2 = D.new_vertex(ps.antipode());
  D.mark(v1) = D.mark(v2) = mark(e);
  bool first = true;
  SHalfedge_around_svertex_const_circulator ec1(E.out_edges(e)), ee(ec1);
  SHalfedge_handle e1,e2;
  CGAL_For_all(ec1,ee) {
    if (first) e1 = D.new_edge_pair(v1,v2);
    else       e1 = D.new_edge_pair(e1, e2, SM_decorator::AFTER, 
				            SM_decorator::BEFORE);
    e2 = D.twin(e1); 
    first = false;
  }
  ec1 = E.out_edges(e);
  SHalfedge_around_svertex_circulator ec2(D.out_edges(v1));
  CGAL_For_all(ec1,ee) {
    D.mark(ec2) = E.mark(ec1);
    D.circle(ec2) = E.circle(ec1);
    D.circle(D.twin(ec2)) = E.circle(E.twin(ec1));
    SFace_handle f = D.new_face();
    D.link_as_new_face_cycle(ec2,f);
    D.mark(f) = E.mark(E.face(ec1));
    ++ec2;
  }
  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
  return v;
}

// ----------------------------------------------------------------------------
// pair_up_halfedges()
// Starting from all local graphs of all vertices of a nef polyhedron
// we pair up all halfedges to halfedge pairs. 

template <typename SNC_>
void SNC_constructor<SNC_>::
pair_up_halfedges() const
{ TRACEN(">>>>>pair_up_halfedges");
  typedef Halfedge_key< Point_3, Halfedge_handle, SNC_decorator>
    Halfedge_key;
  typedef Halfedge_key_lt< Point_3, Halfedge_handle, SNC_decorator> 
    Halfedge_key_lt;
  typedef std::list<Halfedge_key>  Halfedge_list;

  typedef CGAL::Pluecker_line_3<Kernel> Pluecker_line_3;
  typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
  typedef std::map< Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
    Pluecker_line_map;

  SNC_decorator D(*this);
  Pluecker_line_map M;
  Halfedge_iterator e;
  CGAL_nef3_forall_halfedges(e,*sncp()) {
    Point_3 p = point(vertex(e));
    Pluecker_line_3 l(p, p + tmp_point(e));
    TRACEN("  "<<p<<" "<<tmp_point(e)<<" "<<&*e<<" "<<l);
    int inverted;
    l = categorize(l,inverted);
    M[l].push_back(Halfedge_key(p,inverted,e,D));
  }

  typename Pluecker_line_map::iterator it;
  CGAL_nef3_forall_iterators(it,M) {
    it->second.sort(Halfedge_key_lt());
    TRACEN("  "<<it->first<<"\n   "
	   <<(debug_container(it->second),""));
    typename Halfedge_list::iterator itl;
    CGAL_nef3_forall_iterators(itl,it->second) {
      Halfedge_handle e1 = itl->e;
      ++itl; CGAL_nef3_assertion(itl != it->second.end());
      Halfedge_handle e2 = itl->e;
      // TRACEN("e1="<<tmp_point(e1)<<"@"<<point(vertex(e1))<<
      //     " & e2="<<tmp_point(e2)<<"@"<<point(vertex(e2)));
      CGAL_nef3_assertion(tmp_point(e1)==tmp_point(e2).antipode());
      make_twins(e1,e2);
      CGAL_nef3_assertion(mark(e1)==mark(e2));

      // discard temporary sphere_point ?
    }
  }
}

// ----------------------------------------------------------------------------
// link_shalfedges_to_facet_cycles()
// links all edge-uses to facets cycles within the corresponding planes

template <typename SNC_>
void SNC_constructor<SNC_>::
link_shalfedges_to_facet_cycles() const
{
  TRACEN(">>>>>link_shalfedges_to_facet_cycles");
  Halfedge_iterator e;
  CGAL_nef3_forall_edges(e,*sncp()) {
    Halfedge_iterator et = twin(e);
    SM_decorator D(vertex(e)), Dt(vertex(et));
    if ( D.is_isolated(e) ) continue;
    SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
    SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);
    CGAL_For_all(cet,cete) 
      if ( Dt.circle(cet) == D.circle(ce).opposite() ) break;

    /* DEBUG 
    if( Dt.circle(cet) != D.circle(ce).opposite() )
      TRACEN("assertion failed!");

      SHalfedge_around_svertex_circulator sc(D.first_out_edge(e));
      SHalfedge_around_svertex_circulator sct(Dt.first_out_edge(et));
      CGAL_For_all(sc,cee)
	TRACEN("sseg@E addr="<<&*sc<<
	       " src="<<D.point(D.source(sc))<<
	       " tgt="<<D.point(D.target(sc))<<endl<<
	       " circle=" << D.circle(sc));
      CGAL_For_all(sct,cete)
	TRACEN("sseg@ET addr="<<&*sct<<
	       " src="<<Dt.point(Dt.source(sct))<<
	       " tgt="<<Dt.point(Dt.target(sct))<<endl<<
	       " circle=" << Dt.circle(sct));
#ifdef SM_VISUALIZOR
      typedef CGAL::SNC_SM_visualizor<SNC_structure> SMV;
      CGAL::OGL::add_sphere();
      SMV V(vertex(e), CGAL::OGL::spheres_.back());
      V.draw_map();
      SMV Vt(vertex(et), CGAL::OGL::spheres_.back());
      Vt.draw_map();
      CGAL::OGL::start_viewer();
      char c;
      cin >> c;
#endif
      DEBUG */
    CGAL_nef3_assertion( Dt.circle(cet) == D.circle(ce).opposite() ); 
    CGAL_For_all(ce,cee) { 
      CGAL_nef3_assertion(ce->tmp_mark()==cet->tmp_mark());
      link_as_prev_next_pair(Dt.twin(cet),ce);
      link_as_prev_next_pair(D.twin(ce),cet);
      --cet; // ce moves ccw, cet moves cw
    }
  }
}

// ----------------------------------------------------------------------------
// categorize_facet_cycles_and_create_facets()
// sweeping all edge-uses we categorize facet cycle incidence, create
// the facet objects and assign the facet cycles.

template <typename SNC_>
void SNC_constructor<SNC_>::
categorize_facet_cycles_and_create_facets() const
{ TRACEN(">>>>>categorize_facet_cycles_and_create_facets");

  typedef std::list<SObject_handle> SObject_list;
  typedef std::map<Plane_3, SObject_list, Plane_lt> 
    Map_planes;

  Map_planes M;
  SHalfedge_iterator e;
  CGAL_nef3_forall_shalfedges(e,*sncp()) {
    Sphere_circle c(tmp_circle(e));
    Plane_3 h = c.plane_through(point(vertex(e))); 
    if ( sign_of(h)<0 ) continue;
    M[normalized(h)].push_back(SObject_handle(twin(e)));
  }
  SHalfloop_iterator l;
  CGAL_nef3_forall_shalfloops(l,*sncp()) {
    Sphere_circle c(tmp_circle(l));
    Plane_3 h = c.plane_through(point(vertex(l))); 
    if ( sign_of(h)<0 ) continue;
    M[normalized(h)].push_back(SObject_handle(twin(l)));
  }
  typename Map_planes::iterator it;
  CGAL_nef3_forall_iterators(it,M) { TRACEN("  plane "<<it->first);
    FM_decorator D(*sncp());
    D.create_facet_objects(it->first,it->second.begin(),it->second.end());
  }
}

// ----------------------------------------------------------------------------
// create_volumes()
// categorizes all shells and creates volume objects.

template <typename SNC_>
void SNC_constructor<SNC_>::
create_volumes() const
{ 
  TRACEN(">>>>>create_volumes");
  Shell_number_hash  Shell(-1);
  Shell_closed_hash Closed(true);
  SFace_visited_hash Done(false);
  Shell_explorer V(*this,Shell,Closed,Done);
  std::vector<Vertex_handle> MinimalVertex;
  std::vector<SFace_handle> EntrySFace;
  SFace_iterator f;
  /* First, we classify all the Shere Faces per Shell.  For each Shell we
     determine its minimum lexicographyly vertex and we check wheter the
     Shell encloses a region (closed surface) or not. */
  CGAL_nef3_forall_sfaces(f,*sncp()) { 
    if ( Done[f] ) 
      continue;
    V.minimal_vertex() = vertex(f);
    visit_shell_objects(f,V);

    MinimalVertex.push_back(V.minimal_vertex());
    EntrySFace.push_back(f);
    V.increment_shell_number();
  }
  /* then, we determine the Shells which correspond to Volumes via a ray
     shotting in the direction (-1,0,0) over the Sphere_map of the minimal 
     vertex.  The Shell corresponds to a Volume if the object hit belongs 
     to another Shell. */
  sncp()->new_volume(); // outermost volume (nirvana)
  for( unsigned int i = 0; i < MinimalVertex.size(); ++i) {
    Vertex_handle v = MinimalVertex[i];
    TRACEN( "Shell #" << i << " minimal vertex: " << point(v));
    SM_point_locator D(v);
    SObject_handle o = D.locate(Sphere_point(-1,0,0));
    SFace_const_handle sfc;
    if( !assign(sfc, o) || Shell[sfc] != i) { /*UNTESTED CASE: !assign(sfc,o)*/
      SFace_handle f = EntrySFace[i];
      CGAL_assertion( Shell[EntrySFace[i]] == i );
      if( Closed[f] ) {
	SM_decorator SD(v);
	Volume_handle c = sncp()->new_volume();
	mark(c) = SD.mark(f);
	link_as_outer_shell(f, c );
	TRACE( "Shell #" << i <<" linked as outer shell");
	TRACEN( "(sface" << (assign(sfc,o)?"":" not") << " hit case)");
      }
    }
  }
  /* finaly, we go through all the Shells which do not correspond to a Volume 
     and we assign them to its enclosing Volume determined via a facet below
     check. */
  CGAL_nef3_forall_sfaces(f,*sncp()) {
    if ( volume(f) != Volume_handle() ) 
      continue;
    TRACEN( "Inner shell #" << Shell[f] << " volume?");
    Volume_handle c = determine_volume( f, MinimalVertex, Shell );
    link_as_inner_shell( f, c );
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONSTRUCTOR_H

