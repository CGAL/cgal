#ifndef CGAL_SNC_CONSTRUCTOR_H
#define CGAL_SNC_CONSTRUCTOR_H

#include <CGAL/basic.h>
#include <CGAL/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <map>
#include <list>
#undef _DEBUG
#define _DEBUG 43
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename Point, typename Edge>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge> Self;
  Point p; int i; Edge e;
  Halfedge_key(Point pi, int ii, Edge ei) : p(pi), i(ii), e(ei) {}
  Halfedge_key(const Self& k) { p=k.p; i=k.i; e=k.e; }
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge> Key;
  bool operator()(const Key& k1, const Key& k2) const
  { if ( k1.p == k2.p ) return (k1.i < k2.i);
    return CGAL::lexicographically_xyz_smaller(k1.p,k2.p); }
};

template <typename Point, typename Edge>
std::ostream& operator<<(std::ostream& os, 
                         const Halfedge_key<Point,Edge>& k )
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
  typedef typename SNC_structure_::Sphere_kernel Sphere_kernel;
  typedef typename SNC_structure_::Kernel        Kernel;
  typedef SNC_constructor<SNC_structure_>        Self;
  typedef SNC_decorator<SNC_structure_>        Base;
  typedef SNC_decorator<SNC_structure_>        SNC_decorator;
  typedef SNC_FM_decorator<SNC_structure_>     FM_decorator;
  typedef SNC_SM_decorator<SNC_structure_>     SM_decorator;
  typedef SNC_SM_overlayer<SNC_structure_>     SM_overlayer;
  typedef SNC_SM_point_locator<SNC_structure_> SM_point_locator;
  typedef SNC_SM_const_decorator<SNC_structure_> SM_const_decorator;

  #define USING(t) typedef typename SNC_structure_::t t
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Volume_iterator);

  USING(Vertex_const_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_handle);

  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SFace_iterator);
  USING(SHalfloop_iterator);

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

  USING(SFace_cycle_iterator);
  USING(SFace_cycle_const_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Halffacet_cycle_const_iterator);
  USING(Shell_entry_iterator);
  USING(Shell_entry_const_iterator);

  USING(Point_3);
  USING(Vector_3);
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

  typedef CGAL::Unique_hash_map<SFace_handle,int>  Shell_number_hash;
  typedef CGAL::Unique_hash_map<SFace_handle,bool> SFace_visited_hash;
  typedef CGAL::Unique_hash_map<SFace_handle,bool> Shell_closed_hash;

  struct Shell_explorer {
    const SNC_decorator& D;
    Shell_number_hash&  Shell;
    SFace_visited_hash& Done;
    Shell_closed_hash& Closed;
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

    Vertex_handle& minimal_vertex() { 
      return v_min; 
    }
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

  Volume_handle determine_volume( SFace_handle& f, 
    const std::vector<Vertex_handle>& MinimalVertex, 
    const Shell_number_hash&  Shell ) const {
    Vertex_handle v_min = MinimalVertex[Shell[f]];
    Halffacet_handle f_below = get_facet_below(v_min);
    if( f_below == Halffacet_handle() )
      // return volumes_begin(); // is volumes_begin() a const iterator?
      return SNC_decorator(*this).volumes_begin();
    Volume_handle c = volume(f_below);
    if(c != Volume_handle())
      return c;
    SFace_handle sf = sface_custom(f_below);
    c = determine_volume(sf, MinimalVertex, Shell);
    link_as_inner_shell(sf, c);
    return c;
  }

  Halffacet_handle get_facet_below(const Vertex_handle& vi) const {
    //Segment_3 s(point(v), point(vertices_begin())), // vertices_begin()?
    Segment_3 s(point(vi), point(Base(*this).vertices_begin())),
      so(s.opposite());
    Object_handle closest;
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    CGAL_forall_vertices(v,*sncp()) {
      SFace_iterator sf = v->sfaces_begin();
      if (sf != 0 && volume(sf) != Volume_handle() 
	  && contains_internally(s,point(v))) {
	shorten(s,point(v));
	closest = Object_handle(v);
      }
    }
    CGAL_forall_edges(e,*sncp()) {
      Point_3 q;
      SFace_handle sf = sface(e);
      if (sf != SFace_handle() && volume(sf) != Volume_handle()
	  && do_intersect(s,e,q) ) { 
	shorten(s,q); 
	closest = Object_handle(e);
      }
    }
    CGAL_forall_halffacets(f,*sncp()) {
      Point_3 q;
      SFace_handle sf = sface_custom(f);
      if (sf != SFace_handle() && volume(sf) != Volume_handle()
	  && do_intersect(s,f,q) ) { 
	shorten(s,q); 
	closest = Object_handle(f); 
      }
    }
    if( assign(closest, v) ) {
      SHalfedge_handle se = v->shalfedges_begin();
      if( se != SHalfedge_handle() )
	return facet(se);
      SHalfloop_handle sl = v->shalfloop();
      if( sl != SHalfloop_handle() )
	return facet(sl);
      CGAL_assertion_msg(0, "Empty local map.");
    } 
    else if( assign( closest, e ) ) {
      SM_decorator SD(vertex(e));
      CGAL_assertion( SD.first_out_edge(e) != SHalfedge_handle() );
      return facet(SD.first_out_edge(e));
    }
    else if( assign( closest, f ) )
      return f;
    else
      TRACEN("no facet below found");
    return Halffacet_handle();
  }

  /* following 4 functions were copied from SNC_point_locator.h */
  void shorten(Segment_3& s, const Point_3& p) const
  { s = Segment_3(s.source(),p); }

  bool contains_internally(const Segment_3& s, const Point_3& p) const
  { if(!s.has_on(p))
      return false;
    Comparison_result r1 = compare_xyz(s.source(),p); 
    Comparison_result r2 = compare_xyz(s.target(),p); 
    return (r1 == opposite(r2));
  }

  bool do_intersect(const Segment_3& s,
		    Halfedge_handle e,
		    Point_3& p) const {
    bool intersect; //= do_intersect(s,segment(e),p);
    return intersect;
  }

#ifdef LINE3_LINE3_INTERSECTION
  bool do_intersect( const Segment_3& s, 
		     const Segment_3& q, 
		     Point_3& p) const 
  {
    Object o = intersection(Line_3(s),Line_3(q)); 
    if ( !assign(p,o) ) 
      return false;
    if( !contains_internally(s, p) || !contains_internally(q, p) )
      return false;
    return true;
  }

#else // LINE3_LINE3_INTERSECTION

  bool do_intersect( const Segment_3& s, 
		     const Segment_3& r, 
		     Point_3& p) const 
  {
    if(s.is_degenerate() || r.is_degenerate())
      return false;
    /* at least one of the segments is degenerate so 
       there is not internal intersection */
    if(orientation(s.source(),s.target(),r.source(),r.target()) != COPLANAR)
      return false;
    /* the segments doesn't define a plane */
    if(collinear(s.source(),s.target(),r.source()) &&
       collinear(s.source(),s.target(),r.target()) )
      return false;
    /* the segments are collinear */
    Line_3 ls(s), lr(r);
    if(ls.direction() ==  lr.direction() ||
       ls.direction() == -lr.direction() )
      return false;
    /* the segments are parallel */

    Oriented_side os1, os2;
    Vector_3 vs(s.direction()), vr(r.direction()), vt(cross_product(vs, vr)), 
      ws(cross_product(vt, vs)), wr(cross_product(vt, vr));
    Plane_3 hs(s.source(),ws);
    /* hs is a plane which contains line(s) and is perpendicular to the
       plane defined by s and r */
    os1 = hs.oriented_side(r.source());
    os2 = hs.oriented_side(r.target());
    if(os1 != opposite(os2))
      return false;
    Plane_3 hr(r.source(),wr);
    /* hr is a plane which contains line(r) and is perpendicular to the
       plane defined by s and r */
    os1 = hr.oriented_side(s.source());
    os2 = hr.oriented_side(s.target());
    if(os1 != opposite(os2))
      return false;

    Object o = intersection(hs, lr);
    CGAL_assertion(assign(p,o));
    /* since line(s) and line(r) are not parallel they intersects in only
       one point */
    assign(p,o);
    if( !contains_internally(s, p) || !contains_internally(r, p) )
      return false;
    return true;
  }
#endif // LINE3_LINE3_INTERSECTION

  bool do_intersect(const Segment_3& s,
		    Halffacet_handle f,
		    Point_3& p) const
  { Plane_3 h(plane(f));
    Object o = intersection(h, Line_3(s));
    if ( !CGAL::assign(p,o) ) 
      return false;
    Oriented_side os1 = h.oriented_side(s.source());
    Oriented_side os2 = h.oriented_side(s.target());
    if (os1 != opposite(os2))
      return false;
    return (locate_point_in_halffacet(p, f) == CGAL::ON_BOUNDED_SIDE);
  }
  
  Bounded_side locate_point_in_halffacet( const Point_3& p, 
					  const Halffacet_handle& f) const {
    Plane_3 h(plane(f));
    CGAL_assertion(h.has_on(p));
    Halffacet_cycle_iterator fc = f->facet_cycles_begin();
    SHalfedge_handle e;
    Bounded_side outer_bound_pos;
    if ( assign(e,fc) ) {
      vector<Point_3> verts;
      SHalfedge_handle e0(e);
      CGAL_For_all(e,e0)
	verts.push_back(point(vertex(e)));
      outer_bound_pos = bounded_side_3(verts.begin(), verts.end(), p, h);
      verts.clear();
    } else CGAL_assertion_msg(0, "facet's first cycle is a SHalfloop?");
    if( outer_bound_pos != CGAL::ON_BOUNDED_SIDE )
    /* if point p is not in the relative interior of the outer face cycle
       is not necesary to know the possition of p with respect of the 
       inner face cycles */
      return outer_bound_pos;
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
      } else if ( assign(e,fc) ) {
	vector<Point_3> verts;
	SHalfedge_handle e0(e);
	CGAL_For_all(e,e0)
	  verts.push_back(point(vertex(e)));
	inner_bound_pos = bounded_side_3(verts.begin(), verts.end(), p, 
					 h.opposite());
	verts.clear();
      } else CGAL_assertion_msg(0, "Damn wrong handle.");
      if( inner_bound_pos != CGAL::ON_UNBOUNDED_SIDE )
	/* the relative position of the point p in the relative interior
	   of the outer cycle of a face f is known when p belongs to the
	   clousure of any inner face cycle */
	return opposite(inner_bound_pos);
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
create_box_corner(int x, int y, int z,
                  bool boundary=true) const
{ 
  CGAL_nef3_assertion(x*y*z != 0);
  Vertex_handle v = sncp()->new_vertex();
  int R=3; point(v) = Point_3(x*R,y*R,z*R);
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
  CGAL_forall_sedges_of(e,v)
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
  typedef Halfedge_key<Point_3,Halfedge_handle>    Halfedge_key;
  typedef Halfedge_key_lt<Point_3,Halfedge_handle> Halfedge_key_lt;
  typedef std::list<Halfedge_key>                  Halfedge_list;

  typedef CGAL::Pluecker_line_3<Kernel> Pluecker_line_3;
  typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
  typedef std::map<Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
    Pluecker_line_map;

  Pluecker_line_map M;
  Halfedge_iterator e;
  CGAL_forall_halfedges(e,*sncp()) {
    Point_3 p = point(vertex(e));
    Pluecker_line_3 l(p, p + tmp_point(e));
    TRACEN("  "<<p<<" "<<tmp_point(e)<<" "<<&*e<<" "<<l);
    int inverted;
    l = categorize(l,inverted);
    M[l].push_back(Halfedge_key(p,inverted,e));
  }

  Pluecker_line_map::iterator it;
  CGAL_forall_iterators(it,M) {
    it->second.sort(Halfedge_key_lt());
    TRACEN("  "<<it->first<<"\n   "
	   <<(debug_container(it->second),""));
    Halfedge_list::iterator itl;
    CGAL_forall_iterators(itl,it->second) {
      Halfedge_handle e1 = itl->e;
      ++itl; CGAL_nef3_assertion(itl != it->second.end());
      Halfedge_handle e2 = itl->e;
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
  CGAL_forall_edges(e,*sncp()) {
    Halfedge_iterator et = twin(e);
    SM_decorator D(vertex(e)), Dt(vertex(et));
    if ( D.is_isolated(e) ) continue;
    SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
    SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);
    /*debug*/
    SHalfedge_around_svertex_circulator ceti, cetj;
    cetj = cet; ceti = cetj++;
    if(ceti != cetj)
      CGAL_For_all(ceti,cete) 
      {
	TRACEN(" et   ss "<<Dt.point(Dt.source(ceti))<<" "
	       <<Dt.point(Dt.target(ceti))<<" "<<&(*ceti));
	TRACEN(" et+1 ss "<<Dt.point(Dt.source(cetj))<<" "
	       <<Dt.point(Dt.target(cetj))<<" "<<&(*cetj));
	//CGAL_nef3_assertion( Dt.circle(ceti) != D.circle(cetj));
	cetj++;
      }
    /*debug*/
    CGAL_For_all(cet,cete) 
      {
	TRACEN(" e ss circle "<<D.circle(ce));
	TRACEN(" e vertices  "<<D.point(D.source(ce))<<" "
	       <<D.point(D.target(ce)));
	TRACEN(" et ss circle "<<Dt.circle(cet));
	TRACEN(" et vertices  "<<Dt.point(Dt.source(cet))<<" "
	       <<Dt.point(Dt.target(cet)));
	TRACEN("are circles opossite? "
	       <<(D.circle(ce).opposite() == Dt.circle(cet)));
	if ( Dt.circle(cet) == D.circle(ce).opposite() ) break;
      }
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
  CGAL_forall_shalfedges(e,*sncp()) {
    Sphere_circle c(tmp_circle(e));
    Plane_3 h = c.plane_through(point(vertex(e))); 
    if ( sign_of(h)<0 ) continue;
    M[normalized(h)].push_back(SObject_handle(twin(e)));
  }
  SHalfloop_iterator l;
  CGAL_forall_shalfloops(l,*sncp()) {
    Sphere_circle c(tmp_circle(l));
    Plane_3 h = c.plane_through(point(vertex(l))); 
    if ( sign_of(h)<0 ) continue;
    M[normalized(h)].push_back(SObject_handle(twin(l)));
  }
  Map_planes::iterator it;
  CGAL_forall_iterators(it,M) { TRACEN("  plane "<<it->first);
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
  CGAL_forall_sfaces(f,*sncp()) { 
    if ( Done[f] ) 
      continue;
    V.minimal_vertex() = vertex(f);
    visit_shell_objects(f,V);
    MinimalVertex.push_back(V.minimal_vertex());
    EntrySFace.push_back(f);
    V.increment_shell_number();
  }
  Volume_handle outer_volume = sncp()->new_volume();
  for (unsigned i = 0; i < MinimalVertex.size(); ++i) {
    TRACEN("minimal vertex "<<i<<" "<<point(MinimalVertex[i]));
    Vertex_handle v = MinimalVertex[i];
    if(true) { // v is not a bounding box vertex, TODO
      SM_point_locator D(v);
      SObject_handle o = D.locate_no_const(Sphere_point(-1,0,0));
      SFace_handle sf;
      if ( !assign(sf,o) || Shell[sf] != i ) {
	TRACEN("outer shell "<<i);
	SVertex_handle   sv;
	SHalfedge_handle se;
	SHalfloop_handle sl;
	if( !assign(sf,o) ) {
	  if( assign(sv,o) )
	    f = sface(sv);
	  else if( assign(se,o) )
	    f = sface(se);
	  else if( assign(sl,o) )
	    f = sface(sl);
	  else CGAL_assertion_msg(0,"Damn wrong handle");
	}
	TRACEN("closed shell? "<<Closed[sf]);
	if( true ) { // Closed[f] ) {
	  Volume_handle c = sncp()->new_volume();
	  link_as_outer_shell(sf, c );
	}
      }
    }
  }
  
  CGAL_forall_sfaces(f,*sncp()) {
    if ( volume(f) != Volume_handle() ) 
      continue;
    Volume_handle c = determine_volume( f, MinimalVertex, Shell );
    link_as_inner_shell( f, c );
  }
  return;
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONSTRUCTOR_H

