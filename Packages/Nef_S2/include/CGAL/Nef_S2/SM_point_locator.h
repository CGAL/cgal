// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_S2/SM_point_locator.h
// package       : Nef_S2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Point location module
// ============================================================================
#ifndef CGAL_SM_POINT_LOCATOR_H
#define CGAL_SM_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Object.h>
#include <CGAL/Nef_2/geninfo.h>

#undef _DEBUG
#define _DEBUG 143
#include <CGAL/Nef_S2/debug.h>


CGAL_BEGIN_NAMESPACE

/*{\Moptions print_title=yes }*/ 

/*{\Manpage {SM_point_locator}{Decorator}
{Naive point location in plane maps}{PL}}*/
/*{\Mdefinition An instance |\Mvar| of data type |\Mname|
encapsulates naive point location queries within a sphere map |M|.  The
two template parameters are specified via concepts. |PM_decorator_|
must be a model of the concept |PMDecorator| as described in the
appendix.  |Geometry_| must be a model of the concept
|AffineGeometryTraits_2| as described in the appendix. For a
specification of plane maps see also the concept of
|PMConstDecorator|.}*/

/*{\Mgeneralization Decorator}*/

template <typename Decorator_>
class SM_point_locator : public Decorator_ {
protected:
  typedef Decorator_ Base;
  typedef SM_point_locator<Decorator_> Self;

public:
  /*{\Mtypes 5}*/
  typedef Decorator_ Decorator;
  /*{\Mtypemember equals |Decorator_|.}*/
  typedef typename Decorator::Sphere_map Sphere_map;
  /*{\Mtypemember the sphere map type decorated by |Decorator|.}*/
  typedef typename Decorator::Mark      Mark;
  /*{\Mtypemember the attribute of all objects (vertices, edges, loops,
  faces).}*/

  typedef typename Decorator::Kernel Kernel;
  /*{\Mtypemember the sphere kernel.}*/
  typedef typename Kernel::Sphere_point Sphere_point;
  /*{\Mtypemember points.}*/
  typedef typename Kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember segments.}*/
  typedef typename Kernel::Sphere_circle Sphere_circle;
  /*{\Mtypemember circles.}*/
  typedef typename Kernel::Sphere_direction Sphere_direction;
  /*{\Mtypemember directions.}*/

  /*{\Mtext Local types are handles, iterators and circulators of the 
  following kind: |Vertex_const_handle|, |Vertex_const_iterator|, 
  |Halfedge_const_handle|, |Halfedge_const_iterator|, 
  |Halfloop_const_handle|, |Halfloop_const_iterator|, 
  |Face_const_handle|, |Face_const_iterator|.}*/

  typedef typename Sphere_map::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying plane
  map. The kind of the object |(vertex, halfedge,face)| can be determined and
  the object assigned by the three functions:\\ 
  |bool assign(Vertex_const_handle& h, Object_handle o)|\\ 
  |bool assign(Halfedge_const_handle& h, Object_handle o)|\\ 
  |bool assign(Face_const_handle& h, Object_handle o)|\\ where each 
  function returns |true| iff the assignment of |o| to |h| was valid.}*/

  #define USING(t) typedef typename Decorator_::t t
  USING(Vertex_handle);   
  USING(Halfedge_handle);   
  USING(Halfloop_handle);   
  USING(Face_handle);     
  USING(Vertex_const_handle); 
  USING(Halfedge_const_handle); 
  USING(Halfloop_const_handle); 
  USING(Face_const_handle); 
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halfloop_iterator);
  USING(Face_iterator);
  USING(Vertex_const_iterator);
  USING(Halfedge_const_iterator);
  USING(Halfloop_const_iterator);
  USING(Face_const_iterator);
  USING(Halfedge_around_vertex_circulator);
  USING(Halfedge_around_vertex_const_circulator);
  USING(Halfedge_around_face_circulator);
  USING(Halfedge_around_face_const_circulator);
  #undef USING

  Sphere_segment segment(Halfedge_const_handle e) const
  { return Sphere_segment(point(source(e)), point(target(e)), circle(e)); }

  Sphere_direction direction(Halfedge_const_handle e) const
  { return Sphere_direction(circle(e)); }

  Halfedge_const_handle out_wedge(Vertex_const_handle v, 
    const Sphere_direction& d, bool& collinear) const
  /*{\Xop returns a halfedge |e| bounding a wedge in between two
  neighbored edges in the adjacency list of |v| which contains |d|.
  If |d| extends along a edge then |e| is this edge. If |d| extends
  into the interior of such a wedge then |e| is the first edge hit
  when |d| is rotated clockwise. \precond |v| is not isolated.}*/
  { TRACEN("out_wedge "<<PH(v));
    assert(!is_isolated(v));
    collinear=false;
    Sphere_point p = point(v);
    Halfedge_const_handle e_res = first_out_edge(v);
    Sphere_direction d_res = direction(e_res);
    Halfedge_around_vertex_const_circulator el(e_res),ee(el);
    CGAL_For_all(el,ee) {
      if ( strictly_ordered_ccw_at(p,d_res, direction(el), d) )
        e_res = el; d_res = direction(e_res);
    }
    TRACEN("  determined "<<PH(e_res)<<" "<<d_res);
    if ( direction(cyclic_adj_succ(e_res)) == d ) {
      e_res = cyclic_adj_succ(e_res);
      collinear=true;
    }
    TRACEN("  wedge = "<<PH(e_res)<<" "<<collinear);
    return e_res;
  }


  template <class Handle>
  Object_handle make_object(Handle h) const
  { return CGAL::make_object(h); }

  /*{\Mcreation 3}*/

  SM_point_locator() : Base() {}

  /*{\Moptions constref=yes}*/
  SM_point_locator(const Sphere_map& M) : 
    Base(const_cast<Sphere_map&>(M)) {}
  /*{\Mcreate constructs a point locator working on |P|.}*/
  /*{\Moptions constref=no}*/
  /*{\Moperations 2.5 0.5}*/

  const Mark& mark(Object_handle h) const
  /*{\Mop returns the mark associated to the object |h|.}*/
  { Vertex_const_handle v; 
    Halfedge_const_handle e; 
    Face_const_handle f;
    if ( assign(v,h) ) return mark(v);
    if ( assign(e,h) ) return mark(e);
    if ( assign(f,h) ) return mark(f);
    CGAL_assertion_msg(0,
    "PM_point_locator::mark: Object_handle holds no object.");
  }


  enum SOLUTION { is_vertex_, is_edge_, is_loop_ };
  // enumeration for internal use
  
  Object_handle locate(const Sphere_point& p) const
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of the underlying plane map |P| which contains the point |p =
  s.source()| in its relative interior. |s.target()| must be a point
  such that |s| intersects the $1$-skeleton of |P|.}*/
  { TRACEN("locate naivly "<<p);
    Vertex_const_iterator v;
    CGAL_forall_vertices(v,*this) {
      if ( p == point(v) ) return make_object(v);
    }

    Halfedge_const_iterator e;
    CGAL::Unique_hash_map<Halfedge_const_handle,bool> visited(false);
    CGAL_forall_halfedges(e,*this) {
      if ( visited[e] ) continue;
      if ( segment(e).has_on(p) ) return make_object(e);
      visited[e]=visited[twin(e)]=true;
    }
    if ( has_loop() && circle(halfloop()).has_on(p) )
      return make_object(Halfloop_const_handle(halfloop()));

    // now in face:

    if ( number_of_vertices() == 0 && ! has_loop() ) 
      return make_object(faces_begin());

    Vertex_const_handle v_res;
    Halfedge_const_handle e_res;
    Halfloop_const_handle l_res;
    SOLUTION solution;

    Sphere_segment s; // we shorten the segment iteratively
    if ( has_loop() ) {
      Sphere_circle c(circle(halfloop()),p); // orthogonal through p
      s = Sphere_segment(p,intersection(c,circle(halfloop())));
      l_res = circle(halfloop()).has_on_positive_side(p) ? 
	halfloop() : twin(halfloop());
      solution = is_loop_;
    } else { // has vertices !
      CGAL_assertion( number_of_vertices()!=0 );
      Vertex_const_handle vt = vertices_begin();
      Sphere_point pvt = point(vt);
      if ( p != pvt.opposite() ) s = Sphere_segment(p,pvt);
      else s = Sphere_segment(p,pvt,Sphere_circle(p,pvt));
      v_res = vt;
      solution = is_vertex_;
    }

    // s now initialized
    
    Sphere_direction dso(s.sphere_circle().opposite()), d_res;
    visited.clear(false);
    CGAL_forall_vertices(v,*this) {
      Sphere_point p_res, vp = point(v);
      if ( s.has_on(vp) ) {
        TRACEN(" location via vertex at "<<vp);
        s = Sphere_segment(p,vp,s.sphere_circle()); // we shrink the segment
        if ( is_isolated(v) ) {
          v_res = v; solution = is_vertex_;
        } else { // not isolated
          bool dummy;
          e_res = out_wedge(v,dso,dummy);
          Halfedge_around_vertex_const_circulator el(e_res),ee(el);
          CGAL_For_all(el,ee) 
            visited[el] = visited[twin(el)] = true;
          /* e_res is now the counterclockwise maximal halfedge out
             of v just before s */
          if ( circle(e_res).has_on_negative_side(p) )
            e_res = previous(e_res);
            // correction to make e_res visible from p
	  solution = is_edge_;
          TRACEN("  determined "<<PH(e_res));
        }
      }
    }

    CGAL_forall_halfedges(e,*this) {
      if ( visited[e] ) continue;
      Sphere_segment se = segment(e);
      Sphere_point p_res;
      if ( do_intersect_internally(se,s,p_res) ) {
          TRACEN(" location via halfedge "<<se);
        s = Sphere_segment(p,p_res,s.sphere_circle()); 
        e_res = ( circle(e).has_on_positive_side(p) ? e : twin(e) );
        visited[e] = visited[twin(e)] = true;
	solution = is_edge_;
        TRACEN("  determined "<<PH(e_res)<<" "<<mark(face(e_res)));
      }
    }

    switch ( solution ) {
      case is_edge_: 
        return make_object((Face_const_handle)(face(e_res)));
      case is_loop_:
        return make_object((Face_const_handle)(face(l_res)));
      case is_vertex_:
        return make_object((Face_const_handle)(face(v_res)));
      default: CGAL_assertion_msg(0,"missing solution.");
    }
    return Object_handle(); // never reached!
  }

  
  template <typename Object_predicate>
  Object_handle ray_shoot(const Sphere_point& p, 
			  const Sphere_direction& d,
			  const Object_predicate& M) const
  /*{\Mop returns an |Object_handle o| which can be converted to a
  |Vertex_const_handle|, |Halfedge_const_handle|, |Face_const_handle|
  |h| as described above.  The object predicate |M| has to have
  function operators \\ |bool operator() (const
  Vertex_/Halfedge_/Halfloop_/Face_const_handle&)|.\\ The object
  returned is intersected by |d.circle()|, has minimal distance to
  |p|, and |M(h)| holds on the converted object. The operation returns
  the null handle |NULL| if the ray shoot along |s| does not hit any
  object |h| of |M| with |M(h)|.}*/
  { 
    Sphere_circle c(d.circle());
    Sphere_segment s;
    bool s_init(false);
    Object_handle h = locate(p);
    Vertex_const_handle v; 
    Halfedge_const_handle e; 
    Halfloop_const_handle l; 
    Face_const_handle f;
    if ( assign(v,h) && M(v) ||
         assign(e,h) && M(e) ||
	 assign(l,h) && M(l) ||
         assign(f,h) && M(f) ) return h;
    h = Object_handle(); 
    TRACEN("not contained");
#if 0
    HASEN: s am anfang circle, ab wann segment ?
	   wo loop ?

    CGAL_forall_vertices (v,*this) {
      Point pv = point(v);
      if ( !(s_init && s.has_on(pv) ||
	    !s_init && c.has_on(pv)) ) continue;
      TRACEN("candidate "<<pv);
      if ( M(v) ) {
        h = make_object(v);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
        continue;
      }
      // now we know that v is not marked but on s
      bool collinear;
      Halfedge_const_handle e = out_wedge(v,d,collinear);
      if ( collinear ) { 
        if ( M(e) ) {
          h = make_object(e);
          s = Sphere_segment(p,pv,c);
        }
        continue;
      }
      if ( M(face(e)) ) {
        h = make_object(face(e));
        s = Sphere_segment(p,pv,c);
      }
    } // all vertices

    CGAL::Unique_hash_map<Halfedge_const_handle,bool> visited(false);
    Halfedge_const_iterator e_res;
    CGAL_forall_halfedges(e,*this) {
      Sphere_segment se = segment(e);
      Sphere_point p_res;
      if ( do_intersect_internally(se,s,p_res) ) {
        // internal intersection
        TRACEN("candidate "<<se); 
        e_res = e;
        Sphere_segment s_cand = Sphere_segment(p,p_res,c);  
        if ( s_cand.is_short() && circle(e).has_on_negative_side(p) ||
	     s_cand.is_long() && circle(e).has_on_positive_side(p) ||
	     s_cand.is_halfcircle() && 
	       strictly_ordered_ccw_at(p.opposite(),
				       direction(e),d,direction(twin(e))) )
	  e_res = twin(e);
        if ( M(e_res) ) {
          h = make_object(e_res); s = s_cand;
        } else if ( M(face(twin(e_res))) ) {
          h = make_object(face(twin(e_res))); s = s_cand;
        }
      }
    }
#endif
    CGAL_assertion_msg(0,"not yet correct");
    return h;
  }

  void init_marks_of_halfspheres();

  // C++ is really friendly:
  #define USECMARK(t) const Mark& mark(t h) const { return Base::mark(h); }
  #define USEMARK(t)  Mark& mark(t h) const { return Base::mark(h); }
  USEMARK(Vertex_handle)
  USEMARK(Halfedge_handle)
  USEMARK(Face_handle)
  USECMARK(Vertex_const_handle)
  USECMARK(Halfedge_const_handle)
  USECMARK(Face_const_handle)
  #undef USEMARK
  #undef USECMARK
  /*{\Mimplementation Naive query operations are realized by checking
  the intersection points of the $1$-skeleton of the plane map |P| with
  the query segments $s$. This method takes time linear in the size $n$
  of the underlying plane map without any preprocessing.}*/
}; // SM_point_locator<Decorator_>

template <typename D>
void SM_point_locator<D>::init_marks_of_halfspheres()
{ TRACEN("init_marks_of_halfspheres");
  Sphere_point y_minus(0,-1,0);
  Object_handle h = locate(y_minus);
  Face_const_handle f;
  if ( CGAL::assign(f,h) ) {
    mark_of_halfsphere(-1) = mark_of_halfsphere(+1) = mark(f);
    return;
  }

  Halfedge_const_handle e;
  if ( CGAL::assign(e,h) ) { 
    CGAL_assertion(circle(e).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(e).orthogonal_vector());
    TRACEN("on edge "<<op);
    if ( (op.x() > 0) || (op.x() == 0) && (op.z() < 0) ) e = twin(e);
    // if ( (op.z() < 0) || (op.z() == 0) && (op.x() > 0) ) e = twin(e);
    mark_of_halfsphere(+1) = mark(face(e));
    mark_of_halfsphere(-1) = mark(face(twin(e)));
    return;
  }

  Halfloop_const_handle l;
  if ( CGAL::assign(l,h) ) {
    CGAL_assertion(circle(l).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(l).orthogonal_vector());
    TRACEN("on loop "<<op);
    if ( (op.x() > 0) || (op.x() == 0) && (op.z() < 0) ) l = twin(l);
    // if ( (op.z() < 0) || (op.z() == 0) && (op.x() > 0) ) l = twin(l);
    mark_of_halfsphere(+1) = mark(face(l));
    mark_of_halfsphere(-1) = mark(face(twin(l)));
    return;
  }

  //Sphere_circle c(-1,0,0);
  //Sphere_direction up(c),down(c.opposite());
  Sphere_circle c(0,0,1);
  Sphere_direction right(c),left(c.opposite());
  bool collinear(false);
  Vertex_const_handle v;
  if ( CGAL::assign(v,h) ) {
    CGAL_assertion(point(v)==y_minus);
    e = out_wedge(v,left,collinear); 
    if ( collinear ) mark_of_halfsphere(+1) = mark(face(twin(e)));
    else mark_of_halfsphere(+1) = mark(face(e));
    e = out_wedge(v,right,collinear); 
    if ( collinear ) mark_of_halfsphere(-1) = mark(face(twin(e)));
    else mark_of_halfsphere(-1) = mark(face(e));
    return;
  }
  CGAL_assertion_msg(0,"damn wrong type.");
}

CGAL_END_NAMESPACE
#endif // CGAL_SM_POINT_LOCATOR_H



