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
// file          : include/CGAL/Nef_3/SNC_SM_point_locator.h
// package       : Nef_3 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Point location module of sphere map
// ============================================================================
#ifndef CGAL_SNC_SM_POINT_LOCATOR_H
#define CGAL_SNC_SM_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_3/SNC_SM_decorator.h>
#undef _DEBUG
#define _DEBUG 47
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

/*{\Moptions print_title=yes }*/ 

/*{\Manpage {SNC_SM_point_locator}{Refs_}
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

template <typename Refs_>
class SNC_SM_point_locator : public SNC_SM_decorator<Refs_> {
protected:
  typedef SNC_SM_decorator<Refs_>     Base;
  typedef SNC_SM_point_locator<Refs_> Self;

public:
  /*{\Mtypes 5}*/
  typedef SNC_SM_decorator<Refs_> Decorator;
  /*{\Mtypemember equals |Decorator_|.}*/

  typedef typename Refs_::Mark      Mark;
  /*{\Mtypemember the attribute of all objects (vertices, edges, loops,
  faces).}*/

  typedef typename Refs_::Sphere_kernel Sphere_kernel;
  /*{\Mtypemember the sphere kernel.}*/
  typedef typename Refs_::Sphere_point Sphere_point;
  /*{\Mtypemember points.}*/
  typedef typename Refs_::Sphere_segment Sphere_segment;
  /*{\Mtypemember segments.}*/
  typedef typename Refs_::Sphere_circle Sphere_circle;
  /*{\Mtypemember circles.}*/
  typedef typename Refs_::Sphere_direction Sphere_direction;
  /*{\Mtypemember directions.}*/

  /*{\Mtext Local types are handles, iterators and circulators of the 
  following kind: |SVertex_const_handle|, |SVertex_const_iterator|, 
  |SHalfedge_const_handle|, |SHalfedge_const_iterator|, 
  |SHalfloop_const_handle|, |SHalfloop_const_iterator|, 
  |SFace_const_handle|, |SFace_const_iterator|.}*/

  typedef CGAL::Object_handle SObject_handle;
  /*{\Mtypemember a generic handle to an object of the underlying plane
  map. The kind of the object |(vertex, halfedge,face)| can be determined and
  the object assigned by the three functions:\\ 
  |bool assign(SVertex_const_handle& h, SObject_handle o)|\\ 
  |bool assign(SHalfedge_const_handle& h, SObject_handle o)|\\ 
  |bool assign(SFace_const_handle& h, SObject_handle o)|\\ where each 
  function returns |true| iff the assignment of |o| to |h| was valid.}*/

  #define USING(t) typedef typename Refs_::t t
  USING(Vertex_handle);   
  USING(SVertex_handle);   
  USING(SHalfedge_handle);   
  USING(SHalfloop_handle);   
  USING(SFace_handle);     
  USING(SVertex_const_handle); 
  USING(SHalfedge_const_handle); 
  USING(SHalfloop_const_handle); 
  USING(SFace_const_handle); 
  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SHalfloop_iterator);
  USING(SFace_iterator);
  USING(SVertex_const_iterator);
  USING(SHalfedge_const_iterator);
  USING(SHalfloop_const_iterator);
  USING(SFace_const_iterator);
  #undef USING
  #define DECUSING(t) typedef typename Base::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
  DECUSING(SHalfedge_around_sface_const_circulator);
  #undef DECUSING


  Sphere_segment segment(SHalfedge_const_handle e) const
  { return Sphere_segment(point(source(e)), point(target(e)), circle(e)); }

  Sphere_direction direction(SHalfedge_const_handle e) const
  { return Sphere_direction(circle(e)); }

  SHalfedge_const_handle out_wedge(SVertex_const_handle v, 
    const Sphere_direction& d, bool& collinear) const
  /*{\Xop returns a halfedge |e| bounding a wedge in between two
  neighbored edges in the adjacency list of |v| which contains |d|.
  If |d| extends along a edge then |e| is this edge. If |d| extends
  into the interior of such a wedge then |e| is the first edge hit
  when |d| is rotated clockwise. \precond |v| is not isolated.}*/
  { TRACEN("out_wedge "<<PH(v));
    CGAL_nef3_assertion(!is_isolated(v));
    collinear=false;
    Sphere_point p = point(v);
    SHalfedge_const_handle e_res = first_out_edge(v);
    Sphere_direction d_res = direction(e_res);
    SHalfedge_around_svertex_const_circulator el(e_res),ee(el);
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

  /*{\Mcreation 3}*/

  SNC_SM_point_locator() : Base() {}

  /*{\Moptions constref=yes}*/
  SNC_SM_point_locator(Vertex_handle v) : Base(v) {}
  /*{\Mcreate constructs a point locator working on the sphere map of |v|.}*/
  /*{\Moptions constref=no}*/
  /*{\Moperations 2.5 0.5}*/

  const Mark& mark(SObject_handle h) const
  /*{\Mop returns the mark associated to the object |h|.}*/
  { SVertex_const_handle v; 
    SHalfedge_const_handle e; 
    SFace_const_handle f;
    if ( assign(v,h) ) return mark(v);
    if ( assign(e,h) ) return mark(e);
    if ( assign(f,h) ) return mark(f);
    CGAL_nef3_assertion_msg(0,
    "PM_point_locator::mark: SObject_handle holds no object.");
  }

  enum SOLUTION { is_vertex_, is_edge_, is_loop_ };
  // enumeration for internal use
  
  SObject_handle locate(const Sphere_point& p) const
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of the underlying plane map |P| which contains the point |p =
  s.source()| in its relative interior. |s.target()| must be a point
  such that |s| intersects the $1$-skeleton of |P|.}*/
  { TRACEN("locate naivly "<<p);
    SVertex_const_iterator v;
    CGAL_nef3_forall_svertices(v,*this) {
      if ( p == point(v) ) {
	TRACEN( "  on point"); 
	return SObject_handle(v);
      }
    }

    SHalfedge_const_iterator e;
    CGAL_nef3_forall_sedges(e,*this) {
      if ( segment(e).has_on(p) ) {
	TRACEN( "  on segment"); 
	return SObject_handle(e);
      }
    }
    if ( has_loop() && circle(shalfloop()).has_on(p) ) {
      TRACEN( "  on loop"); 
      SHalfloop_const_handle l = shalfloop();
      return SObject_handle(l);
    }

    // now in face:

    if ( number_of_svertices() == 0 && ! has_loop() ) {
      TRACEN("  on unique face");
      SFace_const_handle f = sfaces_begin();
      return SObject_handle(f);
    }

    SVertex_const_handle v_res;
    SHalfedge_const_handle e_res;
    SHalfloop_const_handle l_res;
    SOLUTION solution;

    TRACEN("  on face...");
    Sphere_segment s; // we shorten the segment iteratively
    if ( has_loop() ) {
      Sphere_circle c(circle(shalfloop()),p); // orthogonal through p
      s = Sphere_segment(p,intersection(c,circle(shalfloop())));
      l_res = circle(shalfloop()).has_on_positive_side(p) ? 
	shalfloop() : twin(shalfloop());
      solution = is_loop_;
      TRACEN("has loop, initial ray "<<s);
    } else { // has vertices !
      CGAL_nef3_assertion( number_of_svertices()!=0 );
      SVertex_const_handle vt = svertices_begin();
      Sphere_point pvt = point(vt);
      if ( p != pvt.antipode() ) s = Sphere_segment(p,pvt);
      else s = Sphere_segment(p,pvt,Sphere_circle(p,pvt));
      /* to verify if it is necesary to determine exactly the sphere circle
	 of (p, pvt) segment based on the Halffacet plane information 
	 when they are antipodal */
      v_res = vt;
      solution = is_vertex_;
      TRACEN("has vertices, initial ray "<<s);
    }

    // s now initialized
    
    Sphere_direction dso(s.sphere_circle().opposite());
    Unique_hash_map<SHalfedge_const_handle,bool> visited(false);
    CGAL_nef3_forall_svertices(v,*this) {
      Sphere_point vp = point(v);
      if ( s.has_on(vp) ) {
        TRACEN(" location via vertex at "<<vp);
        s = Sphere_segment(p,vp,s.sphere_circle()); // we shrink the segment
        if ( is_isolated(v) ) {
          v_res = v; solution = is_vertex_;
        } else { // not isolated
          bool dummy;
          e_res = out_wedge(v,dso,dummy);
          SHalfedge_around_svertex_const_circulator el(e_res),ee(el);
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

    CGAL_nef3_forall_shalfedges(e,*this) {
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
        return SObject_handle(SFace_const_handle(face(e_res)));
      case is_loop_:
        return SObject_handle(SFace_const_handle(face(l_res)));
      case is_vertex_:
        return SObject_handle(SFace_const_handle(face(v_res)));
      default: CGAL_nef3_assertion_msg(0,"missing solution.");
    }
    return SObject_handle(); // never reached!
  }

  template <typename Object_predicate>
  SObject_handle ray_shoot(const Sphere_point& p, 
			  const Sphere_direction& d,
			  const Object_predicate& M) const
  /*{\Mop returns an |SObject_handle o| which can be converted to a
  |SVertex_const_handle|, |SHalfedge_const_handle|, |SFace_const_handle|
  |h| as described above.  The object predicate |M| has to have
  function operators \\ |bool operator() (const
  SVertex_/SHalfedge_/SHalfloop_/SFace_const_handle&)|.\\ The object
  returned is intersected by |d.circle()|, has minimal distance to
  |p|, and |M(h)| holds on the converted object. The operation returns
  the null handle |NULL| if the ray shoot along |s| does not hit any
  object |h| of |M| with |M(h)|.}*/
  { 
    Sphere_circle c(d.circle());
    Sphere_segment s;
    bool s_init(false);
    SObject_handle h = locate(p);
    SVertex_const_handle v; 
    SHalfedge_const_handle e; 
    SHalfloop_const_handle l; 
    SFace_const_handle f;
    if ( assign(v,h) && M(v) ||
         assign(e,h) && M(e) ||
	 assign(l,h) && M(l) ||
         assign(f,h) && M(f) ) return h;
    h = SObject_handle(); 
    TRACEN("not contained");
#if 0
    HASEN: s am anfang circle, ab wann segment ?
	   wo loop ?

    CGAL_nef3_forall_svertices (v,*this) {
      Point pv = point(v);
      if ( !(s_init && s.has_on(pv) ||
	    !s_init && c.has_on(pv)) ) continue;
      TRACEN("candidate "<<pv);
      if ( M(v) ) {
        h = SObject_handle(v);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
        continue;
      }
      // now we know that v is not marked but on s
      bool collinear;
      SHalfedge_const_handle e = out_wedge(v,d,collinear);
      if ( collinear ) { 
        if ( M(e) ) {
          h = SObject_handle(e);
          s = Sphere_segment(p,pv,c);
        }
        continue;
      }
      if ( M(face(e)) ) {
        h = SObject_handle(face(e));
        s = Sphere_segment(p,pv,c);
      }
    } // all vertices

    CGAL::Unique_hash_map<SHalfedge_const_handle,bool> visited(false);
    SHalfedge_const_iterator e_res;
    CGAL_nef3_forall_shalfedges(e,*this) {
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
	       strictly_ordered_ccw_at(p.antipode(),
				       direction(e),d,direction(twin(e))) )
	  e_res = twin(e);
        if ( M(e_res) ) {
          h = SObject_handle(e_res); s = s_cand;
        } else if ( M(face(twin(e_res))) ) {
          h = SObject_handle(face(twin(e_res))); s = s_cand;
        }
      }
    }
#endif
    CGAL_nef3_assertion_msg(0,"not yet correct");
    return h;
  }

  void init_marks_of_halfspheres();
  /*{\Mop initializes the default marks of the sphere map.}*/

  // C++ is really friendly:
  #define USECMARK(t) const Mark& mark(t h) const { return Base::mark(h); }
  #define USEMARK(t)  Mark& mark(t h) const { return Base::mark(h); }
  USEMARK(SVertex_handle)
  USEMARK(SHalfedge_handle)
  USEMARK(SHalfloop_handle)
  USEMARK(SFace_handle)
  USECMARK(SVertex_const_handle)
  USECMARK(SHalfedge_const_handle)
  USECMARK(SHalfloop_const_handle)
  USECMARK(SFace_const_handle)
  #undef USEMARK
  #undef USECMARK

  /*{\Mimplementation Naive query operations are realized by checking
  the intersection points of the $1$-skeleton of the plane map |P| with
  the query segments $s$. This method takes time linear in the size $n$
  of the underlying plane map without any preprocessing.}*/

}; // SNC_SM_point_locator<Refs_>



template <typename D>
void SNC_SM_point_locator<D>::init_marks_of_halfspheres()
{ TRACEN("init_marks_of_halfspheres");
  Sphere_point y_minus(0,-1,0);
  SObject_handle h = locate(y_minus);
  SFace_const_handle f;
  if ( CGAL::assign(f,h) ) { 
    TRACEN("on face ");
    mark_of_halfsphere(-1) = mark_of_halfsphere(+1) = mark(f);
    return;
  }

  SHalfedge_const_handle e;
  if ( CGAL::assign(e,h) ) { 
    CGAL_nef3_assertion(circle(e).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(e).orthogonal_vector());
    TRACEN("on edge "<<op);
    if ( (op.x() > 0) || (op.x() == 0) && (op.z() < 0) ) e = twin(e);
    // if ( (op.z() < 0) || (op.z() == 0) && (op.x() > 0) ) e = twin(e);
    mark_of_halfsphere(+1) = mark(face(e));
    mark_of_halfsphere(-1) = mark(face(twin(e)));
    return;
  }

  SHalfloop_const_handle l;
  if ( CGAL::assign(l,h) ) {
    CGAL_nef3_assertion(circle(l).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(l).orthogonal_vector());
    TRACEN("on loop "<<op);
    if ( (op.x() > 0) || (op.x() == 0) && (op.z() < 0) ) l = twin(l);
    // if ( (op.z() < 0) || (op.z() == 0) && (op.x() > 0) ) l = twin(l);
    mark_of_halfsphere(+1) = mark(face(l));
    mark_of_halfsphere(-1) = mark(face(twin(l)));
    return;
  }

  Sphere_circle c(0,0,1);
  Sphere_direction right(c),left(c.opposite());
  bool collinear(false);
  SVertex_const_handle v;
  if ( CGAL::assign(v,h) ) {
    CGAL_nef3_assertion(point(v)==y_minus);
    e = out_wedge(v,left,collinear); 
    if ( collinear ) mark_of_halfsphere(+1) = mark(face(twin(e)));
    else mark_of_halfsphere(+1) = mark(face(e));
    e = out_wedge(v,right,collinear); 
    if ( collinear ) mark_of_halfsphere(-1) = mark(face(twin(e)));
    else mark_of_halfsphere(-1) = mark(face(e));
    return;
  }
  /*
  TRACEN("1 dimensional object");
  mark_of_halfsphere(-1) = mark_of_halfsphere(+1) = 0;
  */
  CGAL_nef3_assertion_msg(0,"damn wrong type");
  return;
}

CGAL_END_NAMESPACE
#endif // CGAL_SNC_SM_POINT_LOCATOR_H



