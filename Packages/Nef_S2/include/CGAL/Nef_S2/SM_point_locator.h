// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_SM_POINT_LOCATOR_H
#define CGAL_SM_POINT_LOCATOR_H

#include <vector>
#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/SM_decorator_traits.h>
#undef _DEBUG
#define _DEBUG 47
#include <CGAL/Nef_S2/debug.h>


CGAL_BEGIN_NAMESPACE

/*{\Moptions print_title=yes }*/ 

/*{\Manpage {SM_point_locator}{Decorator}
{Naive point location in plane maps}{PL}}*/
/*{\Mdefinition An instance |\Mvar| of data type |\Mname|
encapsulates naive point location queries within a sphere map |M|.  The
two template parameters are specified via concepts. |SM_decorator_|
must be a model of the concept |SMDecorator| as described in the
appendix.  |Geometry_| must be a model of the concept
|AffineGeometryTraits_2| as described in the appendix. For a
specification of plane maps see also the concept of
|SMConstDecorator|.}*/

/*{\Mgeneralization Decorator}*/

template <class SM_decorator>
class SM_point_locator : public SM_decorator {
protected:
  typedef SM_decorator                                  Base;
  typedef SM_point_locator<SM_decorator>                Self;
  typedef typename SM_decorator::Sphere_map             Sphere_map;

public:
  /*{\Mtypes 5}*/

  typedef typename SM_decorator::Decorator_traits  Decorator_traits;

  typedef typename Base::Mark      Mark;
  /*{\Mtypemember the attribute of all objects (vertices, edges, loops,
  faces).}*/

  typedef typename SM_decorator::Sphere_kernel Sphere_kernel;
  /*{\Mtypemember the sphere kernel.}*/
  typedef typename Sphere_kernel::Sphere_point Sphere_point;
  /*{\Mtypemember points.}*/
  typedef typename Sphere_kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember segments.}*/
  typedef typename Sphere_kernel::Sphere_circle Sphere_circle;
  /*{\Mtypemember circles.}*/
  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;
  /*{\Mtypemember directions.}*/

  /*{\Mtext Local types are handles, iterators and circulators of the 
  following kind: |SVertex_const_handle|, |SVertex_const_iterator|, 
  |SHalfedge_const_handle|, |SHalfedge_const_iterator|, 
  |SHalfloop_const_handle|, |SHalfloop_const_iterator|, 
  |SFace_const_handle|, |SFace_const_iterator|.}*/

  typedef CGAL::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying plane
  map. The kind of the object |(vertex, halfedge,face)| can be determined and
  the object assigned by the three functions:\\ 
  |bool assign(SVertex_const_handle& h, Object_handle o)|\\ 
  |bool assign(SHalfedge_const_handle& h, Object_handle o)|\\ 
  |bool assign(SFace_const_handle& h, Object_handle o)|\\ where each 
  function returns |true| iff the assignment of |o| to |h| was valid.}*/

  typedef typename Decorator_traits::SVertex_handle       SVertex_handle;   
  typedef typename Decorator_traits::SHalfedge_handle     SHalfedge_handle;   
  typedef typename Decorator_traits::SHalfloop_handle     SHalfloop_handle;   
  typedef typename Decorator_traits::SFace_handle         SFace_handle;     

  typedef typename Decorator_traits::SVertex_iterator     SVertex_iterator;
  typedef typename Decorator_traits::SHalfedge_iterator   SHalfedge_iterator;
  typedef typename Decorator_traits::SHalfloop_iterator   SHalfloop_iterator;
  typedef typename Decorator_traits::SFace_iterator       SFace_iterator;

  typedef typename Decorator_traits::SHalfedge_around_svertex_circulator
                                     SHalfedge_around_svertex_circulator;
  typedef typename Decorator_traits::SHalfedge_around_sface_circulator
                                     SHalfedge_around_sface_circulator;

  Sphere_segment segment(SHalfedge_handle e) const
  { return Sphere_segment(point(source(e)), point(target(e)), circle(e)); }

  Sphere_direction direction(SHalfedge_handle e) const
  { return Sphere_direction(circle(e)); }

  SHalfedge_handle out_wedge(SVertex_handle v, const Sphere_direction& d, 
			     bool& collinear) const
  /*{\Xop returns a halfedge |e| bounding a wedge in between two
  neighbored edges in the adjacency list of |v| which contains |d|.
  If |d| extends along a edge then |e| is this edge. If |d| extends
  into the interior of such a wedge then |e| is the first edge hit
  when |d| is rotated clockwise. \precond |v| is not isolated.}*/
  { TRACEN("out_wedge "<<PH(v));
    CGAL_assertion(!is_isolated(v));
    collinear=false;
    Sphere_point p = point(v);
    SHalfedge_handle e_res = first_out_edge(v);
    Sphere_direction d_res = direction(e_res);
    SHalfedge_around_svertex_circulator el(e_res),ee(el);
    if(direction(el) == d) {
      collinear = true;
      TRACEN("  determined "<<PH(el) << circle(el));
      return el;
    }
    CGAL_For_all(el,ee) {
      if(direction(cyclic_adj_succ(el)) == d) {
	collinear = true;
	TRACEN("  equal "<<PH(cyclic_adj_succ(el)) << circle(cyclic_adj_succ(el)));
	return cyclic_adj_succ(el);
      }
      else {
	TRACEN("strictly_ordered_ccw " << direction(el) << " ? " << d << " ? " << direction(cyclic_adj_succ(el)));
	//	if ( strictly_ordered_ccw_at(p,d_res, direction(el), d) ) {
	if ( strictly_ordered_ccw_at(p,direction(el), d, direction(cyclic_adj_succ(el))) ) {
	  TRACEN("strictly_ordered_ccw " << direction(el) << " - " << d << " - " << direction(cyclic_adj_succ(el)));
	  e_res = el; d_res = direction(e_res); break;
	}
      }
    }
    TRACEN("  finally determined "<<PH(e_res) << circle(e_res));
    /*
    Sphere_direction d2 = direction(cyclic_adj_succ(e_res));
    d2 = normalized(d2);
    TRACEN(d2 << " =?= " << d);
    if (d2 == d ) {
      e_res = cyclic_adj_succ(e_res);
      collinear=true;
    }
    TRACEN("  wedge = "<<PH(e_res)<<" "<<collinear);
    */
    return e_res;
  }

  /*{\Mcreation 3}*/

  SM_point_locator() : Base() {}

  /*{\Moptions constref=yes}*/
  SM_point_locator(Sphere_map* cp) : Base(cp) {}
  /*{\Mcreate constructs a point locator working on |P|.}*/
  /*{\Moptions constref=no}*/
  /*{\Moperations 2.5 0.5}*/

  const Mark& mark(Object_handle h) const
  /*{\Mop returns the mark associated to the object |h|.}*/
  { SVertex_handle v; 
    SHalfedge_handle e; 
    SHalfloop_handle l;
    SFace_handle f;
    if ( CGAL::assign(v,h) ) return Base::mark(v);
    if ( CGAL::assign(e,h) ) return Base::mark(e);
    if ( CGAL::assign(l,h) ) return Base::mark(l);
    CGAL_assertion_msg(CGAL::assign(f,h),
	 "PM_point_locator::mark: Object_handle holds no object.");
    CGAL::assign(f,h);
    return Base::mark(f);
  }

  enum SOLUTION { is_vertex_, is_edge_, is_loop_ };
  // enumeration for internal use

  Object_handle locate(const Sphere_point& p) const
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of the underlying plane map |P| which contains the point |p =
  s.source()| in its relative interior. |s.target()| must be a point
  such that |s| intersects the $1$-skeleton of |P|.}*/
  { TRACEN("locate naivly "<<p);
    SVertex_iterator v;
    CGAL_forall_svertices(v,*this) {
      if ( p == point(v) ) {
	TRACEN( "  on point"); 
	return Object_handle(v);
      }
    }

    SHalfedge_iterator e;
    CGAL_forall_sedges(e,*this) {
      if ( segment(e).has_on(p) ) {
	TRACEN( "  on segment " << segment(e));
	return Object_handle(e);
      }
    }

    if ( this->has_shalfloop() && circle(this->shalfloop()).has_on(p) ) {
      TRACEN( "  on loop");
      return Object_handle(SHalfloop_handle(this->shalfloop()));
    }

    // now in face:

    if(this->number_of_sfaces() == 1) {
      TRACEN("  on unique face");
      SFace_handle f = this->sfaces_begin();
      return Object_handle(f);
    }

    SVertex_handle v_res;
    SHalfedge_handle e_res;
    SHalfloop_handle l_res(this->shalfloop());
    SOLUTION solution;

    TRACEN("  on face...");
    Sphere_segment s; // we shorten the segment iteratively
    if ( this->has_shalfloop() ) {
      Sphere_circle c(circle(this->shalfloop()),p); // orthogonal through p
      s = Sphere_segment(p,intersection(c,circle(this->shalfloop())));
      l_res = circle(this->shalfloop()).has_on_positive_side(p) ? 
	this->shalfloop() : twin(this->shalfloop());
      solution = is_loop_;
      TRACEN("has loop, initial ray "<<s);
      TRACEN(circle(l_res));
    } else { // has vertices !
      CGAL_assertion( this->number_of_svertices()!=0 );
      SVertex_iterator vi = this->svertices_begin();
      if( p == point(vi).antipode()) {
	++vi;
	CGAL_assertion( vi != this->svertices_end());
      }
      TRACEN("initial segment: "<<p<<","<<point(vi));
      CGAL_assertion( p != point(vi).antipode());
      s = Sphere_segment( p, point(vi));
      v_res = vi;
      solution = is_vertex_;
      TRACEN("has vertices, initial ray "<<s);
    }

    // s now initialized
    
    Sphere_direction dso(s.sphere_circle().opposite());
    Unique_hash_map<SHalfedge_handle,bool> visited(false);
    CGAL_forall_svertices(v,*this) {
      Sphere_point vp = point(v);
      if ( s.has_on(vp) ) {
        TRACEN(" location via vertex at "<<vp);
        s = Sphere_segment(p,vp,s.sphere_circle()); // we shrink the segment
        if ( is_isolated(v) ) {
          TRACEN("is_vertex_");
          v_res = v; solution = is_vertex_;
        } else { // not isolated
          bool dummy;
          e_res = out_wedge(v,dso,dummy);
          SHalfedge_around_svertex_circulator el(e_res),ee(el);
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

    CGAL_forall_sedges(e,*this) {
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
        return Object_handle(SFace_handle(face(e_res)));
      case is_loop_:
        return Object_handle(SFace_handle(face(l_res)));
      case is_vertex_:
        return Object_handle(SFace_handle(face(v_res)));
      default: CGAL_assertion_msg(0,"missing solution.");
    }
    return Object_handle(); // never reached!
  }

  template <typename Object_predicate>
  Object_handle ray_shoot(const Sphere_point& p, 
			  const Sphere_direction& d,
			  const Object_predicate& M) const
  /*{\Mop returns an |Object_handle o| which can be converted to a
  |SVertex_handle|, |SHalfedge_handle|, |SFace_handle|
  |h| as described above.  The object predicate |M| has to have
  function operators \\ |bool operator() (const
  SVertex_/SHalfedge_/SHalfloop_/SFace_handle&)|.\\ The object
  returned is intersected by |d.circle()|, has minimal distance to
  |p|, and |M(h)| holds on the converted object. The operation returns
  the null handle |NULL| if the ray shoot along |s| does not hit any
  object |h| of |M| with |M(h)|.}*/
  { 
    Sphere_circle c(d.circle());
    Sphere_segment s;
    bool s_init(false);
    Object_handle h = locate(p);
    SVertex_handle v; 
    SHalfedge_handle e; 
    SHalfloop_handle l; 
    SFace_handle f;
    if ( CGAL::assign(v,h) && M(v) ||
         CGAL::assign(e,h) && M(e) ||
	 CGAL::assign(l,h) && M(l) ||
         CGAL::assign(f,h) && M(f) ) return h;
    h = Object_handle(); 
    TRACEN("not contained");
#if 0
    HASEN: s am anfang circle, ab wann segment ?
	   wo loop ?

    CGAL_forall_svertices (v,*this) {
      Point pv = point(v);
      if ( !(s_init && s.has_on(pv) ||
	    !s_init && c.has_on(pv)) ) continue;
      TRACEN("candidate "<<pv);
      if ( M(v) ) {
        h = Object_handle(v);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
        continue;
      }
      // now we know that v is not marked but on s
      bool collinear;
      SHalfedge_handle e = out_wedge(v,d,collinear);
      if ( collinear ) { 
        if ( M(e) ) {
          h = Object_handle(e);
          s = Sphere_segment(p,pv,c);
        }
        continue;
      }
      if ( M(face(e)) ) {
        h = Object_handle(face(e));
        s = Sphere_segment(p,pv,c);
      }
    } // all vertices

    CGAL::Unique_hash_map<SHalfedge_handle,bool> visited(false);
    SHalfedge_iterator e_res;
    CGAL_forall_sedges(e,*this) {
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
          h = Object_handle(e_res); s = s_cand;
        } else if ( M(face(twin(e_res))) ) {
          h = Object_handle(face(twin(e_res))); s = s_cand;
        }
      }
    }
#endif
    CGAL_assertion_msg(0,"not yet correct");
    return h;
  }

  Object_handle ray_shoot(const Sphere_point& p, 
			  const Sphere_circle& c,
			  Sphere_point& ip,
			  bool start_inclusive = false) const { 
    //    Sphere_circle c(d.circle());
    CGAL_assertion(c.has_on(p));
    Sphere_segment s;
    bool s_init(false);
    Object_handle h = Object_handle();

    SVertex_iterator vi;
    CGAL_forall_svertices (vi,*this) {
      Sphere_point pv = point(vi);
      if (!(s_init && s.has_on(pv)) ||
	  (!s_init && c.has_on(pv))) continue;
      TRACEN("candidate "<<pv);
      if (start_inclusive || p != pv) {
        h = Object_handle(*vi);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
	ip = pv;
	s_init = true;
      }
    }
 
    SHalfedge_iterator ei;
    CGAL_forall_sedges(ei,*this) {
      Sphere_segment se = segment(ei);
      Sphere_point p_res;
      if (!(s_init && do_intersect_internally(se,s,p_res)) ||
	  (!s_init && do_intersect_internally(c,se,p_res))) continue;
      TRACEN("candidate "<<se); 
      if (start_inclusive || p != p_res) {
	h = Object_handle(*ei); 
	s = Sphere_segment(p,p_res,c);
	ip = p_res;
	s_init = true;
      }
    }

    if(this->has_shalfloop()) {
      Sphere_circle cl(this->shalfloop()->circle());
      if(!s_init)
	s = Sphere_segment(p,p.antipode(),c);
      Sphere_point p_res;
      if(!do_intersect_internally(cl,s,p_res))
	return h;
      if(p_res == p.antipode())
	p_res = p;
      if (start_inclusive || p != p_res) {
	ip = p_res;
	return Object_handle(this->shalfloop());
      }
    }

    return h;
  }

  void marks_of_halfspheres(std::vector<Mark>& mohs, int offset, 
			    int axis=2);
  void marks_of_halfspheres(Mark& unten, Mark& oben, int axis=2);

  /*{\Mimplementation Naive query operations are realized by checking
  the intersection points of the $1$-skeleton of the plane map |P| with
  the query segments $s$. This method takes time linear in the size $n$
  of the underlying plane map without any preprocessing.}*/
}; // SM_point_locator<SM_decorator>

template <typename D>
void SM_point_locator<D>::
marks_of_halfspheres(std::vector<Mark>& mohs, int offset, int axis) {
  Mark lower, upper;
  marks_of_halfspheres(lower, upper, axis);
  mohs[offset] = lower;
  mohs[offset+1] = upper;
}

template <typename D>
void SM_point_locator<D>::
marks_of_halfspheres(Mark& lower, Mark& upper, int axis) {

  TRACEN("marks_of_halfspheres ");

  Sphere_point y_minus;
  if(axis!=1) 
    y_minus = Sphere_point(0,-1,0);
  else
    y_minus = Sphere_point(0,0,1);
  Object_handle h = locate(y_minus);
  SFace_handle f;
  if ( CGAL::assign(f,h) ) { 
    TRACEN("on face " << mark(f));
    lower = upper = mark(f);
    return;
  }

  SHalfedge_handle e;
  if ( CGAL::assign(e,h) ) { 
    CGAL_assertion(circle(e).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(e).orthogonal_vector());
    TRACEN("on edge "<<op);
    if (axis==0 && ((op.z() < 0) || (op.z() == 0) && (op.x() < 0))) e = twin(e);
    if (axis==1 && ((op.x() > 0) || (op.x() == 0) && (op.y() < 0))) e = twin(e);
    if (axis==2 && ((op.x() > 0) || (op.x() == 0) && (op.z() < 0))) e = twin(e);
    upper = mark(face(e));
    lower = mark(face(twin(e)));
    return;
  }

  SHalfloop_handle l;
  if ( CGAL::assign(l,h) ) {
    CGAL_assertion(circle(l).has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+circle(l).orthogonal_vector());
    TRACEN("on loop "<<op);
    if (axis==0 && ((op.z() < 0) || (op.z() == 0) && (op.x() < 0))) l = twin(l);
    if (axis==1 && ((op.x() > 0) || (op.x() == 0) && (op.y() < 0))) l = twin(l);
    if (axis==2 && ((op.x() > 0) || (op.x() == 0) && (op.z() < 0))) l = twin(l);
    upper = mark(face(l));
    lower = mark(face(twin(l)));
    return;
  }

  Sphere_circle c;
  switch(axis) {
  case 0: c = Sphere_circle(1,0,0); break;
  case 1: c = Sphere_circle(0,1,0); break;
  case 2: c = Sphere_circle(0,0,1); break;
  }
  Sphere_direction right(c),left(c.opposite());
  bool collinear(false);
  SVertex_handle v;
  if ( CGAL::assign(v,h) ) {
    CGAL_assertion(point(v)==y_minus);
    if(is_isolated(v))
      upper = lower = mark(face(v));
    else {
      e = out_wedge(v,left,collinear); 
      if ( collinear ) upper = mark(face(twin(e)));
      else upper = mark(face(e));
      e = out_wedge(v,right,collinear); 
      if ( collinear ) lower = mark(face(twin(e)));
      else lower = mark(face(e));
    }
  }
}

CGAL_END_NAMESPACE
#endif // CGAL_SM_POINT_LOCATOR_H



