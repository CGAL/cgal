// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_SM_POINT_LOCATOR_H
#define CGAL_SM_POINT_LOCATOR_H

#include <vector>
#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/SM_decorator_traits.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 47
#include <CGAL/Nef_2/debug.h>


namespace CGAL {

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

  using Base::cyclic_adj_succ;
  using Base::is_isolated;
  using Base::first_out_edge;

  Sphere_segment segment(SHalfedge_handle e) const
  { return Sphere_segment(e->source()->point(), e->twin()->source()->point(), e->circle()); }

  Sphere_direction direction(SHalfedge_handle e) const
  { return Sphere_direction(e->circle()); }

  SHalfedge_handle out_wedge(SVertex_handle v, const Sphere_direction& d, 
			     bool& collinear) const
  /*{\Xop returns a halfedge |e| bounding a wedge in between two
  neighbored edges in the adjacency list of |v| which contains |d|.
  If |d| extends along a edge then |e| is this edge. If |d| extends
  into the interior of such a wedge then |e| is the first edge hit
  when |d| is rotated clockwise. \precond |v| is not isolated.}*/
  { CGAL_NEF_TRACEN("out_wedge "<<PH(v));
    CGAL_assertion(!is_isolated(v));
    collinear=false;
    Sphere_point p = v->point();
    SHalfedge_handle e_res = first_out_edge(v);
    Sphere_direction d_res = direction(e_res);
    SHalfedge_around_svertex_circulator el(e_res),ee(el);
    if(direction(el) == d) {
      collinear = true;
      CGAL_NEF_TRACEN("  determined "<<PH(el) << el->circle());
      return el;
    }
    CGAL_For_all(el,ee) {
      if(direction(cyclic_adj_succ(el)) == d) {
	collinear = true;
	CGAL_NEF_TRACEN("  equal "<<PH(cyclic_adj_succ(el)) << cyclic_adj_succ(el)->circle());
	return cyclic_adj_succ(el);
      }
      else {
	CGAL_NEF_TRACEN("strictly_ordered_ccw " << direction(el) << " ? " << d << " ? " << direction(cyclic_adj_succ(el)));
	//	if ( strictly_ordered_ccw_at(p,d_res, direction(el), d) ) {
	if ( strictly_ordered_ccw_at(p,direction(el), d, direction(cyclic_adj_succ(el))) ) {
	  CGAL_NEF_TRACEN("strictly_ordered_ccw " << direction(el) << " - " << d << " - " << direction(cyclic_adj_succ(el)));
	  e_res = el; d_res = direction(e_res); break;
	}
      }
    }
    CGAL_NEF_TRACEN("  finally determined "<<PH(e_res) << e_res->circle());
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
    if ( CGAL::assign(v,h) ) return v->mark();
    if ( CGAL::assign(e,h) ) return e->mark();
    if ( CGAL::assign(l,h) ) return l->mark();
    CGAL_assertion_msg(CGAL::assign(f,h),
	 "PM_point_locator::mark: Object_handle holds no object.");
    CGAL::assign(f,h);
    return f->mark();
  }

  enum SOLUTION { is_vertex_, is_edge_, is_loop_ };
  // enumeration for internal use

  Object_handle locate(const Sphere_point& p, bool skipVEL = false)
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of the underlying plane map |P| which contains the point |p =
  s.source()| in its relative interior. |s.target()| must be a point
  such that |s| intersects the $1$-skeleton of |P|.}*/
  { CGAL_NEF_TRACEN("locate naivly "<< normalized(p));
    SVertex_iterator v;
    SHalfedge_iterator e;

    if(!skipVEL) {
      CGAL_forall_svertices(v,*this) {
	if ( p == v->point() ) {
	  CGAL_NEF_TRACEN( "  on point"); 
	  return make_object(v);
	}
      }
      
      CGAL_forall_sedges(e,*this) {
	if ( segment(e).has_on(p) || 
	     (e->source() == e->twin()->source() && e->circle().has_on(p))) {
	  CGAL_NEF_TRACEN( "  on segment " << segment(e));
	  return make_object(e);
	}
      }
      
      if ( this->has_shalfloop() && this->shalfloop()->circle().has_on(p)) {
	CGAL_NEF_TRACEN( "  on loop");
	return make_object(SHalfloop_handle(this->shalfloop()));
      }
    }
    

    // now in face:

    if(this->number_of_sfaces() == 1) {
      CGAL_NEF_TRACEN("  on unique face");
      SFace_handle f = this->sfaces_begin();
      return make_object(f);
    }

    SVertex_handle v_res;
    SHalfedge_handle e_res;
    SHalfloop_handle l_res(this->shalfloop());
    SOLUTION solution;

    CGAL_NEF_TRACEN("  on face...");
    Sphere_segment s; // we shorten the segment iteratively
    if ( this->has_shalfloop() ) {
      Sphere_circle c(this->shalfloop()->circle(),p); // orthogonal through p
      s = Sphere_segment(p,intersection(c,this->shalfloop()->circle()));
      l_res = this->shalfloop()->circle().has_on_positive_side(p) ? 
	this->shalfloop() : this->shalfloop()->twin();
      solution = is_loop_;
      CGAL_NEF_TRACEN("has loop, initial ray "<<s);
      CGAL_NEF_TRACEN(l_res->circle());
    } else { // has vertices !
      CGAL_assertion( this->number_of_svertices()!=0 );
      SVertex_iterator vi = this->svertices_begin();
      if( p == vi->point().antipode()) {
	++vi;
	CGAL_assertion( vi != this->svertices_end());
      }
      CGAL_NEF_TRACEN("initial segment: "<<p<<","<<vi->point());
      CGAL_assertion( p != vi->point().antipode());
      s = Sphere_segment( p, vi->point());
      v_res = vi;
      solution = is_vertex_;
      CGAL_NEF_TRACEN("has vertices, initial ray "<<s);
    }

    // s now initialized
    
    Sphere_direction dso(s.sphere_circle().opposite());
    Unique_hash_map<SHalfedge_handle,bool> visited(false);
    CGAL_forall_svertices(v,*this) {
      Sphere_point vp = v->point();
      if ( s.has_on(vp) ) {
        CGAL_NEF_TRACEN(" location via vertex at "<<vp);
        s = Sphere_segment(p,vp,s.sphere_circle()); // we shrink the segment
        if ( is_isolated(v) ) {
          CGAL_NEF_TRACEN("is_vertex_");
          v_res = v; solution = is_vertex_;
        } else { // not isolated
          bool dummy;
          e_res = out_wedge(v,dso,dummy);
          SHalfedge_around_svertex_circulator el(e_res),ee(el);
          CGAL_For_all(el,ee) 
            visited[el] = visited[el->twin()] = true;
          /* e_res is now the counterclockwise maximal halfedge out
             of v just before s */
          if ( e_res->circle().has_on_negative_side(p) )
            e_res = e_res->sprev();
            // correction to make e_res visible from p
	  solution = is_edge_;
          CGAL_NEF_TRACEN("  determined "<<PH(e_res));
        }
      }
    }

    CGAL_forall_sedges(e,*this) {
      if ( visited[e] ) continue;
      Sphere_segment se = segment(e);
      Sphere_point p_res;
      if(e->source() == e->twin()->source()) {
	Sphere_point p_res = intersection(e->circle(), s.sphere_circle());
	if(!s.has_in_relative_interior(p_res)) {
	  p_res = p_res.antipode();
	  if(!s.has_in_relative_interior(p_res))
	    continue;
	}
        s = Sphere_segment(p,p_res,s.sphere_circle()); 
        e_res = ( e->circle().has_on_positive_side(p) ? e : e->twin() );
        visited[e] = visited[e->twin()] = true;
	solution = is_edge_;
        CGAL_NEF_TRACEN("  determined "<<PH(e_res)<<" "<< e_res->incident_sface()->mark());	
      }
      else if ( do_intersect_internally(se,s,p_res) ) {
          CGAL_NEF_TRACEN(" location via halfedge "<<se);
        s = Sphere_segment(p,p_res,s.sphere_circle()); 
        e_res = ( e->circle().has_on_positive_side(p) ? e : e->twin() );
        visited[e] = visited[e->twin()] = true;
	solution = is_edge_;
        CGAL_NEF_TRACEN("  determined "<<PH(e_res)<<" "<< e_res->incident_sface()->mark());
      }
    }

    switch ( solution ) {
      case is_edge_: 
        return make_object(SFace_handle(e_res->incident_sface()));
      case is_loop_:
        return make_object(SFace_handle(l_res->incident_sface()));
      case is_vertex_:
        return make_object(SFace_handle(v_res->incident_sface()));
      default: CGAL_error_msg("missing solution.");
    }
    CGAL_error();
    return Object_handle(); // never reached!
  }

#if 0 //THIS CODE DOES NOT SEEM TO BE USED
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
    Object_handle h = locate(p);
    SVertex_handle v; 
    SHalfedge_handle e; 
    SHalfloop_handle l; 
    SFace_handle f;
    if ( ( CGAL::assign(v,h) && M(v) ) ||
         ( CGAL::assign(e,h) && M(e) ) ||
         ( CGAL::assign(l,h) && M(l) ) ||
         ( CGAL::assign(f,h) && M(f) ) ) return h;
    h = Object_handle(); 
    CGAL_NEF_TRACEN("not contained");
#if 0
    HASEN: s am anfang circle, ab wann segment ?
	   wo loop ?
    bool s_init(false);
    CGAL_forall_svertices (v,*this) {
      Point pv = v->point();
      if ( !(s_init && s.has_on(pv) ||
	    !s_init && c.has_on(pv)) ) continue;
      CGAL_NEF_TRACEN("candidate "<<pv);
      if ( M(v) ) {
        h = make_object(v);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
        continue;
      }
      // now we know that v is not marked but on s
      bool collinear;
      SHalfedge_handle e = out_wedge(v,d,collinear);
      if ( collinear ) { 
        if ( M(e) ) {
          h = make_object(e);
          s = Sphere_segment(p,pv,c);
        }
        continue;
      }
      if ( M(e->incident_sface()) ) {
        h = make_object(e->incident_sface());
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
        CGAL_NEF_TRACEN("candidate "<<se); 
        e_res = e;
        Sphere_segment s_cand = Sphere_segment(p,p_res,c);  
        if ( s_cand.is_short() && e->circle().has_on_negative_side(p) ||
	     s_cand.is_long() && e->circle().has_on_positive_side(p) ||
	     s_cand.is_halfcircle() && 
	       strictly_ordered_ccw_at(p.antipode(),
				       direction(e),d,direction(e->twin())) )
	  e_res = e->twin();
        if ( M(e_res) ) {
          h = make_object(e_res); s = s_cand;
        } else if ( M(face(twin(e_res))) ) {
          h = make_object(face(twin(e_res))); s = s_cand;
        }
      }
    }
#endif
    CGAL_error_msg("not yet correct");
    return h;
  }
#endif
  
  Object_handle ray_shoot(const Sphere_point& p, 
			  const Sphere_circle& c,
			  Sphere_point& ip,
			  bool start_inclusive = false) { 
    Sphere_segment seg(p, p.antipode(), c);
    return ray_shoot(seg, ip, start_inclusive);
  }

  Object_handle ray_shoot(const Sphere_segment& d, 
			  Sphere_point& ip,
			  bool start_inclusive = false,
			  bool beyond_end = true,
			  bool end_inclusive = false) { 

    // TODO: end_inclusive=true does not work properly for sedges and sloops

    CGAL_NEF_TRACEN("ray shoot");
    Sphere_circle c(d.sphere_circle());
    Sphere_point p(d.source());
    Sphere_segment s;
    bool s_init(false);
    
    if(!beyond_end) {
      s = d;
      s_init = true;
    }

    Object_handle h = Object_handle();

    if(s_init) {
      CGAL_NEF_TRACEN(" at begin " << s_init << ":" << s);
    } else {
      CGAL_NEF_TRACEN(" at begin " << s_init << ":" << c);
    }

    SVertex_iterator vi;
    CGAL_forall_svertices (vi,*this) {
      Sphere_point pv = vi->point();
      if ((s_init && !s.has_on(pv)) ||
	  (!s_init && !c.has_on(pv))) continue;
      CGAL_NEF_TRACEN("candidate "<<pv);
      CGAL_NEF_TRACEN("p =?= pv: " << p << ", " << pv); 
      if ((start_inclusive || p != pv) && 
	  (end_inclusive || !s_init || s.target() != pv)) {
        h = make_object(vi);     // store vertex
        s = Sphere_segment(p,pv,c); // shorten
	ip = pv;
	s_init = true;
      }
    }

    // TODO: edges on the ray.

    SHalfedge_iterator ei;
    CGAL_forall_sedges(ei,*this) {
      Sphere_segment se = segment(ei);
      CGAL_NEF_TRACEN("ray_shoot " << s_init);
      if(s_init) {
	CGAL_NEF_TRACEN("  " << s.source() << "->" << s.target() << " | " << s.sphere_circle() << " is long " << s.is_long());
      }
      CGAL_NEF_TRACEN("  " << se.source() << "->" << se.target() << " | " << se.sphere_circle() << " is long " << se.is_long());

      // TODO: start-end point of s on se or c

      Sphere_point p_res;
      if(se.source() == se.target()) {

	if(s_init) {
	  if(s.is_long()) {
	    Sphere_segment first_half(p,p.antipode(),c);
	    Sphere_segment second_part(p.antipode(), s.target(), c);
	    if(!do_intersect_internally(ei->circle(), first_half, p_res) &&
	       !do_intersect_internally(ei->circle(), second_part, p_res)) {
	      if(se.has_on(p.antipode())) {
		p_res = p.antipode();
	      } else {
		continue;
              }
            }
	  } else {
	    if(!do_intersect_internally(ei->circle(), s, p_res)) continue;
	  }
	} else {
	  Sphere_segment first_half(p,p.antipode(),c);
	  Sphere_segment second_part(p.antipode(),p,c);
	  if(!do_intersect_internally(ei->circle(), first_half, p_res)) {
	    do_intersect_internally(ei->circle(), second_part, p_res);
	  }
	}

      } else {
	if(s_init) {	  
	  if(s.is_long() && se.is_long()) {
	    Sphere_segment first_half(p,p.antipode(),c);
	    Sphere_segment second_part(p.antipode(), s.target(), c);
	    if(!do_intersect_internally(se, first_half, p_res) &&
	       !do_intersect_internally(se, second_part, p_res)) {
	      if(se.has_on(p.antipode()))
		p_res = p.antipode();
	      else
		continue;
	    }	    
	  } else {
	    if(!do_intersect_internally(se, s, p_res)) continue;
	  }  
	} else {
	  if(se.is_long()) {
	    Sphere_segment first_half(p,p.antipode(),c);
	    Sphere_segment second_half(p.antipode(),p,c); 
	    if(!do_intersect_internally(se, first_half, p_res)) {
	      if(!do_intersect_internally(se, second_half, p_res)) {
		if(start_inclusive)
		  p_res = p;
		else
		  p_res = p.antipode();
	      }
	    }
	  } else {
	    if(!do_intersect_internally(c, se, p_res)) continue;
	  }
	}
      }
      
      CGAL_NEF_TRACEN("candidate "<<se); 
      if (start_inclusive || p != p_res) {
	h = make_object(ei); 
	s = Sphere_segment(p,p_res,c);
	ip = p_res;
	s_init = true;
      }
    }

    // TODO: start-end point of s on cl

    if(this->has_shalfloop()) {
      Sphere_circle cl(this->shalfloop()->circle());
      if(!s_init || s.is_long()) {
	if(cl.has_on(p)) {
	  ip = p.antipode();
	  return make_object(SHalfloop_handle(this->shalfloop()));
	} else 	  
	  s = Sphere_segment(p,p.antipode(),c);
      }
      Sphere_point p_res;
      CGAL_NEF_TRACEN("do intersect " << cl << ", " << s);
      if(!do_intersect_internally(cl,s,p_res))
	return h;
      /*
      if(p_res == p.antipode()) // does this happen ? test has_on for p/p.antipode ?
	p_res = p;
      */
      CGAL_NEF_TRACEN("found intersection point " << p_res);
      CGAL_assertion_code(Sphere_segment testseg(p,p_res,c));
      CGAL_assertion(!testseg.is_long());
      if (start_inclusive || p != p_res) {
	ip = p_res;
	return make_object(SHalfloop_handle(this->shalfloop()));
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

  CGAL_NEF_TRACEN("marks_of_halfspheres ");

  Sphere_point y_minus;
  if(axis!=1) 
    y_minus = Sphere_point(0,-1,0);
  else
    y_minus = Sphere_point(0,0,1);
  Object_handle h = locate(y_minus);
  SFace_handle f;
  if ( CGAL::assign(f,h) ) { 
    CGAL_NEF_TRACEN("on face " << mark(make_object(f)));
    lower = upper = mark(make_object(f));
    return;
  }

  SHalfedge_handle e;
  if ( CGAL::assign(e,h) ) { 
    CGAL_assertion(e->circle().has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+e->circle().orthogonal_vector());
    CGAL_NEF_TRACEN("on edge "<<op);
    if (axis==0 && ((op.z() < 0) || ((op.z() == 0) && (op.x() < 0)))) e = e->twin();
    if (axis==1 && ((op.x() > 0) || ((op.x() == 0) && (op.y() < 0)))) e = e->twin();
    if (axis==2 && ((op.x() > 0) || ((op.x() == 0) && (op.z() < 0)))) e = e->twin();
    upper = e->incident_sface()->mark();
    lower = e->twin()->incident_sface()->mark();
    return;
  }

  SHalfloop_handle l;
  if ( CGAL::assign(l,h) ) {
    CGAL_assertion(l->circle().has_on(y_minus));
    Sphere_point op(CGAL::ORIGIN+l->circle().orthogonal_vector());
    CGAL_NEF_TRACEN("on loop "<<op);
    if (axis==0 && ((op.z() < 0) || ((op.z() == 0) && (op.x() < 0)))) l = l->twin();
    if (axis==1 && ((op.x() > 0) || ((op.x() == 0) && (op.y() < 0)))) l = l->twin();
    if (axis==2 && ((op.x() > 0) || ((op.x() == 0) && (op.z() < 0)))) l = l->twin();
    upper = l->incident_sface()->mark();
    lower = l->twin()->incident_sface()->mark();
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
    CGAL_assertion(v->point()==y_minus);
    if(is_isolated(v))
      upper = lower = mark(make_object(v->incident_sface()));
    else {
      e = out_wedge(v,left,collinear); 
      if ( collinear ) upper = e->twin()->incident_sface()->mark();
      else upper = e->incident_sface()->mark();
      e = out_wedge(v,right,collinear); 
      if ( collinear ) lower = e->twin()->incident_sface()->mark();
      else lower = e->incident_sface()->mark();
    }
  } else {
    CGAL_error();
  }
}

} //namespace CGAL
#endif // CGAL_SM_POINT_LOCATOR_H
