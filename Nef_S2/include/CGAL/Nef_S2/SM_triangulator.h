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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_SM_TRIANGULATOR_H
#define CGAL_SM_TRIANGULATOR_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_S2/SM_constrained_triang_traits.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 137
#include <CGAL/Nef_2/debug.h>

#define CGAL_USING(t) typedef typename Base::t t
#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif
namespace CGAL {

template <typename Decorator_, typename IT, typename INFO>
struct SM_subdivision {
  typedef Decorator_ Triangulator;
  typedef typename Decorator_::SVertex_handle Vertex_handle;
  typedef typename Decorator_::SHalfedge_handle   Halfedge_handle;
  typedef typename Decorator_::Sphere_point   Point;
  typedef typename Decorator_::Sphere_segment Segment;
  Triangulator T;
  CGAL::Unique_hash_map<IT,INFO>& M;
  /* M stores the object that supports the segment that
     is input object of the sweep */

  SM_subdivision(Triangulator Ti, 
                 CGAL::Unique_hash_map<IT,INFO>& Mi) : T(Ti), M(Mi) {}

Vertex_handle new_vertex(const Point& p)
{ Vertex_handle v = T.new_svertex(p); T.assoc_info(v);
  return v;
}

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e)
{ T.link_as_target_and_append(v,e); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
{ Halfedge_handle e = 
  T.new_shalfedge_pair_at_source(v,Decorator_::BEFORE); 
  T.assoc_info(e);
  return e;
}

void halfedge_below(Vertex_handle v, Halfedge_handle e) const 
{ T.halfedge_below(v) = e; }

/* the following operation associates segment support with
   halfedges, we only update if non-NULL; this prevents 
   artificial sphere subdivision segments that have NULL 
   support to overwrite non-NULL support */

void supporting_segment(Halfedge_handle e, IT it) const
{ T.is_forward(e) = true; 
  if ( ! M[it].empty() ) T.support(e) = M[it]; }

/* the following operation associate segment support with
   vertices, we only update if non-NULL; this prevents 
   artificial segments that have NULL support to overwrite
   non-NULL support */

void trivial_segment(Vertex_handle v, IT it) const
{ if ( ! M[it].empty() ) T.support(v) = M[it]; }

void starting_segment(Vertex_handle v, IT it) const
{ if ( ! M[it].empty() ) T.support(v) = M[it]; }

void ending_segment(Vertex_handle v, IT it) const
{ if ( ! M[it].empty() ) T.support(v) = M[it]; }

void passing_segment(Vertex_handle v, IT it) const
{ if ( ! M[it].empty() ) T.support(v) = M[it]; }


}; // SM_subdivision



/*{\Manpage {SM_triangulator}{Decorator_}{Overlay in the sphere}{O}}*/

template <typename Decorator_>
class SM_triangulator : public Decorator_ {
public:
  /*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
  decorator object offering sphere map triangulation calculation.}*/

  typedef Decorator_                            Base;
  typedef typename Decorator_::Sphere_map       Sphere_map;
  typedef CGAL::SM_const_decorator<Sphere_map>  SM_const_decorator;
  typedef SM_const_decorator                    Explorer;
  typedef Decorator_                            Decorator;
  typedef SM_triangulator<Decorator_>           Self;
  typedef CGAL::SM_point_locator<SM_const_decorator>  SM_point_locator;

  typedef typename SM_const_decorator::SVertex_const_handle SVertex_const_handle;
  typedef typename SM_const_decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SM_const_decorator::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename SM_const_decorator::SFace_const_handle SFace_const_handle;
  typedef typename SM_const_decorator::SVertex_const_iterator SVertex_const_iterator;
  typedef typename SM_const_decorator::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename SM_const_decorator::SFace_const_iterator SFace_const_iterator;

  typedef typename Base::SVertex_handle SVertex_handle;
  typedef typename Base::SHalfedge_handle SHalfedge_handle;
  typedef typename Base::SHalfloop_handle SHalfloop_handle;
  typedef typename Base::SFace_handle SFace_handle;
  typedef typename Base::SVertex_iterator SVertex_iterator;
  typedef typename Base::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Base::SFace_iterator SFace_iterator;
  typedef typename Base::Object_handle Object_handle;

  typedef typename Base::SHalfedge_around_svertex_circulator 
                         SHalfedge_around_svertex_circulator;
  typedef typename Base::SHalfedge_around_sface_circulator 
                         SHalfedge_around_sface_circulator;

  typedef std::pair<SHalfedge_handle,SHalfedge_handle> SHalfedge_pair;

  /*{\Mtypes 3}*/

  typedef typename Base::Sphere_kernel Sphere_kernel;

  typedef typename Sphere_kernel::Sphere_point   Sphere_point;
  /*{\Mtypemember the point type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember the segment type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_circle  Sphere_circle;
  /*{\Mtypemember the circle type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_triangle Sphere_triangle;
  /*{\Mtypemember the triangle type of the sphere geometry.}*/

  typedef typename Decorator::Mark        Mark;
  /*{\Mtypemember the mark of sphere map objects.}*/

  /*{\Mgeneralization Decorator_}*/

protected:
  const Explorer* E_;
  const Sphere_kernel& K;

public:

  typedef std::list<Sphere_segment>            Seg_list;
  typedef typename Seg_list::iterator          Seg_iterator;
  typedef std::pair<Seg_iterator,Seg_iterator> Seg_it_pair;
  typedef std::pair<Sphere_segment,Sphere_segment> Seg_pair;
  typedef CGAL::Unique_hash_map<Seg_iterator,Object_handle> Seg_map;

  using Base::info;
  using Base::set_first_out_edge;
  using Base::first_out_edge;
  using Base::last_out_edge;
  using Base::is_isolated;
  using Base::has_outdeg_two;
  using Base::flip_diagonal;
  using Base::out_edges;
  using Base::set_source;
  using Base::set_face;
  using Base::link_as_prev_next_pair;
  using Base::delete_vertex_only;
  using Base::delete_edge_pair_only;
  using Base::is_sm_boundary_object;
  using Base::store_sm_boundary_object;
  using Base::undo_sm_boundary_object;

  // vertex_info stores the origin of vertices
  struct vertex_info {
    Object_handle         o_;
    SHalfedge_handle       e_;
    vertex_info() : o_(),e_() {}
    LEDA_MEMORY(vertex_info)
  };

  void assoc_info(SVertex_handle v) const
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<vertex_info>::create(info(v)); 
    #else
    info(v)=vertex_info();
    #endif
  }

  void discard_info(SVertex_handle v) const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<vertex_info>::clear(info(v)); 
    #else
    info(v)=boost::any();
    #endif
  }

  vertex_info& ginfo(SVertex_handle v) const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<vertex_info>::access(info(v)); 
    #else
    return
      *boost::any_cast<vertex_info>(&info(v));
    #endif
  }

  Object_handle& support(SVertex_handle v) const
  { return ginfo(v).o_; }

  SHalfedge_handle& halfedge_below(SVertex_handle v) const
  { return ginfo(v).e_; }

  // edge_info stores the origin of edges
  struct edge_info {
    Mark m_left_; Object_handle o_; bool forw_;

    edge_info() { m_left_=Mark(); o_=Object_handle(); forw_=false; }
    LEDA_MEMORY(edge_info)
  };

  void assoc_info(SHalfedge_handle e)  const
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<edge_info>::create(info(e)); 
    geninfo<edge_info>::create(info(e->twin()));
    #else
    info(e)=edge_info();
    info(e->twin())=edge_info();
    #endif
  }

  void discard_info(SHalfedge_handle e)  const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<edge_info>::clear(info(e)); 
    geninfo<edge_info>::clear(info(e->twin()));
    #else
    info(e)=boost::any();
    info(e->twin())=boost::any();
    #endif
  }

  edge_info& ginfo(SHalfedge_handle e)  const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<edge_info>::access(info(e)); 
    #else
    return
      *boost::any_cast<edge_info>(&info(e));
    #endif
  }

  Object_handle& support(SHalfedge_handle e) const
  // uedge information we store in the smaller one 
  { if (&*e < &*(e->twin())) return ginfo(e).o_; 
    else                   return ginfo(e->twin()).o_; }

  Mark& incident_mark(SHalfedge_handle e)  const
  // biedge information we store in the edge
  { return ginfo(e).m_left_; }

  const edge_info& ginfo(SHalfedge_const_handle e)  const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<edge_info>::const_access(info(e)); 
    #else
    return 
      *boost::any_cast<edge_info>(&info(e)); 
    #endif
  }
  const Mark& incident_mark(SHalfedge_const_handle e)  const
  { return ginfo(e).m_left_; }

  bool& is_forward(SHalfedge_handle e) const
  // biedge information we store in the edge
  { return ginfo(e).forw_; }

  void assert_equal_marks(SVertex_handle v1, SVertex_handle v2) const 
  { CGAL_assertion(v1->mark()==v2->mark()); CGAL_USE(v1); CGAL_USE(v2); }

  void assert_equal_marks(SHalfedge_handle e1, SHalfedge_handle e2) const
  { CGAL_assertion(e1->mark()==e2->mark()); CGAL_USE(e1); CGAL_USE(e2); }

  Sphere_segment segment(const Explorer* , 
                         SHalfedge_const_handle e) const
  { return Sphere_segment(
	    e->source()->point(),e->twin()->source()->point(),e->circle()); }

  Sphere_segment trivial_segment(const Explorer* , 
                                 SVertex_const_handle v) const
  { Sphere_point p = v->point(); 
    return Sphere_segment(p,p); }

  Seg_pair two_segments(const Explorer* , 
                        SHalfedge_const_handle e) const
  // we know that source(e)==target(e)
  { return e->circle().split_at(e->source()->point()); }

  Seg_pair two_segments(const Explorer* , 
                        SHalfloop_const_handle l) const
  { return l->circle().split_at_xy_plane(); }

  /*{\Mcreation 6}*/
  SM_triangulator(Sphere_map* M, const Explorer* E,
		  const Sphere_kernel& G = Sphere_kernel()) : Base(M), E_(E), K(G) {}
  /*{\Mcreate |\Mvar| is a triangulator object for the map |M|,
     stores the triangulation in |MT|.}*/

  /*{\Moperations 1.1 1}*/

  void triangulate();
  /*{\Mop produces a triangulated sphere map.}*/

  void triangulate_per_hemisphere(SVertex_iterator start, SVertex_iterator end);

  template <typename Iterator, typename T>
  void partition_to_halfsphere(Iterator start, Iterator end,
    Seg_list& L, CGAL::Unique_hash_map<Iterator,T>& M, int pos) const;

  void merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2);
  void merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2);
  void complete_support(SVertex_iterator v_start, SVertex_iterator v_end, 
			Mark mohs) const;

  void correct_triangle_at(SVertex_handle v)
  { CGAL_NEF_TRACEN("correct_triangle_at "<<PH(v));
    if ( !has_outdeg_two(v) ) return;
    SHalfedge_handle e = first_out_edge(v);
    CGAL_assertion(e->snext()->snext()->snext()==e);
    flip_diagonal(e->snext());
  }

  void dump(std::ostream& os = std::cerr) const
  { SM_io_parser<Explorer>::dump(E_,os);
    SM_io_parser<Base>::dump(*this,os); }

  Sphere_triangle incident_triangle(SHalfedge_handle e) const
  { SHalfedge_handle en(e->snext()), enn(en->snext());
    CGAL_assertion(enn->snext()==e);
    return Sphere_triangle(e->source()->point(),
			   en->source()->point(),
			   enn->source()->point(),
			   e->circle(),
			   en->circle(),
			   enn->circle());
  }

  Sphere_triangle incident_triangle(SHalfedge_const_handle e) const
  { SHalfedge_const_handle en(e->snext()), enn(en->snext());
    CGAL_assertion(enn->snext()==e);
    return Sphere_triangle(e->source()->point(),
			   en->source()->point(),
			   enn->source()->point(),
			   e->circle(),
			   en->circle(),
			   enn->circle());
  }

  void discard_info()
  {
    SVertex_iterator v;
    SHalfedge_iterator e;
    CGAL_forall_svertices(v,*this) discard_info(v);
    CGAL_forall_shalfedges(e,*this) discard_info(e);
  }


}; // SM_triangulator<Decorator_>


template <typename Decorator_>
void SM_triangulator<Decorator_>::triangulate()
{ CGAL_NEF_TRACEN("triangulate");
  // first create sphere segments from isoverts, edges, loops
  Seg_list L;
  Seg_map From;
  SVertex_const_iterator v;
  CGAL_forall_svertices(v,*E_) {
    if ( !E_->is_isolated(v) ) continue;
    L.push_back(trivial_segment(E_,v));
    From[--L.end()] = make_object(v);
  }
  SHalfedge_const_iterator e;
  CGAL_forall_sedges(e,*E_) {
    if ( e->source() == e->twin()->source() ) {
      Seg_pair p = two_segments(E_,e);
      L.push_back(p.first); L.push_back(p.second);
      From[--L.end()] = From[--(--L.end())] = make_object(e);
    } else {
      L.push_back(segment(E_,e));
      From[--L.end()] = make_object(e);
    }
  }
  if ( E_->has_shalfloop() ) {
    Seg_pair p = two_segments(E_,E_->shalfloop());
    L.push_back(p.first); L.push_back(p.second);
    From[--L.end()] = From[--(--L.end())] = make_object(E_->shalfloop());
  }

  // partition segments from L to positive and negative hemisphere
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From, -1);

  //  typename Seg_list::iterator it;
  //    std::cerr << "L_pos" << std::endl;
  //    CGAL_forall_iterators(it,L_pos) std::cerr << *it << std::endl;
  //    std::cerr << "L_neg" << std::endl;
  //    CGAL_forall_iterators(it,L_neg) std::cerr << *it << std::endl;

  // sweep the hemispheres to create two half sphere maps
  typedef SM_subdivision<Self,Seg_iterator,Object_handle> SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_handle v_sep;
  SHalfedge_handle e_sep;
  SM_output O(*this,From); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    PH_geometry());
  SP.sweep();
  v_sep=--this->svertices_end(); e_sep=--this->shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    NH_geometry());
  SM.sweep();
  ++v_sep; ++e_sep;
  // now two CCs of sphere graph are calculated
  // v_sep = first vertex of CC in negative x-sphere
  // e_sep = first edge of CC in negative x-sphere
   
  SHalfedge_iterator u;
  CGAL_forall_sedges(u,*this) {
    Sphere_segment s(u->source()->point(),u->twin()->source()->point());
    u->circle() = s.sphere_circle(); 
    u->twin()->circle() = s.sphere_circle().opposite();
  }

  Mark lower, upper;
  SM_point_locator PL(E_->sphere_map());
  PL.marks_of_halfspheres(lower,upper);
  complete_support(this->svertices_begin(), v_sep, lower);
  complete_support(v_sep, this->svertices_end(), upper);

  /*
  CGAL_forall_sedges(u,*this) {
    std::cerr << point(source(u)) << "->" << point(target(u)) << std::endl;
  }
  */

  // triangulate per hemisphere
  typedef SM_constrained_triang_traits<Self,PH_geometry>  PCT_traits;
  typedef CGAL::generic_sweep<PCT_traits> Positive_halfsphere_ct_sweep;
  typedef SM_constrained_triang_traits<Self,NH_geometry>  NCT_traits;
  typedef CGAL::generic_sweep<NCT_traits> Negative_halfsphere_ct_sweep;
  typedef std::pair<SVertex_iterator,SVertex_iterator> SVertex_pair;

  SVertex_pair vpp(this->svertices_begin(),v_sep);
  Positive_halfsphere_ct_sweep PCTS(vpp, *this,
    PH_geometry());
  PCTS.sweep();
  SVertex_pair vpn(v_sep,this->svertices_end());
  Negative_halfsphere_ct_sweep NCTS(vpn, *this,
    NH_geometry());
  NCTS.sweep();

  /*
  std::cerr << std::endl;
  CGAL_forall_sedges(u,*this) {
    std::cerr << point(source(u)) << "->" << point(target(u)) << std::endl;
  }
  */

  /* Note the we divide the world along the xy equator and 
     split the equator at y- and y+. We treat the halfcircle
     at x+ as if perturbed slightly up. This makes triangles
     that have y- or y+ as a vertex degenerate. if such triangles
     appear we repair it by flipping the edge opposite to the
     vertex y-(y+).
  */

  correct_triangle_at(this->svertices_begin());
  correct_triangle_at(--SVertex_iterator(v_sep));
  correct_triangle_at(v_sep);
  correct_triangle_at(--this->svertices_end());

  CGAL_forall_sedges(u,*this) {
    Sphere_segment s(u->source()->point(),u->twin()->source()->point());
    u->circle() = s.sphere_circle();
    u->twin()->circle() = s.sphere_circle().opposite();
  }

  // merge the hemisphere maps into one sphere map
  merge_halfsphere_maps(this->svertices_begin(),v_sep);
  this->check_integrity_and_topological_planarity(false);
}


template <typename Decorator_>
template <typename Iterator, typename T>
void SM_triangulator<Decorator_>::
partition_to_halfsphere(Iterator start, Iterator beyond, Seg_list& L, 
  CGAL::Unique_hash_map<Iterator,T>& M, int pos) const
{ CGAL_NEF_TRACEN("partition_to_halfsphere ");
  CGAL_assertion(pos!=0);
  bool add_cross = true;
  Sphere_segment s1,s2;
  Sphere_circle xycircle(0,0,pos);
  while ( start != beyond ) { 
    if(start->source().hz() * pos > 0 || start->target().hz() * pos > 0)
      add_cross = false;
    int i = start->intersection(xycircle,s1,s2);
    if (i>1) { L.push_back(s2); M[--L.end()] = M[start]; }
    if (i>0) { L.push_back(s1); M[--L.end()] = M[start]; }
    ++start;
  }
  // now all segments are split into halfspheres
  // we still have to:
  // - split segments containing our special poles y^-, y^+
  // - split halfcircles
  // - add four equator segments 
  Sphere_point S(0,-1,0),N(0,1,0);
  Sphere_circle yzcircle(1,0,0);
  typename Seg_list::iterator it, itl;

  bool part_in_hemisphere(false);
  CGAL_forall_iterators(it,L) { CGAL_NEF_TRACEN("  "<<*it);
    if ( equal_as_sets(it->sphere_circle(),xycircle) ) {
      CGAL_NEF_TRACEN("  splitting xy seg "<<*it);
      int n1 =  it->intersection(yzcircle,s1,s2);
      if (n1 > 1 && !s2.is_degenerate()) 
      { M[ L.insert(it,s2) ] = M[it]; }
      if (n1 > 0 && !s1.is_degenerate()) 
      { M[ L.insert(it,s1) ] = M[it]; }
      int n2 =  it->intersection(yzcircle.opposite(),s1,s2);
      if (n2 > 1 && !s2.is_degenerate()) 
      { M[ L.insert(it,s2) ] = M[it]; }
      if (n2 > 0 && !s1.is_degenerate()) 
      { M[ L.insert(it,s1) ] = M[it]; }
      itl = it; --it; L.erase(itl); M[itl] = T();
      // at least one item was appended
    } else {
      part_in_hemisphere = true;
    }
  }
  CGAL_forall_iterators(it,L) {
    if ( it->is_halfcircle() ) {
      CGAL_NEF_TRACEN("  splitting halfcircle "<<*it);
      Sphere_segment s1,s2;
      it->split_halfcircle(s1,s2);
      *it = s2; 
      M[ L.insert(it,s1) ] = M[it];
    }
  }
  // append 4 xy-equator segments:
  Sphere_segment sp(S,N,xycircle);
  Sphere_segment sm(S,N,xycircle.opposite());
  Sphere_segment s[4];
  sp.split_halfcircle(s[0],s[1]);
  sm.split_halfcircle(s[2],s[3]);
  L.insert(L.end(),s,s+4);
  /* if no segment is covering the interior of the hemisphere
     we have to add a trivial segment to allow for a correct
     triangulation */
  if ( !part_in_hemisphere || add_cross) {
    Sphere_point p(0,0,pos);
    Sphere_circle c(1,0,0);
    L.push_back(Sphere_segment(p,p,c));
  }
}

template <typename Decorator_>
void SM_triangulator<Decorator_>::
merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2)
{
  SVertex_handle v1 = e1->source(), v2 = e2->twin()->source();
  CGAL_NEF_TRACEN("merge_nodes "<<PH(v1)<<PH(v2));
  CGAL_assertion(v1->point()==v2->point());
  SHalfedge_handle ep1 = e1->sprev(), en2 = e2->snext();
  SHalfedge_around_svertex_circulator eav(out_edges(v2)),ee(eav);
  CGAL_For_all(eav,ee) { set_source(eav,v1); }
  link_as_prev_next_pair(e2,e1);  
  link_as_prev_next_pair(ep1,en2); 
  assert_equal_marks(v1,v2);
  discard_info(v2);
  delete_vertex_only(v2);
}


template <typename Decorator_>
void SM_triangulator<Decorator_>::
merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2)
{ CGAL_NEF_TRACEN("merging halfspheres "<<PH(v1)<<PH(v2));
  CGAL_assertion(v1->point()==v2->point());
  std::list<SHalfedge_pair> L_equator;
  SHalfedge_around_sface_circulator 
    ep(last_out_edge(v1)), en(first_out_edge(v2)->twin());
  do { 
   L_equator.push_back(SHalfedge_pair(ep,en));
   merge_nodes(ep,en); ++ep; --en; 
  } while ( ep->source() != v1 );
  
  typename std::list<SHalfedge_pair>::iterator it;
  CGAL_forall_iterators(it,L_equator) { 
    SHalfedge_handle e1 = it->first, e2 = it->second;
    SHalfedge_handle e1t = e1->twin(), e2t = e2->twin();
    CGAL_NEF_TRACEV(PH(e1));CGAL_NEF_TRACEV(PH(e2));
    SHalfedge_handle e2tp = e2t->sprev();
    SHalfedge_handle e2tn = e2t->snext();
    link_as_prev_next_pair(e2tp,e1);
    link_as_prev_next_pair(e1,e2tn);
    SFace_handle f = e2t->incident_sface();
    if ( is_sm_boundary_object(e2t) )
    { undo_sm_boundary_object(e2t,f); store_sm_boundary_object(e1,f); }
    set_face(e1,f);
    if ( e2 == first_out_edge(e2->source()) )
      set_first_out_edge(e2->source(),e1t);
    discard_info(e2);
    delete_edge_pair_only(e2);
  }
}

template <typename Decorator_>
void SM_triangulator<Decorator_>::
complete_support(SVertex_iterator v_start, SVertex_iterator v_end,
		 Mark mohs) const
{ CGAL_NEF_TRACEN("complete_support");
  Mark m_buffer(mohs); 
  for (SVertex_iterator v = v_start; v != v_end; ++v) { 
    CGAL_NEF_TRACEN(" vertex = "<<PH(v));
    SHalfedge_handle e_below = halfedge_below(v);
    if ( v != v_start ) {
      if ( e_below != SHalfedge_handle() ) {
	m_buffer = incident_mark(e_below); 
      } else { // e_below does not exist
	/* this is only the case for a vertex v on the final equatorial
	   halfcircle; there we take the mark from an inedge edge into v */
	//	CGAL_assertion( point(v).z() == 0 && 
	//	  ( pos > 0 ? (point(v).x() >= 0) : (point(v).x()<=0)) );
	m_buffer = incident_mark(first_out_edge(v)->sprev());
      }
    }
    CGAL_NEF_TRACEN(" face mark below "<<m_buffer);

    Object_handle o = support(v);
    SVertex_const_handle vs;
    SHalfedge_const_handle es;
    SHalfloop_const_handle ls;
    if ( o.empty() ) { v->mark() = m_buffer; }
    else if ( CGAL::assign(vs,o) ) { v->mark() = vs->mark(); }
    else if ( CGAL::assign(es,o) ) {
      if ( es->source()->point() == v->point() ) 
	{ v->mark() = es->source()->mark(); }
      else if ( es->twin()->source()->point() == v->point() ) 
	{ v->mark() = es->twin()->source()->mark(); }
      else { v->mark() = es->mark(); }
    }
    else if ( CGAL::assign(ls,o) ) { v->mark() = ls->mark(); }
    else CGAL_error_msg("damn wrong support.");
    CGAL_NEF_TRACEN(" face mark at "<<v->mark());

    if ( is_isolated(v) ) continue;

    SHalfedge_around_svertex_circulator e(first_out_edge(v)), hend(e);
    CGAL_For_all(e,hend) {
      CGAL_NEF_TRACEN("  edge "<<PH(e));
      if ( !is_forward(e) ) break;
      if ( ! support(e).empty() ) {
        SHalfedge_const_handle ei;
        if ( CGAL::assign(ei,support(e)) ) { 
          if ( ei->circle() != e->circle() ) { ei = ei->twin(); }
          CGAL_assertion( ei->circle() == e->circle() ); 
          CGAL_NEF_TRACEN("  supporting edge "<<PH(ei));
          incident_mark(e->twin()) = ei->twin()->incident_sface()->mark();
          e->mark() = ei->mark();
          incident_mark(e) = m_buffer = ei->incident_sface()->mark(); 
        }
        SHalfloop_const_handle li;
        if ( CGAL::assign(li,support(e)) ) { 
          if ( li->circle() != e->circle() ) { li = li->twin(); }
          CGAL_assertion( li->circle() == e->circle() ); 
          CGAL_NEF_TRACEN("  supporting loop "<<PH(li));
          incident_mark(e->twin()) = li->twin()->incident_sface()->mark();
          e->mark() = li->mark();
          incident_mark(e) = m_buffer = li->incident_sface()->mark();
        }
      } else { CGAL_NEF_TRACEN("  support from face below ");
        incident_mark(e->twin()) = e->mark() = 
        incident_mark(e) = m_buffer;
      }
      CGAL_NEF_TRACEN("  new face mark "<<m_buffer);

    } // CGAL_For_all(e,hend)

    CGAL_NEF_TRACEN(" mark of "<<PH(v));
  }

}




} //namespace CGAL
#undef CGAL_USING
#endif //CGAL_SM_TRIANGULATOR_H
