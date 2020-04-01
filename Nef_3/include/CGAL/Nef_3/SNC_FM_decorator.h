// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_SNC_FM_DECORATOR_H
#define CGAL_SNC_FM_DECORATOR_H

#include <CGAL/license/Nef_3.h>

#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Lazy_kernel.h>
#include <list>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 31
#include <CGAL/Nef_2/debug.h>
#include <CGAL/use.h>

namespace CGAL {

//--------------------------------------------------------------------------
/* We introduce a template that is a point that refers to a vertex.
   We derive from the point type and add the vertex handle */

template <typename P, typename H>
class Vertex_point : public H {
  P p_;
public:
  Vertex_point() {}
  Vertex_point(const P& p, H v) : H(v), p_(p) {}
  Vertex_point(const Vertex_point& vp) : H(vp) { p_=vp.p_; }
  Vertex_point& operator=(const Vertex_point& vp)
  { H::operator=(vp); p_=vp.p_; return *this; }
  H vertex() const { return *this; }
  P point() const { return p_; }
};

/*
bool operator==(const Vertex_point& vp1, const Vertex_point& vp2) {
    return vp1.vertex()==vp2.vertex();
}
*/

template <typename P, typename H>
std::ostream& operator<<(std::ostream& os, const Vertex_point<P,H>& p)
{ os << p.point(); return os; }

template <typename P, typename H>
std::ostream& operator<<(std::ostream& os,
  const std::pair< Vertex_point<P,H>, Vertex_point<P,H> > & s)
{ os << s.first << s.second; return os; }

template <typename V, typename SE>
struct Sort_sedges {
  bool operator()(SE e1, SE e2) {
    V v1[2], v2[2];
    v1[0] = e1->source()->center_vertex();
    v1[1] = e1->next()->source()->center_vertex();
    v2[0] = e2->source()->center_vertex();
    v2[1] = e2->next()->source()->center_vertex();
    int i1(0), i2(0);
    if(CGAL::lexicographically_xyz_smaller(v1[1]->point(),v1[0]->point())) i1 = 1;
    if(CGAL::lexicographically_xyz_smaller(v2[1]->point(),v2[0]->point())) i2 = 1;
    if(v1[i1] != v2[i2]) return CGAL::lexicographically_xyz_smaller(v1[i1]->point(),v2[i2]->point());
    if(v1[1-i1] != v2[1-i2]) return CGAL::lexicographically_xyz_smaller(v1[1-i1]->point(),v2[1-i2]->point());
    return i1<i2;
  }
};

template <typename P, typename SE>
struct Sort_sedges2 {
  bool operator()(SE e1, SE e2) {
    P p1[2], p2[2];
    p1[0] = e1->source()->center_vertex()->point();
    p1[1] = e1->twin()->source()->twin()->center_vertex()->point();
    p2[0] = e2->source()->center_vertex()->point();
    p2[1] = e2->twin()->source()->twin()->center_vertex()->point();
    int i1(0), i2(0);
    if(CGAL::lexicographically_xyz_smaller(p1[1],p1[0])) i1 = 1;
    if(CGAL::lexicographically_xyz_smaller(p2[1],p2[0])) i2 = 1;
    if(p1[i1] != p2[i2]) return CGAL::lexicographically_xyz_smaller(p1[i1],p2[i2]);
    if(p1[1-i1] != p2[1-i2]) return CGAL::lexicographically_xyz_smaller(p1[1-i1],p2[1-i2]);
    return i1<i2;
  }
};


//--------------------------------------------------------------------------
/* The following type is an output model for our generic segment
   sweep module |Segment_overlay_traits|. We use that code to
   sweep all sedges of a plane to finally associate facet cycles
   to facets. Note that we do not use the segment intersection
   functionality of the code as the sedges only touch in endpoints.
   There is room for optimization here. */
//--------------------------------------------------------------------------

template <typename P, typename V, typename E, typename I>
struct Halffacet_output {

Halffacet_output(CGAL::Unique_hash_map<I,E>& F, std::vector<E>& S)
  : From(F), Support(S) { edge_number=0; Support[0]=E(); }

typedef P         Point;
typedef V         Vertex_handle;
typedef unsigned  Halfedge_handle;

CGAL::Unique_hash_map<I,E>& From;
std::vector<E>& Support;
unsigned edge_number;

Vertex_handle new_vertex(const Point& p) const
{
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<unsigned>::create(p.vertex()->info());
  #else
  p.vertex()->info()=unsigned();
  #endif
  return p.vertex(); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
{
  CGAL_USE(v);
  CGAL_NEF_TRACEN("new_edge "<<&*v<<" "<<edge_number+1);
  return ++edge_number;
}

void supporting_segment(Halfedge_handle e, I it)
{ if ( From[it] != E() ) Support[e] = From[it]; }

void halfedge_below(Vertex_handle v, Halfedge_handle e)
{ CGAL_NEF_TRACEN("halfedge_below point "<< v->point() <<": " << e);
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<unsigned>::access(v->info()) = e;
  #else
  v->info() = e;
  #endif
}

// all empty, no update necessary
void link_as_target_and_append(Vertex_handle, Halfedge_handle)
{ /* do nothing */ }

void trivial_segment(Vertex_handle, I) const {}
void starting_segment(Vertex_handle, I) const {}
void passing_segment(Vertex_handle, I) const {}
void ending_segment(Vertex_handle, I) const {}

};

//--------------------------------------------------------------------------
// the following class is the geometry kernel of the generic sweep
// module |Segment_overlay_traits|.
//--------------------------------------------------------------------------

template <typename Point, typename Plane, typename Handle>
class Halffacet_geometry {

 public:
  typedef Point Point_3;
  typedef Plane Plane_3;
  typedef Vertex_point<Point,Handle>  Point_2;
  typedef std::pair<Point_2,Point_2>  Segment_2;

 private:
  // the underlying plane:
  Plane_3 h;

  Point_3 above(const Point_3& p) const
  { return p + h.orthogonal_vector(); }

 public:

Halffacet_geometry(const Plane_3& hi) : h(hi) {}

Point_2 source(const Segment_2& s) const  { return s.first; }
Point_2 target(const Segment_2& s) const  { return s.second; }

bool is_degenerate(const Segment_2& s) const {
  return source(s).vertex()==target(s).vertex();
}

Segment_2 construct_segment(const Point_2& p1, const Point_2& p2) const
{ return Segment_2(p1,p2); }

int orientation(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
{ return static_cast<int>(
    CGAL::orientation(p1.point(),p2.point(),p3.point(),above(p1.point()))); }

int orientation(const Segment_2& s, const Point_2& p) const {
    if(source(s).vertex() == p.vertex() ||
       target(s).vertex() == p.vertex())
        return 0;
    return orientation(source(s),target(s),p);
}

int compare_xy(const Point_2& p1, const Point_2& p2) const {
    if(p1.vertex()==p2.vertex())
        return 0;
    return static_cast<int>(
        CGAL::compare_xyz(p1.point(),p2.point()));
}

Point_2 intersection(const Segment_2& s1, const Segment_2& s2) const
{
  CGAL_USE(s2);
  CGAL_assertion(target(s1).vertex()==target(s2).vertex());
  return target(s1); }

bool left_turn(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
{ return CGAL::orientation(p1,p2,p3,above(p1)) == CGAL::POSITIVE; }

}; // Halffacet_geometry


template<class Kernel, typename SHalfedge_handle, typename Halffacet_geometry>
class SmallerXYZ {

  typedef typename Kernel::Point_3  Point_3;
  Halffacet_geometry& G;
 public:
  SmallerXYZ(Halffacet_geometry& Gin) : G(Gin) {}

  bool operator()(SHalfedge_handle se, SHalfedge_handle min, bool /*init*/) {
    if(se->twin()->source()->twin()->source() == min->twin()->source()->twin()->source()) {
      Point_3 p1 = se->source()->source()->point(),
        p2 = se->twin()->source()->twin()->source()->point(),
        p3 = se->next()->twin()->source()->twin()->source()->point();
      return !G.left_turn(p1,p2,p3);
    }
    return CGAL::lexicographically_xyz_smaller(se->twin()->source()->twin()->source()->point(),
                                               min->twin()->source()->twin()->source()->point());
  }
};

/*
template<typename SHalfedge_handle, typename EK>
class SmallerXYZ<CGAL::Lazy_kernel<EK>, SHalfedge_handle> {

  typedef CGAL::Lazy_kernel<EK> Kernel;
  typedef typename Kernel::Point_3  Point_3;
 public:
  SmallerXYZ() {}

  bool point_in_positive_direction_3(const Point_3& p) const {
    if(p.x() < 0) return false;
    if(p.x() > 0) return true;
    if(p.y() < 0) return false;
    if(p.y() > 0) return true;
    return p.z() > 0;
  }

  bool operator()(const SHalfedge_handle se, const Point_3 min, bool) {
    return
      //      (point_in_positive_direction_3(se->next()->source()->point()) &&
      //            point_in_positive_direction_3(se->next()->twin()->source()->point()) &&
      //            (!init ||
             CGAL::lexicographically_xyz_smaller(se->twin()->source()->twin()->source()->point(),
                                                 min);
  }
};
*/

/*
template<class K2, typename EK, typename SHalfedge_handle>
class SmallerXYZ<CGAL::Lazy_kernel<EK>, K2, SHalfedge_handle> {

  typedef typename K2::Point_3  Point_3;
 public:
  SmallerXYZ() {}

  bool point_in_positive_direction_2(const Point_3& p) const {
    if(p.y() < 0) return false;
    if(p.y() > 0) return true;
    return p.z() > 0;
  }

  bool point_in_positive_direction_3(const Point_3& p) const {
    if(p.x() < 0) return false;
    if(p.x() > 0) return true;
    if(p.y() < 0) return false;
    if(p.y() > 0) return true;
    return p.z() > 0;
  }

  bool operator()(const SHalfedge_handle se, const Point_3 min, bool init) {
    return (point_in_positive_direction_3(se->next()->source()->point()) &&
            point_in_positive_direction_3(se->next()->twin()->source()->point()) &&
            (!init ||
             CGAL::lexicographically_xyz_smaller(se->twin()->source()->twin()->source()->point(),
                                                 min)));
  }
};
*/

//--------------------------------------------------------------------------
// SNC_FM_decorator
// Note that we interpret sedges as edge uses between vertices.  We
//   overwrite some operations from the base class for that semantics.
//--------------------------------------------------------------------------

template <typename S> class SNC_decorator;

template <typename SNC_structure_>
class SNC_FM_decorator : public SNC_decorator<SNC_structure_> {
public:
  typedef SNC_structure_ SNC_structure;
  typedef SNC_decorator<SNC_structure> Base;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Kernel  Kernel;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Mark Mark;

  typedef typename Base::SHalfedge_around_facet_circulator
    SHalfedge_around_facet_circulator;
#if 0
  typedef typename Base::SHalfedge_around_facet_const_circulator
    SHalfedge_around_facet_const_circulator;
#endif

  typedef typename std::list<Object_handle>::iterator
    Object_list_iterator;

  typedef CGAL::Vertex_point<Point_3,Vertex_handle>  Vertex_point;
  typedef std::pair<Vertex_point,Vertex_point> Vertex_segment;
  typedef std::list<Vertex_segment>            Segment_list;
  typedef typename Segment_list::iterator      Segment_iterator;

  typedef CGAL::Halffacet_geometry<Point_3,Plane_3,Vertex_handle>
    Halffacet_geometry;
  using Base::debug; using Base::link_as_facet_cycle; using Base::info; using Base::link_as_interior_loop;
protected:
  Halffacet_handle f_;
public:

  SNC_FM_decorator(SNC_structure& W) : Base(W), f_() {}
  SNC_FM_decorator(SNC_structure& W, Halffacet_handle f)
   : Base(W),f_(f) {}


  Halffacet_cycle_iterator facet_cycles_begin() const
  { return f_->facet_cycles_begin(); }
  Halffacet_cycle_iterator facet_cycles_end()   const
  { return f_->facet_cycles_end(); }

  void create_facet_objects(const Plane_3& h,
    Object_list_iterator start, Object_list_iterator end) const;

protected:

//--------------------------------------------------------------------------
/* We provide some information on determine_facet. To understand its
   functionality please refer to the Nef_2 implementation report where
   there is similar function determine_face with more documentation.

   When we call |determine_facet|$(e,\ldots)$ we know that the
   shalfedge |e| is not linked to a face object yet and thus no
   shalfedge in its face cycle is linked. Thus we jump to the minimal
   shalfedge and look down (along the previously used sweep line in
   the plane of the facet). If we see an shalfedge we ask for its
   face. If it does not have one we recurse.  Note that the target
   vertex of the minimal halfedge actually has a view downwards as we
   examine a hole facet cycle. The method |link_as_facet_cycle| does
   the linkage between the facet object and all sedges of the facet
   cycle. Its cost is linear in the size of the facet cycle. Note also
   that we do the linking bottom up along the recursion stack for all
   visited hole (facet) cycles.  Thus we visit each hole facet cycle
   only once as afterwards each edge of the facet cycle is incident to
   a facet. */
//--------------------------------------------------------------------------

  Halffacet_handle determine_facet(SHalfedge_handle e,
    const std::vector<SHalfedge_handle>& MinimalEdge,
    const CGAL::Unique_hash_map<SHalfedge_handle,int>& FacetCycle,
    const std::vector<SHalfedge_handle>& Edge_of) const
  { CGAL_NEF_TRACEN("  determine_facet "<<debug(e));
    int fc = FacetCycle[e];
    SHalfedge_handle e_min = MinimalEdge[fc];
    SHalfedge_handle e_below =
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      Edge_of[geninfo<unsigned>::access(info(e_min->twin()->source()->twin()->source()))];
    #else
      Edge_of[ boost::any_cast<unsigned>(info(e_min->twin()->source()->twin()->source())) ];
    #endif
    CGAL_assertion( e_below != SHalfedge_handle() );
    CGAL_NEF_TRACEN("  edge below " << debug(e_below));
    Halffacet_handle f = e_below->facet();
    if ( f != Halffacet_handle() ) return f; // has already a facet
    // e_below also has no facet
    f = determine_facet(e_below, MinimalEdge, FacetCycle, Edge_of);
    CGAL_NEF_TRACEN("  edge below " << debug(e_below));
    link_as_facet_cycle(e_below,f);
    link_as_facet_cycle(e_below->twin(),f->twin());
    return f;
  }

//--------------------------------------------------------------------------
/* segment(...)
   an shalfedge |e| can be interpreted as an edge-use extending along
   the halfedge that leaves the local graph in |target(e)|. Only
   edge-uses allow us to code the non-unique character of an edge as a
   boundary object of several facets. An edge-use |e| represents the
   edge |target(e)| in the boundary structure of a facet. */
//--------------------------------------------------------------------------

  Vertex_segment segment(SHalfedge_handle e) const
  { Vertex_handle vs = e->source()->source(),
      vt = e->twin()->source()->twin()->source();
    Vertex_point  ps(vs->point(),vs), pt(vt->point(),vt);
    return Vertex_segment(ps,pt); }

  Vertex_segment segment(SHalfloop_handle l) const
  { Vertex_handle v = l->incident_sface()->center_vertex();
    Vertex_point  p(v->point(),v);
    return Vertex_segment(p,p); }


}; // SNC_FM_decorator<SNC_structure_>



//--------------------------------------------------------------------------
/* create_facet_objects()
   In this method we use the visibility along a sweep line in the
   facet plane to create facet objects. |SEdge_below[v]| either
   provides the sedge that is hit first by a vertical ray downwards or
   an uninitialized edge if there is none. */
//--------------------------------------------------------------------------


template <typename SNC_>
void SNC_FM_decorator<SNC_>::
create_facet_objects(const Plane_3& plane_supporting_facet,
  Object_list_iterator start, Object_list_iterator end) const
{ CGAL_NEF_TRACEN(">>>>>create_facet_objects "
                  << normalized(plane_supporting_facet));
  CGAL::Unique_hash_map<SHalfedge_handle,int> FacetCycle(-1);
  std::vector<SHalfedge_handle> MinimalEdge;
  std::list<SHalfedge_handle> SHalfedges;
  std::list<SHalfloop_handle> SHalfloops;

  CGAL::Unique_hash_map<Vertex_handle,SHalfedge_handle> SHalfedgeBelow;
  CGAL::Unique_hash_map<Segment_iterator,SHalfedge_handle>  From;

  Segment_list Segments;
  SHalfedge_handle e; SHalfloop_handle l;
  typename std::list<SHalfedge_handle>::iterator eit,epred;
  typename std::list<SHalfloop_handle>::iterator lit;

  // the output decorator for the facet plane sweep
  typedef CGAL::Halffacet_output<Vertex_point, Vertex_handle,
    SHalfedge_handle, Segment_iterator> Halffacet_output;

  // the sweep traits class instantiated with the input, output and
  // geometry models
  typedef CGAL::Segment_overlay_traits
    <Segment_iterator, Halffacet_output, Halffacet_geometry>
    Halffacet_sweep_traits;
  typedef CGAL::generic_sweep<Halffacet_sweep_traits>   Halffacet_sweep;
  Halffacet_geometry G(plane_supporting_facet);

  /* We first separate sedges and sloops, and fill a list of segments
     to trigger a sweep. */

  for ( ; start != end; ++start ) {
    if ( CGAL::assign(e,*start) ) {
      SHalfedges.push_back(e);
    } else if ( CGAL::assign(l,*start) ) {
      SHalfloops.push_back(l);
    } else
      CGAL_error_msg("Damn wrong handle.");
  }

  /* We iterate all shalfedges and assign a number for each facet
     cycle.  After that iteration for an edge |e| the number of its
     facet cycle is |FacetCycle[e]| and for a facet cycle |c| we know
     |MinimalEdge[c]|. */
  int i=0;
  //  bool xyplane = plane_supporting_facet.b() != 0 || plane_supporting_facet.c() != 0;
  SmallerXYZ<Kernel, SHalfedge_handle, Halffacet_geometry> smallerXYZ(G);
  CGAL_forall_iterators(eit,SHalfedges) {
    e = *eit;
    if ( FacetCycle[e] >= 0 ) continue; // already assigned
    SHalfedge_around_facet_circulator hfc(e),hend(hfc);
    FacetCycle[hfc]=i;
    SHalfedge_handle e_min = e;
    bool init=false;
    CGAL_NEF_TRACEN("\n  facet cycle numbering (up) "<<i);
    CGAL_For_all(hfc,hend) {
      FacetCycle[hfc]=i; // assign face cycle number
      if(smallerXYZ(hfc, e_min, init)) {
        init = true;
        e_min = hfc;
      }

      CGAL_NEF_TRACEN(hfc->twin()->source()->twin()->source()->point() << " lex xyz smaller " <<
                      e_min->twin()->source()->twin()->source()->point() << "=" <<
                      CGAL::lexicographically_xyz_smaller(hfc->twin()->source()->twin()->source()->point(),
                                                          e_min->twin()->source()->twin()->source()->point()));

    } CGAL_NEF_TRACEN("");
    MinimalEdge.push_back(e_min);
    ++i;
  }



  /* We now know the number of facet cycles |i| and we have a minimal
     edge |e| for each facet cycle. We just check the geometric
     embedding of |e| and |e->next()| to characterize the facet cycle
     (outer or hole). Note that the two edges cannot be collinear due
     to the minimality of |e| (the lexicographic minimality of the
     embedding of its target vertex). Outer facet cycles obtain facet
     objects right away. */

  for (int j=0; j<i; ++j) {
    SHalfedge_handle e = MinimalEdge[j];
    CGAL_NEF_TRACEN("  facet cycle "<<j<<" minimal halfedge "<<debug(e));
    Point_3 p1 = e->source()->source()->point(),
      p2 = e->twin()->source()->twin()->source()->point(),
      p3 = e->next()->twin()->source()->twin()->source()->point();

    //      std::cerr << "minimal shalfedge " << e->source()->source()->point() << ":"
    //                << e->source()->point() << "->" << e->twin()->source()->point() << std::endl;


    if ( G.left_turn(p1,p2,p3) ) {
      Halffacet_handle f = this->sncp()->new_halffacet_pair(plane_supporting_facet);
      link_as_facet_cycle(e,f); link_as_facet_cycle(e->twin(),f->twin());
      f->mark() = f->twin()->mark() = e->mark();
      CGAL_NEF_TRACEN("  creating new facet object "<<&*f<<" bd "<<&*e);
    }
  }

  /* Now the only shalfedges not linked are those on hole facet
     cycles.  We use a recursive scheme to find the bounding cycle
     providing the facet object and finally iterate over all isolated
     vertices to link them accordingly to their containing facet
     object. Note that in that final iteration all shalfedges already
     have facet links. Thus that ensures termination. The recursive
     operation $|determine_facet|(e,\ldots)$ returns the facet
     containing the hole cycle of |e|. As a postcondition of this code
     part we have all shalfedges and isolated shalfloops linked to
     facet objects, and all facet objects know their bounding facet
     cycles. */

  bool do_sweep = false;
  if(SHalfloops.size() > 0)
    do_sweep = true;

  CGAL_forall_iterators(eit,SHalfedges) {
    //    std::cerr << "fc " << FacetCycle[*eit] << std::endl;
    if ( (*eit)->facet() == Halffacet_handle() ) {
      //      std::cerr << "nicht verlinkte shalfedge " << (*eit)->source()->source()->point() << ":"
      //                << (*eit)->source()->point() << "->" << (*eit)->twin()->source()->point() << std::endl;
      do_sweep = true;
      break;
    }
  }


  //  std::cerr << std::endl;
#ifndef CGAL_NEF3_PLANE_SWEEP_OPTIMIZATION_OFF
  if(!do_sweep) return;
#endif

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
  number_of_plane_sweeps++;
  timer_plane_sweeps.start();
#endif

  //  Note that we only allow those edges that are
  //  directed from lexicographically smaller to larger vertices.
  //  Insertion of SHalfedges into Segments is shifted below in order
  //  to guarantee that there are no gaps in the overlay.


  //  SHalfedges.sort(Sort_sedges2<Point_3,SHalfedge_handle>());
  SHalfedges.sort(Sort_sedges<Vertex_handle,SHalfedge_handle>());
  for(eit = SHalfedges.begin();eit != SHalfedges.end();) {
    CGAL_NEF_TRACEN("  appending edge "<< debug(*eit));
    Segments.push_front(segment(*eit));
    From[Segments.begin()] = *eit;
    epred=eit;
    ++eit;
    if(eit != SHalfedges.end()) {
      CGAL_NEF_TRACEN("test " << std::endl << "  " << debug(*epred)
             << std::endl << "  " << debug(*eit));
    }
    if(eit != SHalfedges.end() &&
       (*epred)->source()->source() ==(*eit)->next()->source()->source() &&
       (*eit)->source()->source() == (*epred)->next()->source()->source())
      ++eit;
  }

  CGAL_forall_iterators(lit,SHalfloops) {
    CGAL_NEF_TRACEN("  appending loop " << (*lit)->incident_sface()->center_vertex()->point());
    Segments.push_back(segment(*lit));
  }

  std::vector<SHalfedge_handle> Edge_of(Segments.size()+1);
  Halffacet_output O(From,Edge_of);
  Halffacet_sweep f_s(typename Halffacet_sweep::INPUT(
    Segments.begin(),Segments.end()), O, G); f_s.sweep();

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
  timer_plane_sweeps.stop();
#endif

  CGAL_forall_iterators(eit,SHalfedges) { e=*eit;
  if ( e->facet() != Halffacet_handle() ) continue;
    CGAL_NEF_TRACEN("  linking hole "<<debug(e));
    Halffacet_handle f = determine_facet(e,MinimalEdge,FacetCycle,Edge_of);
    link_as_facet_cycle(e,f); link_as_facet_cycle(e->twin(),f->twin());
  }

  CGAL_forall_iterators(lit,SHalfloops) { l=*lit;
    SHalfedge_handle e_below =
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      Edge_of[geninfo<unsigned>::access(info(l->incident_sface()->center_vertex()))];
    #else
      Edge_of[ boost::any_cast<unsigned>(info(l->incident_sface()->center_vertex())) ];
    #endif

    CGAL_assertion( e_below != SHalfedge_handle() );
    CGAL_NEF_TRACEN("link sloop at vertex "<< l->incident_sface()->center_vertex()->point());
    CGAL_NEF_TRACEN("e_below "  << debug(e_below));
    CGAL_NEF_TRACEN("next    "  << debug(e_below->next()));
    CGAL_NEF_TRACEN("next    "  << debug(e_below->next()->next()));
    CGAL_NEF_TRACEN("next    "  << debug(e_below->next()->next()->next()));
    CGAL_NEF_TRACEN("next    "  << debug(e_below->next()->next()->next()->next()));
    link_as_interior_loop(l,e_below->facet());
    link_as_interior_loop(l->twin(),e_below->facet()->twin());
  }

  CGAL_NEF_TRACEN("exit FM");
}

} //namespace CGAL
#endif //CGAL_SNC_FM_DECORATOR_H
