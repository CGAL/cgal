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
// file          : include/CGAL/Nef_3/SNC_FM_decorator.h
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
// implementation: operations on a facet substructure of an SNC
// =========================================================================== 

#ifndef CGAL_SNC_FM_DECORATOR_H
#define CGAL_SNC_FM_DECORATOR_H

#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <list>
#undef _DEBUG
#define _DEBUG 31
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------------
/* We introduce a template that is a point that refers to a vertex.
   We derive from the point type and add the vertex handle */

template <typename P, typename H>
class Vertex_point : public P {
  H v_;
public:
  Vertex_point() {}
  Vertex_point(const P& p, H v) : P(p)  { v_=v; }
  Vertex_point(const Vertex_point& vp) : P(vp) { v_=vp.v_; }
  Vertex_point& operator=(const Vertex_point& vp)
  { P::operator=(vp); v_=vp.v_; return *this; } 
  H vertex() const { return v_; }
};

template <typename P, typename H>
std::ostream& operator<<(std::ostream& os, const Vertex_point<P,H>& p)
{ os << static_cast<P>(p); return os; }

template <typename P, typename H>
std::ostream& operator<<(std::ostream& os, 
  const std::pair< Vertex_point<P,H>, Vertex_point<P,H> > & s)
{ os << s.first << s.second; return os; }


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
{ geninfo<unsigned>::create(p.vertex()->info());
  return p.vertex(); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v) 
{ TRACEN("new_edge "<<&*v<<" "<<edge_number+1);
  return ++edge_number; }

void supporting_segment(Halfedge_handle e, I it)
{ if ( From[it] != E() ) Support[e] = From[it]; }

void halfedge_below(Vertex_handle v, Halfedge_handle e)
{ TRACEN("halfedge_below "<<&*v<<" "<<e); 
  geninfo<unsigned>::access(v->info()) = e; }

// all empty, no update necessary
void link_as_target_and_append(Vertex_handle v, Halfedge_handle e)
{ /* do nothing */ }

void trivial_segment(Vertex_handle v, I it) const {}
void starting_segment(Vertex_handle v, I it) const {}
void passing_segment(Vertex_handle v, I it) const {}
void ending_segment(Vertex_handle v, I it) const {}

};

//--------------------------------------------------------------------------
// the following class is the geometry kernel of the generic sweep
// module |Segment_overlay_traits|. 
//--------------------------------------------------------------------------

template <typename Point, typename Plane, typename Handle>
class Halffacet_geometry { 

typedef Point Point_3;
typedef Plane Plane_3;
typedef Vertex_point<Point,Handle>  Point_2;
typedef std::pair<Point_2,Point_2>  Segment_2;

  // the underlying plane:
  Plane_3 h;

  Point_3 above(const Point_3& p) const
  { return p + h.orthogonal_vector(); }

public:

Halffacet_geometry(const Plane_3& hi) : h(hi) {}

Point_2 source(const Segment_2& s) const  { return s.first; }
Point_2 target(const Segment_2& s) const  { return s.second; }

bool is_degenerate(const Segment_2& s) const
{ return source(s)==target(s); }

Segment_2 construct_segment(const Point_2& p1, const Point_2& p2) const
{ return Segment_2(p1,p2); }

int orientation(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
{ return static_cast<int>(
    CGAL::orientation(p1,p2,p3,above(p1))); }

int orientation(const Segment_2& s, const Point_2& p) const
{ return orientation(source(s),target(s),p); }

int compare_xy(const Point_2& p1, const Point_2& p2) const
{ return static_cast<int>(
    CGAL::compare_lexicographically_xyz(p1,p2)); }

Point_2 intersection(const Segment_2& s1, const Segment_2& s2) const
{ CGAL_nef3_assertion(target(s1)==target(s2)); 
  return target(s1); }

bool leftturn(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
{ return CGAL::orientation(p1,p2,p3,above(p1)) == CGAL::POSITIVE; }

}; // Halffacet_geometry


//--------------------------------------------------------------------------
// SNC_FM_decorator
/* Note that we interpret sedges as edge uses between vertices.  We
   overwrite some operations from the base class for that semantics. */
//--------------------------------------------------------------------------


template <typename SNC_structure_>
class SNC_FM_decorator : public SNC_decorator<SNC_structure_> {
public:
  typedef SNC_structure_ SNC_structure;
  typedef SNC_decorator<SNC_structure> Base;

#define USING(t) typedef typename SNC_structure_::t t
  USING(Vertex_iterator); USING(Vertex_handle);
  USING(Halfedge_iterator); USING(Halfedge_handle);
  USING(Halffacet_iterator); USING(Halffacet_handle);
  USING(Volume_iterator); USING(Volume_handle);
  USING(SVertex_iterator); USING(SVertex_handle);
  USING(SHalfedge_iterator); USING(SHalfedge_handle);
  USING(SHalfloop_iterator); USING(SHalfloop_handle);
  USING(SFace_iterator); USING(SFace_handle);
  USING(SFace_cycle_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Shell_entry_iterator);
  USING(SObject_handle);
  USING(Point_3);
  USING(Plane_3);
  USING(Mark);
#undef USING

  typedef typename Base::SHalfedge_around_facet_circulator
    SHalfedge_around_facet_circulator;
#if 0
  typedef typename Base::SHalfedge_around_facet_const_circulator
    SHalfedge_around_facet_const_circulator;
#endif

  typedef typename std::list<SObject_handle>::iterator 
    SObject_list_iterator;

  typedef Vertex_point<Point_3,Vertex_handle>  Vertex_point;
  typedef std::pair<Vertex_point,Vertex_point> Vertex_segment;
  typedef std::list<Vertex_segment>            Segment_list;
  typedef typename Segment_list::iterator      Segment_iterator;

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
    SObject_list_iterator start, SObject_list_iterator end) const;


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
  { TRACEN("  determine_facet "<<debug(e));
    int fc = FacetCycle[e];
    SHalfedge_handle e_min = MinimalEdge[fc];
    SHalfedge_handle e_below = 
      Edge_of[geninfo<unsigned>::access(info(target(e_min)))];
    CGAL_nef3_assertion( e_below != SHalfedge_handle() );
    Halffacet_handle f = facet(e_below);
    if ( f != Halffacet_handle() ) return f; // has already a facet 
    // e_below also has no facet
    f = determine_facet(e_below, MinimalEdge, FacetCycle, Edge_of);
    link_as_facet_cycle(e_below,f); 
    link_as_facet_cycle(twin(e_below),twin(f)); 
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
  { Vertex_handle vs = source(e), vt = target(e); 
    Vertex_point  ps(point(vs),vs), pt(point(vt),vt);
    return Vertex_segment(ps,pt); }

  Vertex_segment segment(SHalfloop_handle l) const
  { Vertex_handle v = vertex(l); 
    Vertex_point  p(point(v),v);
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
  SObject_list_iterator start, SObject_list_iterator end) const
{ TRACEN(">>>>>create_facet_objects");

  CGAL::Unique_hash_map<SHalfedge_handle,int> FacetCycle(-1);
  CGAL::Unique_hash_map<Vertex_handle,SHalfedge_handle> SHalfedgeBelow;
  CGAL::Unique_hash_map<Segment_iterator,SHalfedge_handle>  From;
  std::vector<SHalfedge_handle> MinimalEdge;
  std::list<SHalfedge_handle> SHalfedges; 
  std::list<SHalfloop_handle> SHalfloops; 
  Segment_list Segments;
  SHalfedge_handle e; SHalfloop_handle l;
  typename std::list<SHalfedge_handle>::iterator eit;
  typename std::list<SHalfloop_handle>::iterator lit;

  // the output decorator for the facet plane sweep
  typedef CGAL::Halffacet_output<Vertex_point, Vertex_handle, 
    SHalfedge_handle, Segment_iterator> Halffacet_output;

  // the geometry kernel for the facet plane sweep
  typedef CGAL::Halffacet_geometry<Point_3,Plane_3,Vertex_handle> 
    Halffacet_geometry;

  // the sweep traits class instantiated with the input, output and
  // geometry models
  typedef CGAL::Segment_overlay_traits
    <Segment_iterator, Halffacet_output, Halffacet_geometry>
    Halffacet_sweep_traits;
  typedef CGAL::generic_sweep<Halffacet_sweep_traits>   Halffacet_sweep;
  Halffacet_geometry G(plane_supporting_facet);

  /* We first separate sedges and sloops, and fill a list of segments
     to trigger a sweep. Note that we only allow those edges that are
     directed from lexicographically smaller to larger vertices.  */

  for ( ; start != end; ++start ) {
    if ( CGAL::assign(e,*start) ) 
    { SHalfedges.push_back(e); TRACEN("  appending "<<debug(e));
      Segments.push_back(segment(e)); From[--Segments.end()] = e; }
    else if ( CGAL::assign(l,*start) ) 
    { SHalfloops.push_back(l); Segments.push_back(segment(l)); }
    else CGAL_nef3_assertion_msg(0,"Damn wrong handle.");
  }

  std::vector<SHalfedge_handle> Edge_of(Segments.size()+1);
  Halffacet_output O(From,Edge_of);
  Halffacet_sweep FS(Halffacet_sweep::INPUT(
    Segments.begin(),Segments.end()), O, G); FS.sweep();

  /* We iterate all shalfedges and assign a number for each facet
     cycle.  After that iteration for an edge |e| the number of its
     facet cycle is |FacetCycle[e]| and for a facet cycle |c| we know
     |MinimalEdge[c]|. */

  int i=0; 
  CGAL_nef3_forall_iterators(eit,SHalfedges) { e = *eit;
    if ( FacetCycle[e] >= 0 ) continue; // already assigned
    SHalfedge_around_facet_circulator hfc(e),hend(hfc);
    SHalfedge_handle e_min = e;
    TRACEN("  facet cycle numbering "<<i<<"\n    ");
    CGAL_For_all(hfc,hend) {
      FacetCycle[hfc]=i; // assign face cycle number
      if ( CGAL::lexicographically_xyz_smaller(
             point(target(hfc)), point(target(e_min))))
        e_min = hfc;
      TRACE(debug(hfc));
    } TRACEN("");
    MinimalEdge.push_back(e_min);
    ++i;
  }

  /* We now know the number of facet cycles |i| and we have a minimal
     edge |e| for each facet cycle. We just check the geometric
     embedding of |e| and |next(e)| to characterize the facet cycle
     (outer or hole). Note that the two edges cannot be collinear due
     to the minimality of |e| (the lexicographic minimality of the
     embedding of its target vertex). Outer facet cycles obtain facet
     objects right away. */

  for (int j=0; j<i; ++j) {
    SHalfedge_handle e = MinimalEdge[j];
    TRACEN("  facet cycle "<<j<<" minimal halfedge "<<debug(e));
    Point_3 p1 = point(source(e)), 
            p2 = point(target(e)), 
            p3 = point(target(next(e)));
    if ( G.leftturn(p1,p2,p3) ) { 
      Halffacet_handle f = sncp()->new_halffacet_pair(plane_supporting_facet);
      link_as_facet_cycle(e,f); link_as_facet_cycle(twin(e),twin(f)); 
      mark(f) = mark(e); 
      TRACEN("  creating new facet object "<<&*f<<" bd "<<&*e);
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

  CGAL_nef3_forall_iterators(eit,SHalfedges) { e=*eit;
    if ( facet(e) != Halffacet_handle() ) continue;
    TRACEN("  linking hole "<<debug(e));
    Halffacet_handle f = determine_facet(e,MinimalEdge,FacetCycle,Edge_of);
    link_as_facet_cycle(e,f); link_as_facet_cycle(twin(e),twin(f));
  }
  CGAL_nef3_forall_iterators(lit,SHalfloops) { l=*lit;
    SHalfedge_handle e_below = 
      Edge_of[geninfo<unsigned>::access(info(vertex(l)))];
    TRACEN("vertex "<<&*vertex(l));
    CGAL_nef3_assertion( e_below != SHalfedge_handle() );
    link_as_interior_loop(l,facet(e_below));
    link_as_interior_loop(twin(l),twin(facet(e_below)));
  }

}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_FM_DECORATOR_H


