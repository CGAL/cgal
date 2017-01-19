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

#ifndef CGAL_PM_OVERLAYER_H
#define CGAL_PM_OVERLAYER_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Union_find.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 13
#include <CGAL/Nef_2/debug.h>

#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <boost/type_traits/is_same.hpp>

#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t)
#else
#include <LEDA/system/basic.h>
#endif

namespace CGAL {

template <typename PMD, typename I, typename DA>
struct PMO_from_segs {
  typedef PMD Decorator;
  typedef typename Decorator::Vertex_handle   Vertex_handle;
  typedef typename Decorator::Halfedge_handle Halfedge_handle;
  typedef typename Decorator::Point           Point;
  const Decorator& G;
  DA& D; 
  PMO_from_segs(const Decorator& Gi, DA& Di) : 
    G(Gi),D(Di) {}

  Vertex_handle new_vertex(const Point& p)
  { Vertex_handle v = G.new_vertex(p); 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<Halfedge_handle>::create(G.info(v));
    #else
    G.info(v)=Halfedge_handle();
    #endif
    return v;
  }

  void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) 
  { G.link_as_target_and_append(v,e); }

  Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
  { Halfedge_handle e = 
    G.new_halfedge_pair_at_source(v,Decorator::BEFORE); 
    return e;
  }

  void supporting_segment(Halfedge_handle e, I it) const
  { D.supporting_segment(e,it); }

  void trivial_segment(Vertex_handle v, I it) const
  { D.trivial_segment(v,it); }

  void starting_segment(Vertex_handle v, I it) const
  { D.starting_segment(v,it); }

  void passing_segment(Vertex_handle v, I it) const
  { D.passing_segment(v,it); }

  void ending_segment(Vertex_handle v, I it) const
  { D.ending_segment(v,it); }

  void halfedge_below(Vertex_handle v, Halfedge_handle e) const
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<Halfedge_handle>::access(G.info(v)) = e;
    #else
    *boost::any_cast<Halfedge_handle>(&G.info(v)) = e;
    #endif
  }

  Halfedge_handle halfedge_below(Vertex_handle v) const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<Halfedge_handle>::access(G.info(v)); 
    #else
    return 
      boost::any_cast<Halfedge_handle>(G.info(v)); 
    #endif
  }

  void clear_temporary_vertex_info() const
  { Vertex_handle v;
    for(v = G.vertices_begin(); v!= G.vertices_end(); ++v)
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<Halfedge_handle>::clear(G.info(v));
    #else
      G.info(v)=boost::any();
    #endif
  }


}; // PMO_from_segs


template <typename PMD, typename IT, typename INFO>
struct PMO_from_pm {
  typedef PMD Decorator;
  typedef typename PMD::Const_decorator Const_decorator;
  typedef typename Decorator::Vertex_handle Vertex_handle;
  typedef typename Decorator::Halfedge_handle Halfedge_handle;
  typedef typename Decorator::Vertex_const_handle Vertex_const_handle;
  typedef typename Decorator::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Decorator::Point Point;

  const Decorator& G;
  const Const_decorator* pGI[2];
  CGAL::Unique_hash_map<IT,INFO>& M;
  PMO_from_pm(const Decorator& Gi, 
              const Const_decorator* pG0, 
              const Const_decorator* pG1,
              CGAL::Unique_hash_map<IT,INFO>& Mi) : G(Gi),M(Mi) 
 { pGI[0]=pG0; pGI[1]=pG1; }

 Vertex_handle new_vertex(const Point& p) const
 { Vertex_handle v = G.new_vertex(p);
   G.assoc_info(v);
   return v;
 }

 void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) const
 { G.link_as_target_and_append(v,e); }

 Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v) const
 { Halfedge_handle e = 
   G.new_halfedge_pair_at_source(v,Decorator::BEFORE); 
   G.assoc_info(e);
   return e;
 }

 void halfedge_below(Vertex_handle v, Halfedge_handle e) const
 { G.halfedge_below(v) = e; }

 void supporting_segment(Halfedge_handle e, IT it) const
 { INFO& si = M[it];
   CGAL_assertion( si.e != Halfedge_const_handle() );
   G.supp_halfedge(e,si.i) = si.e;
   G.is_forward(e) = true;
 }


 void trivial_segment(Vertex_handle v, IT it) const
 { INFO& si = M[it];
   CGAL_assertion( si.v != Vertex_const_handle() );
   G.supp_vertex(v,si.i) = si.v; 
 }

 void starting_segment(Vertex_handle v, IT it) const
 { INFO& si = M[it];
   G.supp_vertex(v,si.i) = pGI[si.i]->source(si.e);
 }

 void ending_segment(Vertex_handle v, IT it) const
 { INFO& si = M[it];
   G.supp_vertex(v,si.i) = pGI[si.i]->target(si.e);
 }

 void passing_segment(Vertex_handle v, IT it) const
 { INFO& si = M[it];
   G.supp_halfedge(v,si.i) = si.e; 
 }

 Halfedge_handle halfedge_below(Vertex_handle v) const
 { return G.halfedge_below(v); }


}; // PMO_from_pm

/*{\Moptions print_title=yes }*/ 
/*{\Msubst 
PM_decorator_#PMD
Geometry_#GEO
}*/
/*{\Manpage {PM_overlayer}{PMD,GEO}{Plane Map Overlay}{O}}*/
template <typename PM_decorator_, typename Geometry_>
class PM_overlayer : public PM_decorator_ {
  typedef PM_decorator_ Base;
  typedef PM_overlayer<PM_decorator_,Geometry_>  Self;
  const Geometry_& K; // geometry reference

/*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
decorator object offering plane map overlay calculation. Overlay is
either calculated from two plane maps or from a set of segments.  The
result is stored in a plane map |P| that carries the geometry and the
topology of the overlay.

The two template parameters allow to adapt the overlay calculation
to different scenarios.  The template parameter |PM_decorator_| has to
be a model conforming to our plane map decorator concept
|PMDecorator|.  The concept describes the interface how the
topological information stored in |P| can be extracted.  The geometry
|Geometry_| has to be a model conforming to the concept 
|OverlayerGeometry_2|.

The overlay of a set of segments $S$ is stored in a plane map $P =
(V,E,F)$. Vertices are either the endpoints of segments (trivial
segments are allowed) or the result of a non-degenerate internal
intersection of two segments. Between two vertices there is an edge if
there is a segment that supports the straight line embedding of $e$ and
if there is no vertex in the relative interior of the embedding of $e$.

The faces refer to the maximal connected open point sets of the
planar subdivision implied by the embedding of the vertices and edges.
Faces are bounded by possibly several face cycles\cgalFootnote{For the
definition of plane maps and their concepts see the manual page of
|PMConstDecorator|.} including isolated vertices. The overlay process
in the method |create| creates the objects, the topology of the result
and allows to link the plane map objects to input segments by means of
a data accessor. The method starts from zero- and one-dimensional
geometric objects in $S$ and produces a plane map |P| where each point
of the plane can be assigned to an object (vertex, edge, or face) of
|P|.

The overlay of two plane maps $P_i = (V_i, E_i, F_i)$ has the
additional aspect that we already start from two planar subdivisions.
We use the index $i=0,1$ defining the reference to $P_i$, unindexed
variables refer to the resulting plane map $P$.  The $1$-skeleta of
the two maps subdivide the edges and faces of the complementary
structure into smaller units. This means vertices and edges of $P_i$
can split edges of $P_{1-i}$ and face cycles of $P_i$ subdivide faces
of $P_{1-i}$. The 1-skeleton $P'$ of $P$ is defined by the overlay of
the embedding of the 1-skeleta of $P_0$ and $P_1$ (Take a trivial
segment for each vertex and a segment for each edge and use the
overlay definition of a set of segments above). The faces of $P$ refer
to the maximal connected open point sets of the planar subdivision
implied by the embedding of $P'$. Each object from the output tuple
$(V,E,F)$ has a \emph{supporting} object $u_i$ in each of the two
input structures.  Imagine the two maps to be transparencies, which we
stack. Then each point of the plane is covered by an object from each
of the input structures.  This support relation from the input
structures to the output structure defines an information flow. Each
supporting object $u_i$ of $u$ $(i=0,1)$ carries an attribute
$|mark|(u_i)$. After the subdivision operation this attribute
is associated to the output object $u$ by $|mark|(u,i)$.}*/

/*{\Mgeneralization PM_decorator_}*/

public:
/*{\Mtypes 8}*/
  typedef PM_decorator_                 Decorator;
  /*{\Mtypemember the plane map decorator |PM_decorator_|.}*/
  typedef typename Decorator::Plane_map Plane_map;
  /*{\Mtypemember the plane map type decorated by |PM_decorator_|.}*/
  typedef Geometry_                     Geometry;
  /*{\Mtypemember the geometry kernel |Geometry_|.}*/
  typedef typename Geometry::Point_2    Point;
  /*{\Mtypemember the point type of the geometric kernel, 
     \precond |Point| equals |Plane_map::Point|.}*/
  typedef typename Geometry::Segment_2  Segment;
  /*{\Mtypemember the segment type of the geometric kernel.}*/
  typedef typename Decorator::Mark      Mark;
  /*{\Mtypemember the attribute type of plane map objects.}*/


  typedef typename Decorator::Base Const_decorator;
  typedef typename Decorator::Halfedge_handle Halfedge_handle;
  typedef typename Decorator::Vertex_handle Vertex_handle;
  typedef typename Decorator::Face_handle Face_handle;
  typedef typename Decorator::Vertex_iterator Vertex_iterator;
  typedef typename Decorator::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator::Face_iterator Face_iterator;
  typedef typename Decorator::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Decorator::Vertex_const_handle Vertex_const_handle;
  typedef typename Decorator::Face_const_handle Face_const_handle;
  typedef typename Decorator::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Decorator::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Decorator::Face_const_iterator Face_const_iterator;
  typedef typename Decorator::Halfedge_around_vertex_circulator 
    Halfedge_around_vertex_circulator;
  typedef typename Decorator::Halfedge_around_face_circulator 
    Halfedge_around_face_circulator;
  typedef typename Decorator::Hole_iterator Hole_iterator;
  typedef typename Decorator::Isolated_vertex_iterator Isolated_vertex_iterator;

  using Base::clear;
  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::halfedges_begin;
  using Base::halfedges_end;
  using Base::faces_begin;
  using Base::faces_end;
  using Base::number_of_vertices;
  using Base::number_of_halfedges;
  using Base::number_of_faces;
  using Base::new_vertex;
  using Base::new_face;
  using Base::target;
  using Base::source;
  using Base::point;
  using Base::next;
  using Base::previous;
  using Base::twin;
  using Base::info;
  using Base::link_as_outer_face_cycle;
  using Base::link_as_isolated_vertex;
  using Base::link_as_hole;
  using Base::face;
  using Base::set_face;
  using Base::is_isolated;
  using Base::first_out_edge;
  using Base::halfedge;
  using Base::clear_face_cycle_entries;
  using Base::is_closed_at_source;
  using Base::delete_halfedge_pair;
  using Base::delete_face;
  using Base::set_halfedge;
  using Base::set_hole;
  using Base::delete_vertex_only;
  using Base::set_isolated_vertex;
  using Base::has_outdeg_two;
  using Base::merge_halfedge_pairs_at_target;

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

    enum Creation {POLYGON=0, POLYLINE=1};

  /*{\Moperations 1.1 1}*/

  struct Seg_info { // to transport information from input to output
    Halfedge_const_handle e;
    Vertex_const_handle   v;
    int                   i;

    Seg_info() : i(-1) {}
    Seg_info(Halfedge_const_handle e_, int i_) 
    { e=e_; i=i_; }
    Seg_info(Vertex_const_handle v_, int i_) 
    { v=v_; i=i_; }
    Seg_info(const Seg_info& si) 
    { e=si.e; v=si.v; i=si.i; }
    Seg_info& operator=(const Seg_info& si) 
    { e=si.e; v=si.v; i=si.i; return *this; }
    LEDA_MEMORY(Seg_info)
  };

  typedef std::list<Segment>                   Seg_list;
  typedef typename Seg_list::const_iterator    Seg_iterator;
  typedef std::pair<Seg_iterator,Seg_iterator> Seg_it_pair;


/*{\Mcreation 6}*/
PM_overlayer(Plane_map& P, const Geometry& g = Geometry()) : 
/*{\Mcreate |\Mvar| is a decorator object manipulating |P|.}*/
  Base(P), K(g) {}


template <typename Forward_iterator, typename Object_data_accessor>
void create(Forward_iterator start, Forward_iterator end, 
            Object_data_accessor& A, Creation cr = POLYGON) const
/*{\Mop produces in |P| the plane map consistent with the overlay
of the segments from the iterator range |[start,end)|. The data accessor 
|A| allows to initialize created vertices and edges with respect to the
segments in the iterator range. |A| requires the following methods:\\
[[void supporting_segment(Halfedge_handle e, Forward_iterator it)]]\\
[[void trivial_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void starting_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void passing_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void ending_segment(Vertex_handle v, Forward_iterator it)]]\\
where |supporting_segment| is called for each non-trivial segment |*it|
supporting a newly created edge |e|, |trivial_segment| is called for
each trivial segment |*it| supporting a newly created vertex |v|, and
the three last operations are called for each non-trivial segment
|*it| starting at/passing through/ending at the embedding of a newly
created vertex |v|. 
\precond |Forward_iterator| has value type |Segment|.}*/
{
  CGAL_NEF_TRACEN("creating from iterator range");
  CGAL_assertion(cr == POLYGON || cr == POLYLINE);
  typedef PMO_from_segs<Self,Forward_iterator,Object_data_accessor> 
    Output_from_segments;
  typedef Segment_overlay_traits<
    Forward_iterator, Output_from_segments, Geometry> seg_overlay;
  typedef generic_sweep< seg_overlay > seg_overlay_sweep;
  typedef typename seg_overlay::INPUT input_range;
  Output_from_segments Out(*this, A);
  seg_overlay_sweep SOS( input_range(start, end), Out, K);
  SOS.sweep();
  if(cr==POLYGON)
    create_face_objects(Out);
  else
    create_face_objects_pl(Out);
  Out.clear_temporary_vertex_info();
}

void subdivide(const Plane_map& P0, const Plane_map& P1) const
/*{\Mop constructs the overlay of the plane maps |P0| and |P1| in
|P|, where all objects (vertices, halfedges, faces) of |P| are
\emph{enriched} by the marks of the supporting objects of the two
input structures: e.g. let |v| be a vertex supported by a node |v0| in
|P0| and by a face |f1| in |P1| and |D0|, |D1| be decorators of
type |PM_decorator| on |P0|,|P1|. Then |\Mvar.mark(v,0) = D0.mark(v0)|
and |\Mvar.mark(v,1) = D1.mark(f1)|.}*/
{
  Const_decorator PI[2];
  PI[0] = Const_decorator(P0); PI[1] = Const_decorator(P1);
  Seg_list Segments; int i;
  CGAL::Unique_hash_map<Seg_iterator,Seg_info> From;
  for (i=0; i<2; ++i) {
    Vertex_const_iterator v;
    for(v = PI[i].vertices_begin(); v != PI[i].vertices_end(); ++v)
      if ( PI[i].is_isolated(v) ) {
        Segments.push_back(segment(PI[i],v));
        From[--Segments.end()] = Seg_info(v,i);
      }
    Halfedge_const_iterator e;
    for(e = PI[i].halfedges_begin(); e != PI[i].halfedges_end(); ++e)
      if ( is_forward_edge(PI[i],e) ) {
        Segments.push_back(segment(PI[i],e));
        From[--Segments.end()] = Seg_info(e,i);
      }
  }


  typedef PMO_from_pm<Self,Seg_iterator,Seg_info> Output_from_plane_maps;
  typedef Segment_overlay_traits<
    Seg_iterator, Output_from_plane_maps, Geometry> pm_overlay;
  typedef generic_sweep< pm_overlay > pm_overlay_sweep;
  Output_from_plane_maps Out(*this,&PI[0],&PI[1],From);
  pm_overlay_sweep SOS(Seg_it_pair(Segments.begin(),Segments.end()),Out,K);
  SOS.sweep();
  create_face_objects(Out);


  CGAL_NEF_TRACEN("transfering marks");
  Face_iterator f = this->faces_begin(); assoc_info(f);
  for (i=0; i<2; ++i) mark(f,i) = PI[i].mark(PI[i].faces_begin());

  Vertex_iterator v, vend = this->vertices_end();
  for (v = this->vertices_begin(); v != vend; ++v) {
    CGAL_NEF_TRACEN("mark at "<<PV(v));
    Halfedge_handle e_below = halfedge_below(v);
    Mark m_below[2];
    if ( e_below != Halfedge_handle() ) {
      for (int i=0; i<2; ++i) {
        m_below[i] = incident_mark(e_below,i); 
      }
    } else { // e_below does not exist
      for (int i=0; i<2; ++i) 
        m_below[i] = PI[i].mark(PI[i].faces_begin());
    }

    for (i=0; i<2; ++i) 
      if ( supp_halfedge(v,i) != Halfedge_const_handle() ) {
        mark(v,i) = PI[i].mark(supp_halfedge(v,i));
      } else if ( supp_vertex(v,i) != Vertex_const_handle() ) {
        mark(v,i) = PI[i].mark(supp_vertex(v,i));
      } else {
        mark(v,i) = m_below[i];
      }

    if ( is_isolated(v) ) continue;
    Halfedge_around_vertex_circulator 
      e(first_out_edge(v)), hend(e);
    CGAL_For_all(e,hend) {
      if ( is_forward(e) ) {
        CGAL_NEF_TRACEN("   halfedge "<<PE(e));
        Halfedge_const_handle ei;
        bool supported;
        for (int i=0; i<2; ++i) {
          supported = ( supp_halfedge(e,i) != Halfedge_const_handle() );
          if ( supported ) {
            ei = supp_halfedge(e,i);
            CGAL_NEF_TRACEN("   supp halfedge "<<i<<" "<<PE(ei));
            incident_mark(twin(e),i) = 
              PI[i].mark(PI[i].face(PI[i].twin(ei)));
            mark(e,i) = PI[i].mark(ei);
            incident_mark(e,i) = m_below[i] =
              PI[i].mark(PI[i].face(ei));
          } else { // no support from input PI[i]
            incident_mark(twin(e),i) = mark(e,i) = incident_mark(e,i) = 
              m_below[i];
          }
        }
      } else break;
    }

  }
  for (f = ++this->faces_begin(); f != this->faces_end(); ++f) { // skip first face
    assoc_info(f);
    for (i=0; i<2; ++i) mark(f,i) = incident_mark(halfedge(f),i);
  }


}



template <typename Selection>
void select(Selection& predicate) const
/*{\Mop sets the marks of all objects according to the selection
predicate |predicate|. |Selection| has to be a function object type
with a function operator\\ [[Mark operator()(Mark m0, Mark m1)]]\\ For
each object |u| of |P| enriched by the marks of the supporting objects
according to the previous procedure |subdivide|, after this operation
|\Mvar.mark(u) = predicate ( \Mvar.mark(u,0),\Mvar.mark(u,1) )|. The
additional marks are invalidated afterwards. }*/
{ 
  Vertex_iterator vit = this->vertices_begin(),
                  vend = this->vertices_end();
  for( ; vit != vend; ++vit) {
    mark(vit) = predicate(mark(vit,0),mark(vit,1));
    discard_info(vit); 
  }
  Halfedge_iterator hit = this->halfedges_begin(),
                    hend = this->halfedges_end();
  for(; hit != hend; ++(++hit)) {
    mark(hit) = predicate(mark(hit,0),mark(hit,1));
    discard_info(hit);
  }
  Face_iterator fit = this->faces_begin(),
                fend = this->faces_end();
  for(; fit != fend; ++fit) {
    mark(fit) = predicate(mark(fit,0),mark(fit,1));
    discard_info(fit);
  }
}


template <typename Keep_edge>
void simplify(const Keep_edge& keep) const
/*{\Mop simplifies the structure of |P| according to the marks of
its objects. An edge |e| separating two faces |f1| and |f2| and equal
marks |mark(e) == mark(f1) == mark(f2)| is removed and the faces are
unified.  An isolated vertex |v| in a face |f| with |mark(v)==mark(f)|
is removed.  A vertex |v| with outdegree two, two collinear out-edges
|e1|,|e2| and equal marks |mark(v) == mark(e1) == mark(e2)| is removed
and the edges are unified. The data accessor |keep| requires the function
call operator\\[[bool operator()(Halfedge_handle e)]]\\that allows to
avoid the simplification for edge pairs referenced by |e|.}*/
{
  CGAL_NEF_TRACEN("simplifying"); 
  typedef typename CGAL::Union_find<Face_handle>::handle Union_find_handle;
  CGAL::Unique_hash_map< Face_iterator, Union_find_handle> Pitem;
  CGAL::Union_find<Face_handle> unify_faces;

  Face_iterator f, fend = this->faces_end();
  for (f = this->faces_begin(); f!= fend; ++f) { 
     Pitem[f] = unify_faces.make_set(f);
     clear_face_cycle_entries(f);
  }


  Halfedge_iterator e = this->halfedges_begin(), en,
                    eend = this->halfedges_end();
  for(; en=e, ++(++en), e != eend; e=en) { 
    if ( keep(e) ) continue;
    if ( mark(e) == mark(face(e)) &&
         mark(e) == mark(face(twin(e))) ) {
        CGAL_NEF_TRACEN("deleting "<<PE(e));
      if ( !unify_faces.same_set(Pitem[face(e)],
                                 Pitem[face(twin(e))]) ) {
        unify_faces.unify_sets( Pitem[face(e)],
                                Pitem[face(twin(e))] );
        CGAL_NEF_TRACEN("unioning disjoint faces");
      }
      if ( is_closed_at_source(e) )       set_face(source(e),face(e));
      if ( is_closed_at_source(twin(e)) ) set_face(target(e),face(e));
      delete_halfedge_pair(e);
    }
  }
  
  CGAL::Unique_hash_map<Halfedge_handle,bool> linked(false);
  for (e = this->halfedges_begin(); e != eend; ++e) {
    if ( linked[e] ) continue;
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Halfedge_handle e_min = e;
    Face_handle f = *(unify_faces.find(Pitem[face(e)]));
    CGAL_For_all(hfc,hend) {
      set_face(hfc,f);
      if(target(hfc) == target(e_min)) {
	Point p1 = point(source(hfc)), 
          p2 = point(target(hfc)), 
          p3 = point(target(next(hfc)));
	if (!K.left_turn(p1,p2,p3) )
	  e_min = hfc;
      } else if ( K.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      linked[hfc]=true;
    }
    Point p1 = point(source(e_min)),
          p2 = point(target(e_min)),
          p3 = point(target(next(e_min)));
    if ( K.orientation(p1,p2,p3) > 0 ) set_halfedge(f,e_min); // outer
    else set_hole(f,e_min); // store as inner
  }


  Vertex_iterator v, vn, vend = this->vertices_end();
  for(v = this->vertices_begin(); v != vend; v=vn) { CGAL_NEF_TRACEN("at vertex "<<PV(v));
    vn=v; ++vn;
    if ( is_isolated(v) ) {
      if ( mark(v) == mark(face(v)) ) delete_vertex_only(v);
      else set_isolated_vertex(face(v),v); 
    } else { // v not isolated
      Halfedge_handle e2 = first_out_edge(v), e1 = previous(e2);
      Point p1 = point(source(e1)), p2 = point(v), 
            p3 = point(target(e2));
      if ( has_outdeg_two(v) &&
           mark(v) == mark(e1) && mark(v) == mark(e2) &&
           (K.orientation(p1,p2,p3) == 0) ) 
        merge_halfedge_pairs_at_target(e1); 
    }
  }

  Face_iterator fn;
  for (f = this->faces_begin(); f != fend; f=fn) {
    fn=f; ++fn;
    Union_find_handle pit = Pitem[f];
    if ( unify_faces.find(pit) != pit ) delete_face(f);
  }


}

struct vertex_info {
  Mark                  m[2];
  Vertex_const_handle   v_supp[2];
  Halfedge_const_handle e_supp[2];
  Halfedge_handle       e_below;
  vertex_info() 
  { v_supp[0]=v_supp[1]=Vertex_const_handle(); 
    e_supp[0]=e_supp[1]=Halfedge_const_handle(); }
  LEDA_MEMORY(vertex_info)
};

void assoc_info(Vertex_handle v) const
{
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<vertex_info>::create(info(v));
  #else
  info(v)=vertex_info();
  #endif
}

void discard_info(Vertex_handle v) const
{
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<vertex_info>::clear(info(v)); 
  #else
  info(v)=boost::any();
  #endif
}

vertex_info& ginfo(Vertex_handle v) const
{
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  return geninfo<vertex_info>::access(info(v)); 
  #else
  return 
    *boost::any_cast<vertex_info>(&info(v)); 
  #endif
}

Mark& mark(Vertex_handle v, int i) const
{ return ginfo(v).m[i]; }

Vertex_const_handle& supp_vertex(Vertex_handle v, int i) const
{ return ginfo(v).v_supp[i]; }

Halfedge_const_handle& supp_halfedge(Vertex_handle v, int i) const
{ return ginfo(v).e_supp[i]; }

Halfedge_handle& halfedge_below(Vertex_handle v) const
{ return ginfo(v).e_below; }

struct halfedge_info {
  Mark                  m[2];
  Mark                  mf[2];
  Halfedge_const_handle e_supp[2];
  bool                  forw;
  halfedge_info()
  { m[0]=m[1]=mf[0]=mf[1]=Mark(); 
    e_supp[0]=e_supp[1]=Halfedge_const_handle(); 
    forw=false; }
  LEDA_MEMORY(halfedge_info)
};

void assoc_info(Halfedge_handle e)  const
{ 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<halfedge_info>::create(info(e)); 
  geninfo<halfedge_info>::create(info(twin(e)));
  #else
  info(e)=halfedge_info();
  info(twin(e))=halfedge_info();
  #endif
}

void discard_info(Halfedge_handle e)  const
{ 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<halfedge_info>::clear(info(e)); 
  geninfo<halfedge_info>::clear(info(twin(e)));
  #else
  info(e)=boost::any();
  info(twin(e))=boost::any();
  #endif
}

halfedge_info& ginfo(Halfedge_handle e)  const
{
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  return geninfo<halfedge_info>::access(info(e));
  #else
  return 
    *boost::any_cast<halfedge_info>(&info(e));
  #endif
}

Mark& mark(Halfedge_handle e, int i)  const
// uedge information we store in the smaller one 
{ if (&*e < &*(twin(e))) return ginfo(e).m[i]; 
  else                   return ginfo(twin(e)).m[i]; }

Halfedge_const_handle& supp_halfedge(Halfedge_handle e, int i) const
// uedge information we store in the smaller one 
{ if (&*e < &*(twin(e))) return ginfo(e).e_supp[i]; 
  else                   return ginfo(twin(e)).e_supp[i]; }

Mark& incident_mark(Halfedge_handle e, int i)  const
// biedge information we store in the halfedge
{ return ginfo(e).mf[i]; }

bool& is_forward(Halfedge_handle e) const
// biedge information we store in the halfedge
{ return ginfo(e).forw; }

struct face_info {
  Mark m[2];
  face_info() { m[0]=m[1]=Mark(); }
  LEDA_MEMORY(face_info)
};

void assoc_info(Face_handle f)  const
{ 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<face_info>::create(info(f)); 
  #else
  info(f)=face_info();
  #endif
}

void discard_info(Face_handle f)  const
{ 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  geninfo<face_info>::clear(info(f)); 
  #else
  info(f)=boost::any();
  #endif
}

face_info& ginfo(Face_handle f)  const
{ 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  return geninfo<face_info>::access(info(f)); 
  #else
  return 
    *boost::any_cast<face_info>(&info(f)); 
  #endif
}

Mark& mark(Face_handle f, int i)  const
{ return ginfo(f).m[i]; }

void clear_associated_info_of_all_objects() const 
{
  Vertex_iterator vit;
  for (vit = this->vertices_begin(); vit != this->vertices_end(); ++vit)
    discard_info(vit);
  Halfedge_iterator hit;
  for (hit = this->halfedges_begin(); hit != this->halfedges_end(); ++hit) 
    discard_info(hit);
  Face_iterator fit;
  for (fit = this->faces_begin(); fit != this->faces_end(); ++fit) 
    discard_info(fit);
}

template <typename Below_info>
void create_face_objects(const Below_info& D) const
{
  CGAL_NEF_TRACEN("create_face_objects()");
  CGAL::Unique_hash_map<Halfedge_handle,int> FaceCycle(-1);
  std::vector<Halfedge_handle>  MinimalHalfedge;
  int i=0;
  Halfedge_iterator e, eend = this->halfedges_end();
  for (e=this->halfedges_begin(); e != eend; ++e) {
    if ( FaceCycle[e] >= 0 ) continue; // already assigned
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Halfedge_handle e_min = e;
    CGAL_NEF_TRACE("face cycle "<<i<<"\n");
    CGAL_For_all(hfc,hend) {
      FaceCycle[hfc]=i; // assign face cycle number
      if(target(hfc) == target(e_min)) {
	Point p1 = point(source(hfc)), 
          p2 = point(target(hfc)), 
          p3 = point(target(next(hfc)));
	if (!K.left_turn(p1,p2,p3) )
	  e_min = hfc;
      } else if ( K.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      CGAL_NEF_TRACE(PE(hfc));
    } 
    CGAL_NEF_TRACEN("");
    MinimalHalfedge.push_back(e_min); ++i;
  }

  Face_handle f_outer = this->new_face();
  for (int j=0; j<i; ++j) {
    Halfedge_handle e = MinimalHalfedge[j];
      CGAL_NEF_TRACEN("  face cycle "<<j);CGAL_NEF_TRACEN("  minimal halfedge "<<PE(e));
    Point p1 = point(source(e)), 
          p2 = point(target(e)), 
          p3 = point(target(next(e)));
    if ( K.left_turn(p1,p2,p3) ) { // left_turn => outer face cycle
      CGAL_NEF_TRACEN("  creating new face object");
      Face_handle f = this->new_face();
      link_as_outer_face_cycle(f,e);
    }
  }

  for (e = this->halfedges_begin(); e != eend; ++e) {
    if ( face(e) != Face_handle() ) continue;
    CGAL_NEF_TRACEN("linking hole "<<PE(e));
    Face_handle f = determine_face(e,MinimalHalfedge,FaceCycle,D);
    link_as_hole(f,e);
  }
  Vertex_iterator v, v_end = this->vertices_end();
  for (v = this->vertices_begin(); v != v_end; ++v) {
    if ( !is_isolated(v) ) continue;
    Halfedge_handle e_below = D.halfedge_below(v);
    if ( e_below == Halfedge_handle() ) 
      link_as_isolated_vertex(f_outer,v);
    else
      link_as_isolated_vertex(face(e_below),v);    
  }

}

template <typename Below_info>
void create_face_objects_pl(const Below_info& D) const
{
  CGAL_NEF_TRACEN("create_face_objects_pl()");
  CGAL::Unique_hash_map<Halfedge_handle,int> FaceCycle(-1);
  std::vector<Halfedge_handle>  MinimalHalfedge;
  int i=0;
  Halfedge_iterator e, eend = this->halfedges_end();
  for (e=this->halfedges_begin(); e != eend; ++e) {
    if ( FaceCycle[e] >= 0 ) continue; // already assigned
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Halfedge_handle e_min = e;
    CGAL_NEF_TRACE("face cycle "<<i<<"\n");
    CGAL_For_all(hfc,hend) {
      FaceCycle[hfc]=i; // assign face cycle number
      if(target(hfc) == target(e_min)) {
	Point p1 = point(source(hfc)), 
          p2 = point(target(hfc)), 
          p3 = point(target(next(hfc)));
	if (!K.left_turn(p1,p2,p3) )
	  e_min = hfc;
      } else if ( K.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      CGAL_NEF_TRACE(PE(hfc));
    } 
    CGAL_NEF_TRACEN("");
    MinimalHalfedge.push_back(e_min); ++i;
  }

  (void)/* Face_handle f_outer = */ this->new_face();
  for (int j=0; j<i; ++j) {
    Halfedge_handle e = MinimalHalfedge[j];
      CGAL_NEF_TRACEN("  face cycle "<<j);CGAL_NEF_TRACEN("  minimal halfedge "<<PE(e));
    Point p1 = point(source(e)), 
          p2 = point(target(e)), 
          p3 = point(target(next(e)));
    if ( K.left_turn(p1,p2,p3) ) { // left_turn => outer face cycle
      CGAL_NEF_TRACEN("  creating new face object");
      Face_handle f = this->new_face();
      link_as_outer_face_cycle(f,e);
    }
  }

  for (e = this->halfedges_begin(); e != eend; ++e) {
    if ( face(e) != Face_handle() ) continue;
    CGAL_NEF_TRACEN("linking hole "<<PE(e));
    Face_handle f = determine_face(e,MinimalHalfedge,FaceCycle,D);
    link_as_hole(f,e);
  }
}

template <typename Below_info>
Face_handle determine_face(Halfedge_handle e, 
  const std::vector<Halfedge_handle>& MinimalHalfedge,
  const CGAL::Unique_hash_map<Halfedge_handle,int>& FaceCycle,
  const Below_info& D) const
{ CGAL_NEF_TRACEN("determine_face "<<PE(e));
  Halfedge_handle e_min = MinimalHalfedge[FaceCycle[e]];
  Halfedge_handle e_below = D.halfedge_below(target(e_min));
  if ( e_below == Halfedge_handle() ) // below is nirwana
    return this->faces_begin();
  Face_handle f = face(e_below);
  if (f != Face_handle()) return f; // has face already
  f = determine_face(e_below, MinimalHalfedge, FaceCycle,D);
  link_as_hole(f,e_below);
  return f;
}

Segment segment(const Const_decorator& N, 
                Halfedge_const_handle e) const
{ return K.construct_segment(
    N.point(N.source(e)),N.point(N.target(e))); }

Segment segment(const Const_decorator& N, 
                Vertex_const_handle v) const
{ Point p = N.point(v); 
  return K.construct_segment(p,p); }

bool is_forward_edge(const Const_decorator& N, 
                     Halfedge_const_iterator hit) const
{ Point p1 = N.point(N.source(hit));
  Point p2 = N.point(N.target(hit));
  return (K.compare_xy(p1,p2) < 0); }

void assert_type_precondition() const
{ typename PM_decorator_::Point p1; Point p2;
  CGAL_static_assertion((boost::is_same<typename PM_decorator_::Point, Point>::value)); }




}; // PM_overlayer<PM_decorator_,Geometry_>



} //namespace CGAL
#endif // CGAL_PM_OVERLAYER_H
