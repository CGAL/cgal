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
// file          : include/CGAL/Nef_3/SNC_SM_overlayer.h
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
// implementation: Overlay module for sphere maps
// ============================================================================

/* IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT 
   This module is mostly indentical to SM_overlayer.h of Nef_S2 Please
   visit the documentation of Nef_S2/noweb/SM_overlayer.lw for an
   elaborate treatment of sphere map overlay. */

#ifndef CGAL_SNC_SM_OVERLAYER_H
#define CGAL_SNC_SM_OVERLAYER_H

#include <CGAL/basic.h>
#include <CGAL/Union_find.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_3/SNC_SM_decorator.h>
#include <CGAL/Nef_3/SNC_SM_io_parser.h>
#undef _DEBUG
#define _DEBUG 131
#include <CGAL/Nef_3/debug.h>


#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

CGAL_BEGIN_NAMESPACE

// ============================================================================
// an output model for |Segment_overlay_traits| from segments
// beware the SHandle names goto Handle names
// ============================================================================

template <typename Decorator_, typename I>
struct SMO_from_segs {
  typedef Decorator_ Decorator;
  typedef typename Decorator::SVertex_handle    Vertex_handle;
  typedef typename Decorator::SHalfedge_handle  Halfedge_handle;
  typedef typename Decorator::Sphere_point      Point;
  typedef typename Decorator::Sphere_segment    Segment;
  typedef CGAL::Unique_hash_map<I,bool> Iterator_map;
  Decorator G;
  const Iterator_map& M;
  SMO_from_segs(Decorator_ Gi, const Iterator_map& Mi) : G(Gi),M(Mi) {}

  Vertex_handle new_vertex(const Point& p)
  { Vertex_handle v = G.new_vertex(p); 
    geninfo<Halfedge_handle>::create(G.info(v));
    return v;
  }

  void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) 
  { G.link_as_target_and_append(v,e); }

  Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
  { Halfedge_handle e = 
    G.new_edge_pair_at_source(v,Decorator::BEFORE); 
    return e;
  }

  void supporting_segment(Halfedge_handle e, I it) const
  { if ( M[it] ) G.mark(e) = true; }

  void trivial_segment(Vertex_handle v, I it) const
  { if ( M[it] ) G.mark(v) = true; }

  void starting_segment(Vertex_handle v, I it) const
  { if ( M[it] ) G.mark(v) = true; }

  void passing_segment(Vertex_handle v, I it) const
  { if ( M[it] ) G.mark(v) = true; }

  void ending_segment(Vertex_handle v, I it) const
  { if ( M[it] ) G.mark(v) = true; }

  void halfedge_below(Vertex_handle v, Halfedge_handle e) const
  { geninfo<Halfedge_handle>::access(G.info(v)) = e; }

  Halfedge_handle halfedge_below(Vertex_handle v) const
  { return geninfo<Halfedge_handle>::access(G.info(v)); }

  void assert_equal_marks(Vertex_handle v1, Vertex_handle v2) const 
  { CGAL_nef3_assertion(G.mark(v1)==G.mark(v2)); }

  void discard_info(Vertex_handle v) const 
  { geninfo<Halfedge_handle>::clear(G.info(v)); }

  void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
  { CGAL_nef3_assertion(G.mark(e1)==G.mark(e2)); }

  void discard_info(Halfedge_handle e) const {}

  void clear_temporary_vertex_info() const
  { Vertex_handle v;
    CGAL_nef3_forall_svertices(v,G)
      geninfo<Halfedge_handle>::clear(G.info(v));
  }


}; // SMO_from_segs



// ============================================================================
// an output model for |Segment_overlay_traits| from two sphere maps
// beware the SHandle names goto Handle names
// ============================================================================

template <typename Decorator_, typename IT, typename INFO>
struct SMO_from_sm {
  typedef Decorator_ Overlayer;
  typedef typename Decorator_::Base Decorator;
  typedef typename Decorator::SVertex_handle     Vertex_handle;
  typedef typename Decorator::SHalfedge_handle   Halfedge_handle;
  typedef typename Decorator::Sphere_point       Point;
  typedef typename Decorator::Sphere_segment     Segment;
  Overlayer G;
  Decorator* pGI;
  CGAL::Unique_hash_map<IT,INFO>& M;
  SMO_from_sm(Overlayer Gi, 
              Decorator* pGIi, 
              CGAL::Unique_hash_map<IT,INFO>& Mi) : 
    G(Gi), pGI(pGIi), M(Mi) {}

Vertex_handle new_vertex(const Point& p) const
{ Vertex_handle v = G.new_vertex(p);
  G.assoc_info(v);
  return v;
}

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) const
{ G.link_as_target_and_append(v,e); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v) const
{ Halfedge_handle e = 
  G.new_edge_pair_at_source(v,Decorator::BEFORE); 
  G.assoc_info(e);
  return e;
}

void halfedge_below(Vertex_handle v, Halfedge_handle e) const
{ G.halfedge_below(v) = e; }

void supporting_segment(Halfedge_handle e, IT it) const
{ INFO& si = M[it];
  G.is_forward(e) = true;
  if ( si._from == -1 )  return; // equatorial segment
  G.supp_object(e,si._from) = si._o;
  TRACEN("   supporting segment "<<si._from<<" "<<*it);
}


void trivial_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  CGAL_nef3_assertion( si._o != NULL );
  G.supp_object(v,si._from) = si._o; 
}

void starting_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  if ( si._from == -1 ) return;
  G.supp_object(v,si._from) = si._o;
}

void ending_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  if ( si._from == -1 ) return;
  G.supp_object(v,si._from) = si._o;
}

void passing_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  if ( si._from == -1 ) return;
  G.supp_object(v,si._from) = si._o; 
}

Halfedge_handle halfedge_below(Vertex_handle v) const
{ return G.halfedge_below(v); }

void assert_equal_marks(Vertex_handle v1, Vertex_handle v2) const 
{ TRACEV(G.mark(v1,0));TRACEV(G.mark(v1,1));
  TRACEV(G.mark(v2,0));TRACEV(G.mark(v2,1));
  CGAL_nef3_assertion(G.mark(v1,0)==G.mark(v2,0)&&
		      G.mark(v1,1)==G.mark(v2,1)); }
void discard_info(Vertex_handle v) const 
{ G.discard_info(v); }

void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
{ CGAL_nef3_assertion(G.mark(e1,0)==G.mark(e2,0) && 
		      G.mark(e1,1)==G.mark(e2,1)); }

void discard_info(Halfedge_handle e) const 
{ G.discard_info(e); }

}; // SMO_from_sm

// ============================================================================
// ============================================================================

/*{\Manpage {SNC_SM_overlayer}{Refs_}{Overlay in the sphere}{O}}*/

template <typename Refs_>
class SNC_SM_overlayer : public SNC_SM_decorator<Refs_> {
public:

  /*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
  decorator object offering sphere map overlay calculation. Overlay is
  either calculated from two sphere maps or from a set of halfspaces.
  The result is stored in a sphere map |M| that carries the geometry and
  the topology of the overlay.

  The template parameter provides the underlying topological interface
  to sphere maps. The template parameter |Decorator_| has to be a model
  conforming to our map decorator concept |SNC_SM_decorator|.  The concept
  also describes the interface how the topological information stored in
  |M| can be extracted or extended.

  The overlay of a set of sphere segments $S$ is stored in a sphere map
  $M = (V,E,L,F)$. Vertices are either the endpoints of segments (trivial
  segments are allowed) or the result of the internal intersection of
  two segments. Between two vertices there's an edge if there's a
  segment that supports the spherical embedding of $e$ and if there's no
  vertex in the relative interior of the embedding of $e$.

  The faces refer to the maximal connected open point sets of the
  spherical subdivision implied by the embedding of the vertices and
  edges.  Faces are bounded by possibly several face cycles\footnote{For
  the definition of sphere maps and their concepts see the manual page
  of |SNC_SM_decorator|.} including isolated vertices. The overlay process
  in the method |create_from_segments| creates the objects and the
  topology of the result. The method starts from zero- and
  one-dimensional geometric objects in $S$ and produces a spherical
  structure where each point of the sphere can be assigned to an object
  (vertex, edge, loop, or face) of |M|.

  The overlay of two sphere maps $M_i = (V_i, E_i, L_i, F_i)$ has the
  additional aspect that we already start from two spherical
  subdivisions.  We use the index $i=0,1$ defining the reference to
  $M_i$, unindexed variables refer to the resulting sphere map $M$.  The
  $1$-skeleta of the two maps subdivide the edges, loops, and faces of
  the complementary structure into smaller units. This means vertices,
  edges, and loops of $M_i$ can split edges and loops of $M_{1-i}$ and
  face cycles of $M_i$ subdivide faces of $M_{1-i}$. The 1-skeleton $G$
  of $M$ is defined by the overlay of the embedding of the 1-skeleta of
  $M_0$ and $M_1$ (Take a trivial segment for each vertex and a segment
  for each edge, and a circle for a loop, and use the overlay definition
  of a set of segments and loops above). The faces of $M$ refer to the
  maximal connected open point sets of the spherical subdivision implied
  by the embedding of $G$. Each object from the output tuple $(V,E,F)$
  has a \emph{supporting} object $u_i$ in each of the two input
  structures.  Imagine the two maps to be transparent balls, where one
  contains the other. Then each point of the sphere is covered by an
  object from each of the input structures.  This support relationship
  from the input structures to the output structure defines an
  information flow. Each supporting object $u_i$ of $u$ $(i=0,1)$
  carries an associated information $|mark|(u_i)$. After the subdivision
  operation this information is attributed to the output object $u$ by
  $|mark|(u,i)$.}*/

  typedef Refs_ SNC_structure;
  typedef SNC_SM_decorator<Refs_> Base;
  typedef SNC_SM_decorator<Refs_> Decorator;
  typedef SNC_SM_overlayer<Refs_> Self;

#define USING(t) typedef typename Refs_::t t
  USING(Vertex_handle);
  USING(SVertex_handle);
  USING(SHalfedge_handle);
  USING(SHalfloop_handle);
  USING(SFace_handle);
  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SFace_iterator);
  USING(SObject_handle);
#undef USING
#define DECUSING(t) typedef typename Decorator::t t
  DECUSING(SHalfedge_around_svertex_circulator);
  DECUSING(SHalfedge_around_sface_circulator);
#undef DECUSING

  typedef std::pair<SHalfedge_handle,SHalfedge_handle> SHalfedge_pair;

  /*{\Mtypes 3}*/

  typedef typename Refs_::Sphere_kernel          Sphere_kernel;
  /*{\Mtypemember the geometry kernel.}*/
  typedef typename Sphere_kernel::Sphere_point   Sphere_point;
  /*{\Mtypemember the point type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember the segment type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_circle  Sphere_circle;
  /*{\Mtypemember the circle type of the sphere geometry.}*/
  typedef typename Refs_::Mark Mark;
  /*{\Mtypemember the mark of sphere map objects.}*/

  typedef typename Base::GenPtr GenPtr;

  /*{\Mgeneralization SNC_SM_decorator<Refs_>}*/

protected:
  Decorator PI[2];
  const Sphere_kernel& K;

public:

  // ---------------------------------------------------------------

  struct Seg_info { // to transport information from input to output
    SObject_handle _o; int _from;

    Seg_info() : _o(), _from(-1) {}
    Seg_info(SVertex_handle v, int i) 
    { _o=SObject_handle(v); _from=i; }
    Seg_info(SHalfedge_handle e, int i) 
    { _o=SObject_handle(e); _from=i; }
    Seg_info(SHalfloop_handle l, int i) 
    { _o=SObject_handle(l); _from=i; }
    Seg_info(const Seg_info& si) 
    { _o=si._o; _from=si._from; }
    Seg_info& operator=(const Seg_info& si) 
    { _o=si._o; _from=si._from; return *this; }
    LEDA_MEMORY(Seg_info)
  }; // Seg_info

  typedef std::list<Sphere_segment>            Seg_list;
  typedef typename Seg_list::iterator          Seg_iterator;
  typedef std::pair<Seg_iterator,Seg_iterator> Seg_it_pair;
  typedef std::pair<Sphere_segment,Sphere_segment> Seg_pair;
  typedef CGAL::Unique_hash_map<Seg_iterator,Seg_info> Seg_map;

  // ---------------------------------------------------------------

  struct vertex_info {
    Mark m[2];
    SObject_handle o_supp[2];
    SHalfedge_handle e_below;
    vertex_info() 
    { o_supp[0]=o_supp[1]=SObject_handle(); }
    LEDA_MEMORY(vertex_info)
  }; // vertex_info

  void assoc_info(SVertex_handle v) const
  { geninfo<vertex_info>::create(info(v)); }

  void discard_info(SVertex_handle v) const
  { geninfo<vertex_info>::clear(info(v)); }

  vertex_info& ginfo(SVertex_handle v) const
  { return geninfo<vertex_info>::access(info(v)); }

  Mark& mark(SVertex_handle v, int i) const
  { return ginfo(v).m[i]; }

  SObject_handle& supp_object(SVertex_handle v, int i) const
  { return ginfo(v).o_supp[i]; }

  SHalfedge_handle& halfedge_below(SVertex_handle v) const
  { return ginfo(v).e_below; }

  // ---------------------------------------------------------------

  struct edge_info {
    Mark m[2];
    Mark mf[2];
    SObject_handle o_supp[2];
    bool forw;
    edge_info()
    { m[0]=m[1]=mf[0]=mf[1]=Mark(); 
      o_supp[0]=o_supp[1]=SObject_handle(); 
      forw=false; }
    LEDA_MEMORY(edge_info)
  };

  void assoc_info(SHalfedge_handle e)  const
  { geninfo<edge_info>::create(info(e)); 
    geninfo<edge_info>::create(info(twin(e))); }

  void discard_info(SHalfedge_handle e)  const
  { geninfo<edge_info>::clear(info(e)); 
    geninfo<edge_info>::clear(info(twin(e))); }

  edge_info& ginfo(SHalfedge_handle e)  const
  { return geninfo<edge_info>::access(info(e)); }

  Mark& mark(SHalfedge_handle e, int i)  const
  // uedge information we store in the smaller one 
  { if (&*e < &*(twin(e))) return ginfo(e).m[i]; 
    else                   return ginfo(twin(e)).m[i]; }

  SObject_handle& supp_object(SHalfedge_handle e, int i) const
  // uedge information we store in the smaller one 
  { if (&*e < &*(twin(e))) return ginfo(e).o_supp[i]; 
    else                   return ginfo(twin(e)).o_supp[i]; }

  Mark& incident_mark(SHalfedge_handle e, int i)  const
  // biedge information we store in the edge
  { return ginfo(e).mf[i]; }

  bool& is_forward(SHalfedge_handle e) const
  // biedge information we store in the edge
  { return ginfo(e).forw; }

  // ---------------------------------------------------------------

  struct face_info {
    Mark m[2];
    face_info() { m[0]=m[1]=Mark(); }
    LEDA_MEMORY(face_info)
  };

  void assoc_info(SFace_handle f)  const
  { geninfo<face_info>::create(info(f)); }

  void discard_info(SFace_handle f)  const
  { geninfo<face_info>::clear(info(f)); }

  face_info& ginfo(SFace_handle f)  const
  { return geninfo<face_info>::access(info(f)); }

  Mark& mark(SFace_handle f, int i)  const
  { return ginfo(f).m[i]; }

  // ---------------------------------------------------------------
  // for documentation of determine face see the Nef_2 
  // implementation report

  template <typename Below_accessor>
  SFace_handle determine_face(SHalfedge_handle e, 
    const std::vector<SHalfedge_handle>& MinimalHalfedge,
    const CGAL::Unique_hash_map<SHalfedge_handle,int>& FaceCycle,
    const Below_accessor& D) const
  { TRACEN("determine_face "<<PH(e));
    int fc = FaceCycle[e];
    SHalfedge_handle e_min = MinimalHalfedge[fc];
    SHalfedge_handle e_below = D.halfedge_below(target(e_min));
    CGAL_nef3_assertion( e_below != SHalfedge_handle() );
    SFace_handle f = face(e_below);
    if ( f != SFace_handle() ) return f; // has already a face 
    // e_below also has no face
    f = determine_face(e_below, MinimalHalfedge, FaceCycle,D);
    link_as_face_cycle(e_below,f);
    return f;
  }

  Sphere_segment segment(Decorator N, 
                         SHalfedge_handle e) const
  { return K.construct_segment(
      N.point(N.source(e)),N.point(N.target(e)),N.circle(e)); }

  Sphere_segment trivial_segment(Decorator N, 
                                 SVertex_handle v) const
  { Sphere_point p = N.point(v); 
    return K.construct_segment(p,p); }

  Seg_pair two_segments(Decorator N, 
                        SHalfedge_handle e) const
  // we know that source(e)==target(e)
  { return N.circle(e).split_at(N.point(N.source(e))); }

  Seg_pair two_segments(Decorator N, 
                        SHalfloop_handle l) const
  { return N.circle(l).split_at_xy_plane(); }


  Mark& mark(SVertex_handle h) const
  { return Base::mark(h); }
  Mark& mark(SHalfedge_handle h) const
  { return Base::mark(h); }
  Mark& mark(SHalfloop_handle h) const
  { return Base::mark(h); }
  Mark& mark(SFace_handle h) const
  { return Base::mark(h); }


  // ---------------------------------------------------------------
  // INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE INTERFACE 
  // ---------------------------------------------------------------

  /*{\Mcreation 6}*/

  SNC_SM_overlayer(Vertex_handle v, 
    const Sphere_kernel& G = Sphere_kernel()) : Base(v), K(G) {}
  /*{\Mcreate |\Mvar| is a decorator object manipulating the map
  of |v|.}*/

  /*{\Moperations 1.1 1}*/

  template <typename Forward_iterator>
  void create_from_segments(
    Forward_iterator start, Forward_iterator end) const; 
  /*{\Mop produces the sphere map which is the overlay of the
  segments from the iterator range |[start,end)|.  \precond
  |Forward_iterator| has value type |Sphere_segment|.}*/

  template <typename Forward_iterator>
  void create_from_circles(
    Forward_iterator start, Forward_iterator end) const; 
  /*{\Mop produces the sphere map which is the overlay of the
  circles from the iterator range |[start,end)|.  \precond
  |Forward_iterator| has value type |Sphere_circle|.}*/

  void create(const Sphere_circle& c) const;
  /*{\Mop produces the sphere map which consists of one loop
  and the two halfspheres incident to it.}*/

  void subdivide(Vertex_handle v0, Vertex_handle v1);
  /*{\Mop constructs the overlay of the sphere maps of the vertices
  |v0| and |v1| in |M|, where all objects (svertices, shalfedges, sfaces)
  of |M| are \emph{enriched} by the marks of the supporting objects of
  the two input structures.}*/

  template <typename Selection> 
  void select(const Selection& SP) const;
  /*{\Mop sets the marks of all objects according to the selection
  predicate |SP|. |Selection| has to be a function object type with a
  function operator\\
  [[Mark operator()(Mark m0, Mark m1) const]]\\
  For each object |u| of |M| enriched by the marks of the supporting
  objects according to the previous procedure |subdivide|, after this
  operation |\Mvar.mark(u) = SP ( \Mvar.mark(u,0),\Mvar.mark(u,1)
  )|. The additional marks are invalidated afterwards. 
  \precond subdivide() was called before.}*/

  void simplify() const;
  /*{\Mop simplifies the structure of |M| according to the marks of
  its objects. An edge |e| separating two faces |f1| and |f2| and equal
  marks |mark(e) == mark(f1) == mark(f2)| is removed and the faces are
  unified.  An isolated vertex |v| in a face |f| with |mark(v)==mark(f)|
  is removed.  A vertex |v| with outdegree two, two collinear out-edges
  |e1|,|e2| and equal marks |mark(v) == mark(e1) == mark(e2)| is removed
  and the edges are unified.}*/

  template <typename Iterator, typename T>
  void partition_to_halfsphere(Iterator start, Iterator end,
    Seg_list& L, CGAL::Unique_hash_map<Iterator,T>& M, int pos) const;

  template <typename Mark_accessor>
  void merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2,
    const Mark_accessor& D) const;

  template <typename Mark_accessor>
  void merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2,
    const Mark_accessor& D) const;

  template <typename Below_accessor, typename Halfsphere_geometry>
  void create_face_objects(SHalfedge_iterator e_start, 
			   SHalfedge_iterator e_end,
			   SVertex_iterator v_start, SVertex_iterator v_end,
			   const Below_accessor& D, 
			   const Halfsphere_geometry& SG) const;

  template <typename Below_accessor>
  void complete_face_support(SVertex_iterator v_start, SVertex_iterator v_end,
    const Below_accessor& D, int pos) const;

  void dump(std::ostream& os = std::cerr) const
  { SNC_SM_io_parser<Refs_>::dump(center_vertex(),os); }

}; // SNC_SM_overlayer<Refs_>


template <typename Refs_>
template <typename Forward_iterator>
void SNC_SM_overlayer<Refs_>::
create_from_segments(Forward_iterator start, Forward_iterator end) const
{
  TRACEN("creating from segment iterator range");
  Seg_list L(start,end);
  Unique_hash_map<Seg_iterator,bool> From_input(false);
  Seg_iterator it;
  CGAL_nef3_forall_iterators(it,L) From_input[it]=true;
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input, -1);
  //TRACEN("L_pos="<<(MSDEBUG::print_elements(L_pos),""));
  //TRACEN("L_neg="<<(MSDEBUG::print_elements(L_neg),""));

  typedef SMO_from_segs<Self,Seg_iterator> SNC_SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_iterator v;
  SHalfedge_iterator e;
  SNC_SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--svertices_end(); e=--shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(shalfedges_begin(), e, svertices_begin(), v, O,
                      K.get_positive_halfsphere_geometry());
  create_face_objects(e, shalfedges_end(), v, svertices_end(), O,
                      K.get_negative_halfsphere_geometry());

  SHalfedge_iterator u;
  CGAL_nef3_forall_sedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(svertices_begin(),v,O);
  check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Refs_>
template <typename Forward_iterator>
void SNC_SM_overlayer<Refs_>::
create_from_circles(Forward_iterator start, Forward_iterator end) const
{
  TRACEN("creating from circle iterator range");
  Seg_list L;
  Unique_hash_map<Seg_iterator,bool> From_input(false);
  for ( ; start != end; ++start ) {
    std::pair<Sphere_segment,Sphere_segment> spair =
      start->split_at_xy_plane();
    L.push_back(spair.first); L.push_back(spair.second);
  }
  Seg_iterator it;
  CGAL_nef3_forall_iterators(it,L) From_input[it]=true;
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input, -1);


  typedef SMO_from_segs<Self,Seg_iterator> SNC_SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_iterator v;
  SHalfedge_iterator e;
  SNC_SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--svertices_end(); e=--shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(shalfedges_begin(), e, svertices_begin(), v, O,
                      K.get_positive_halfsphere_geometry());
  create_face_objects(e, shalfedges_end(), v, svertices_end(), O,
                      K.get_negative_halfsphere_geometry());

  SHalfedge_iterator u;
  CGAL_nef3_forall_sedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(svertices_begin(),v,O);
  check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Refs_>
void SNC_SM_overlayer<Refs_>::
create(const Sphere_circle& c) const
{ SHalfloop_handle l1 = new_loop_pair();
  SHalfloop_handle l2 = twin(l1);
  circle(l1) = c; circle(l2) = c.opposite();
  SFace_handle f1 = new_face();
  SFace_handle f2 = new_face();
  link_as_loop(l1,f1);
  link_as_loop(l2,f2);
}

template <typename Refs_>
void SNC_SM_overlayer<Refs_>::
subdivide(Vertex_handle v0, Vertex_handle v1)
{
  PI[0] = Decorator(v0); PI[1] = Decorator(v1);
  Seg_list L;
  Seg_map  From;
  for (int i=0; i<2; ++i) {
    SVertex_iterator v;
    CGAL_nef3_forall_svertices(v,PI[i]) {
      if ( !PI[i].is_isolated(v) ) continue;
      TRACEN("isolated " << PH(v));
      L.push_back(trivial_segment(PI[i],v));
      From[--L.end()] = Seg_info(v,i);
    }
    SHalfedge_iterator e;
    CGAL_nef3_forall_sedges(e,PI[i]) {
      if ( source(e) == target(e) ) {
	TRACEN("degenerierte Kante " << PH(e));
        Seg_pair p = two_segments(PI[i],e);
        L.push_back(p.first); 
        L.push_back(p.second);
        From[--L.end()] = From[--(--L.end())] = Seg_info(e,i);
      } else {
	TRACEN("normale Kante " << PH(e));
        L.push_back(segment(PI[i],e));
        From[--L.end()] = Seg_info(e,i);
      }
    }
    if ( PI[i].has_loop() ) {
      TRACEN("halfloop");
      SHalfloop_handle shl = PI[i].shalfloop();
      Seg_pair p = two_segments(PI[i],shl);
      L.push_back(p.first); 
      L.push_back(p.second);
      From[--L.end()] = From[--(--L.end())] = 
        Seg_info(shl,i);
      /*  
      p = two_segments(PI[i],PI[i].twin(shl));
      L.push_back(p.first); 
      L.push_back(p.second);
      From[--L.end()] = From[--(--L.end())] = 
        Seg_info(PI[i].twin(shl),i);
      */
    }
  }
  
  typename Seg_list::iterator it;
  CGAL_nef3_forall_iterators(it,L) TRACEN("  "<<*it);

  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From, -1);
  //TRACEN("L_pos="<<(MSDEBUG::print_elements(L_pos),""));
  //TRACEN("L_neg="<<(MSDEBUG::print_elements(L_neg),""));

  typedef SMO_from_sm<Self,Seg_iterator,Seg_info> SNC_SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SNC_SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_handle v;
  SHalfedge_handle e;
  SNC_SM_output O(*this,PI,From); 

  /* DEBUG CODE: to do: have all svertices a halfedge below associated? */
  TRACEN("Vertex info before swep");
  SVertex_iterator svi;
  for( svi=svertices_begin(); svi!=svertices_end(); svi++) {
    GenPtr i = info(svi);
    TRACEN("vertex "<<point(svi)<<" info "<<i<<
	   " marks "<<mark(svi,0)<<" "<<mark(svi,1));
  }

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  v=--svertices_end(); e=--shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(shalfedges_begin(), e, svertices_begin(), v, O,
                      PH_geometry());
  create_face_objects(e, shalfedges_end(), v, svertices_end(), O,
                      NH_geometry());


  SHalfedge_iterator u;
  CGAL_nef3_forall_sedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle(); 
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  complete_face_support(svertices_begin(), v, O, +1);
  complete_face_support(v, svertices_end(), O, -1);



  /* DEBUG CODE: to do: have all svertices a halfedge below associated? */
  TRACEN("Vertex info after swep");
  for( svi=svertices_begin(); svi!=svertices_end(); svi++) {
    GenPtr i = info(svi);
    TRACEN("vertex "<<point(svi)<<" info "<<i<<
	   " marks "<<mark(svi,0)<<" "<<mark(svi,1));
  }

  merge_halfsphere_maps(svertices_begin(),v,O);
  check_integrity_and_topological_planarity();

  TRACEN("subdivided");
  CGAL_nef3_forall_svertices(v,*this) {
    TRACEN(PH(v));
  }
}


template <typename Refs_>
template <typename Iterator, typename T>
void SNC_SM_overlayer<Refs_>::
partition_to_halfsphere(Iterator start, Iterator beyond, Seg_list& L, 
  CGAL::Unique_hash_map<Iterator,T>& M, int pos) const
{ TRACEN("partition_to_halfsphere ");
  CGAL_nef3_assertion(pos!=0);
  Sphere_segment s1,s2;
  Sphere_circle xycircle(0,0,pos);
  while ( start != beyond ) { 
    int i = start->intersection(xycircle,s1,s2);
    TRACEN("segment " << start->source() << " " << start->target());
    if (i>1) { L.push_back(s2); M[--L.end()] = M[start];    TRACEN(">1 " << s2.source() << " " << s2.target()); }
    if (i>0) { L.push_back(s1); M[--L.end()] = M[start];    TRACEN(">0 " << s1.source() << " " << s1.target()); }
    ++start;
  }
  // now all segments are split into hemispheres
  // we still have to:
  // - split segments containing our special poles y^-, y^+
  // - split halfcircles
  // - add four equator segments 
  Sphere_point S(0,-1,0),N(0,1,0);
  Sphere_circle yzcircle(1,0,0);
  typename Seg_list::iterator it, itl;
  
  CGAL_nef3_forall_iterators(it,L) { TRACEN("  "<<*it);
    if ( equal_as_sets(it->sphere_circle(),xycircle) ) {
      TRACEN("  splitting xy seg "<<*it);
      bool added=false;
      int n1 =  it->intersection(yzcircle,s1,s2);
      if (n1 > 1 && !s2.is_degenerate()) 
      { M[ L.insert(it,s2) ] = M[it]; added=true; TRACEN(">1 " << s2.source() << " " << s2.target()); }
      if (n1 > 0 && !s1.is_degenerate()) 
      { M[ L.insert(it,s1) ] = M[it]; added = true; TRACEN(">1 " << s1.source() << " " << s1.target()); }
      int n2 =  it->intersection(yzcircle.opposite(),s1,s2);
      if (n2 > 1 && !s2.is_degenerate()) 
      { M[ L.insert(it,s2) ] = M[it]; added=true; TRACEN(">1 " << s2.source() << " " << s2.target()); }
      if (n2 > 0 && !s1.is_degenerate()) 
      { M[ L.insert(it,s1) ] = M[it]; added=true; TRACEN(">1 " << s1.source() << " " << s1.target()); }
      if(added) {
	itl = it; --it; L.erase(itl); M[itl] = T();
      }
      // at least one item was appended
    }
  }
  CGAL_nef3_forall_iterators(it,L) {
    if ( it->is_halfcircle() ) {
      TRACEN("  splitting halfcircle "<<*it);
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
}



template <typename Refs_>
template <typename Below_accessor, typename Halfsphere_geometry>
void SNC_SM_overlayer<Refs_>::
create_face_objects(SHalfedge_iterator e_start, SHalfedge_iterator e_end,
  SVertex_iterator v_start, SVertex_iterator v_end,
  const Below_accessor& D, 
  const Halfsphere_geometry& SG) const
{
  TRACEN("create_face_objects()");
  CGAL::Unique_hash_map<SHalfedge_handle,int> FaceCycle(-1);
  std::vector<SHalfedge_handle>  MinimalHalfedge;
  SHalfedge_around_sface_circulator hfc(last_out_edge(v_start)),hend(hfc);
  TRACEN("equator cycle "<<PH(hfc));
  CGAL_For_all(hfc,hend) FaceCycle[hfc]=0; // outer face cycle = 0
  MinimalHalfedge.push_back(twin(first_out_edge(v_start)));
  int i=1; 
  for (SHalfedge_iterator e = e_start; e != e_end; ++e) {
    if ( FaceCycle[e] >= 0 ) continue; // already assigned
    SHalfedge_around_sface_circulator hfc(e),hend(hfc);
    SHalfedge_handle e_min = e;
    TRACEN("  face cycle numbering "<<i);
    CGAL_For_all(hfc,hend) {
      FaceCycle[hfc]=i; // assign face cycle number
      if ( SG.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      TRACE(PH(hfc));
    } TRACEN("");
    MinimalHalfedge.push_back(e_min);
    ++i;
  }

  for (int j=1; j<i; ++j) {
    SHalfedge_handle e = MinimalHalfedge[j];
    TRACEN("  face cycle "<<j<<" minimal halfedge "<<PH(e));
    Sphere_point p1 = point(source(e)), 
                 p2 = point(target(e)), 
                 p3 = point(target(next(e)));
    if ( SG.orientation(p1,p2,p3) > 0 ) { // left_turn => outer face cycle
      SFace_handle f = new_face();
      link_as_face_cycle(e,f);
      TRACEN("  creating new face object "<<&*f<<" bd "<<&*e);
    }
  }

  for (SHalfedge_iterator e = e_start; e != e_end; ++e) {
    if ( face(e) != SFace_handle() ) continue;
    if ( FaceCycle[e] == 0 ) continue;
    TRACEN("linking hole "<<PH(e));
    SFace_handle f = determine_face(e,MinimalHalfedge,FaceCycle,D);
    link_as_face_cycle(e,f);
  }
  for (SVertex_iterator v = v_start; v != v_end; ++v) {
    if ( !is_isolated(v) ) continue;
    SHalfedge_handle e_below = D.halfedge_below(v);
    CGAL_nef3_assertion( e_below != SHalfedge_handle() );
    link_as_isolated_vertex(v,face(e_below));
  }

}

template <typename Refs_>
template <typename Below_accessor>
void SNC_SM_overlayer<Refs_>::
complete_face_support(SVertex_iterator v_start, SVertex_iterator v_end,
  const Below_accessor& D, int pos) const
{ TRACEN("complete_face_support");
  for (SVertex_iterator v = v_start; v != v_end; ++v) { 
    TRACEN("VERTEX = "<<PH(v));
    Mark m_buffer[2];
    SHalfedge_handle e_below = halfedge_below(v);
    if ( v == v_start ) {     
      for (int i=0; i<2; ++i){ 
	SHalfedge_around_sface_circulator e(first_out_edge(v)), end(e);
	CGAL_For_all(e,end) {
	  if(supp_object(e,i) != NULL)
	    break;
	}
	if(supp_object(e,i) != NULL) {
	  SHalfedge_handle ei;
	  if ( assign(ei,supp_object(e,i)) ) { 
	    if ( PI[i].circle(ei) != circle(e) ) { ei = PI[i].twin(ei); }
	    CGAL_nef3_assertion( PI[i].circle(ei) == circle(e) ); 
	    TRACEN("initial " << PH(e) << " " << PH(ei)<< " von Ebene " << i);
	    m_buffer[i] = PI[i].mark(PI[i].face(ei));       
	  }
	  SHalfloop_handle li;
	  if ( assign(li,supp_object(e,i)) ) { 
	    if ( PI[i].circle(li) != circle(e) ) { li = PI[i].twin(li); }
	    m_buffer[i] = PI[i].mark(PI[i].face(li));
	    TRACEN("initial " << PH(li) << " von Ebene " << i);
	  }
	}
	else {
	  m_buffer[i] = PI[i].mark_of_halfsphere(-pos);
	  TRACEN("no initial support");
	}
      }
    } else if ( e_below != SHalfedge_handle() ) {
      for (int i=0; i<2; ++i) {
	TRACEN("edge below "<< PH(e_below) << " " << mark(e_below,i));
        m_buffer[i] = incident_mark(e_below,i); 
      }
    } else { // e_below does not exist
      CGAL_nef3_assertion( point(v).hz() == 0 && 
                      ( pos > 0 ? (point(v).hx() >= 0) : (point(v).hx()<=0)) );
      for (int i=0; i<2; ++i) 
        m_buffer[i] = incident_mark(previous(first_out_edge(v)),i);
    } TRACEN(" faces right and below "<<m_buffer[0]<<" "<<m_buffer[1]);

    for (int i=0; i<2; ++i) {
      SObject_handle o = supp_object(v,i);
      if ( o == NULL ) { TRACEN("no vertex support"); mark(v,i) = m_buffer[i]; continue; }
      SVertex_handle vs;
      SHalfedge_handle es;
      SHalfloop_handle ls;
      if ( assign(vs,o) ) { mark(v,i) = PI[i].mark(vs); continue; }
      if ( assign(es,supp_object(v,i)) ) {
        if ( point(source(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(source(es)); continue; }
        if ( point(target(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(target(es)); continue; }
        mark(v,i) = PI[i].mark(es); continue;
      }
      if ( assign(ls,o) ) { mark(v,i) = PI[i].mark(ls); TRACEN("loop " << PI[i].circle(ls)); continue; }
      CGAL_nef3_assertion_msg(0,"wrong handle");
    } TRACEN(" vertex marks "<<mark(v,0)<<" "<<mark(v,1));

    if ( is_isolated(v) ) continue;
    SHalfedge_around_svertex_circulator e(first_out_edge(v)), hend(e);
    CGAL_For_all(e,hend) {
      if ( !is_forward(e) ) break;
      TRACEN("  forward edge "<<PH(e));
      for (int i=0; i<2; ++i) {
        if ( supp_object(e,i) != NULL ) {
          SHalfedge_handle ei; 
          if ( assign(ei,supp_object(e,i)) ) { 
            if ( PI[i].circle(ei) != circle(e) ) { ei = PI[i].twin(ei); }
            CGAL_nef3_assertion( PI[i].circle(ei) == circle(e) ); 
            TRACEN("  supporting edge "<<i<<" "<<PH(ei));
            incident_mark(twin(e),i) =
              PI[i].mark(PI[i].face(PI[i].twin(ei)));
            mark(e,i) = PI[i].mark(ei);
            incident_mark(e,i) = m_buffer[i] =
              PI[i].mark(PI[i].face(ei)); 
          }
          SHalfloop_handle li;
          if ( assign(li,supp_object(e,i)) ) { 
            if ( PI[i].circle(li) != circle(e) ) { li = PI[i].twin(li); }
            CGAL_nef3_assertion( PI[i].circle(li) == circle(e) ); 
            TRACEN("  supporting loop "<<i<<" "<<PH(li));
            incident_mark(twin(e),i) =
              PI[i].mark(PI[i].face(PI[i].twin(li)));
            mark(e,i) = PI[i].mark(li);
            incident_mark(e,i) = m_buffer[i] =
              PI[i].mark(PI[i].face(li)); 
          }
        } else { TRACEN("  support from face below "<<i);
          incident_mark(twin(e),i) = mark(e,i) = 
          incident_mark(e,i) = m_buffer[i];
        }
      } TRACEN("  face marks "<<m_buffer[0]<<" "<<m_buffer[1]);
    }

    TRACEN(" mark of "<<PH(v)<<" "<<mark(v,0)<<" "<<mark(v,1));
  }

 
  SFace_iterator f;
  for (f = sfaces_begin(); f != sfaces_end(); ++f) {
    assoc_info(f);
    SObject_handle boundary_object = sface_cycles_begin(f);
    CGAL_nef3_assertion(boundary_object != NULL);
    SHalfedge_handle e;
    if ( !CGAL::assign(e,boundary_object) ) 
      CGAL_nef3_assertion_msg(0,"Outer face cycle should be first.");
    for (int i=0; i<2; ++i) mark(f,i) = incident_mark(e,i);
  }

  TRACEN(psm_->point());

  SVertex_handle v;
  CGAL_nef3_forall_svertices(v,*this)
    TRACEN(PH(v) << " " << mark(v,0));
  TRACEN(" ");
  CGAL_nef3_forall_svertices(v,PI[0])
    TRACEN(PH(v));
  TRACEN(" ");
  CGAL_nef3_forall_svertices(v,*this)
    TRACEN(PH(v) << " " << mark(v,1));
  TRACEN(" ");
  CGAL_nef3_forall_svertices(v,PI[1])
    TRACEN(PH(v));
  TRACEN(" ");

  SHalfedge_handle e;
  CGAL_nef3_forall_shalfedges(e,*this)
    TRACEN(PH(e)<< " " << mark(e,0) << " " << incident_mark(e,0));
  TRACEN(" ");
  CGAL_nef3_forall_shalfedges(e,PI[0])
    TRACEN(PH(e)<<  "|" << circle(e) <<"|" << PI[0].mark(e) << " " << PI[0].mark(PI[0].face(e)));
  TRACEN(" ");
  CGAL_nef3_forall_shalfedges(e,*this)
    TRACEN(PH(e) << " " << mark(e,1) << " " << incident_mark(e,1));
  TRACEN(" ");
  CGAL_nef3_forall_shalfedges(e,PI[1])
    TRACEN(PH(e) <<  "|" << circle(e) <<"|" << PI[1].mark(e) << " " << PI[1].mark(PI[1].face(e)));
  TRACEN(" ");

    SFace_handle ff;
  CGAL_nef3_forall_sfaces(ff,*this)
    TRACEN(&*ff << " " << mark(ff,0));
  TRACEN(" ");
    CGAL_nef3_forall_sfaces(ff,PI[0])
    TRACEN(&*ff << " " << PI[0].mark(ff));
  TRACEN(" ");
   CGAL_nef3_forall_sfaces(ff,*this)
    TRACEN(&*ff << " " << mark(ff,1));
  TRACEN(" "); 
    CGAL_nef3_forall_sfaces(ff,PI[1])
    TRACEN(&*ff << " " << PI[1].mark(ff));
    TRACEN(" ");

}

template <typename Refs_>
template <typename Mark_accessor>
void SNC_SM_overlayer<Refs_>::
merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2,
  const Mark_accessor& D) const
{
  SVertex_handle v1 = source(e1), v2 = target(e2);
  TRACEN("merge_nodes "<<PH(v1)<<PH(v2));
  CGAL_nef3_assertion(point(v1)==point(v2));
  SHalfedge_handle ep1 = previous(e1), en2 = next(e2);
  SHalfedge_around_svertex_circulator eav(out_edges(v2)),ee(eav);
  CGAL_For_all(eav,ee) { set_source(eav,v1); }
  link_as_prev_next_pair(e2,e1);  
  link_as_prev_next_pair(ep1,en2); 
  //D.assert_equal_marks(v1,v2);
  mark(v1,0) = mark(v1,0) || mark(v2,0);
  mark(v1,1) = mark(v1,1) || mark(v2,1);
  D.discard_info(v2);
  delete_vertex_only(v2);
}

template <typename Refs_>
template <typename Mark_accessor>
void SNC_SM_overlayer<Refs_>::
merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2,
  const Mark_accessor& D) const
{ TRACEN("merging halfspheres "<<PH(v1)<<PH(v2));
  CGAL_nef3_assertion(point(v1)==point(v2));
  std::list<SHalfedge_pair> L_equator;
  SHalfedge_around_sface_circulator 
    ep(last_out_edge(v1)), en(twin(first_out_edge(v2)));
  do { 
   L_equator.push_back(SHalfedge_pair(ep,en));
   merge_nodes(ep,en,D); ++ep; --en; 
  } while ( source(ep) != v1 );
  
  typename std::list<SHalfedge_pair>::iterator it;
  CGAL_nef3_forall_iterators(it,L_equator) { 
    SHalfedge_handle e1 = it->first, e2 = it->second;
    SHalfedge_handle e1t = twin(e1), e2t = twin(e2);
    TRACEV(PH(e1));TRACEV(PH(e2));
    SHalfedge_handle e2tp = previous(e2t);
    SHalfedge_handle e2tn = next(e2t);
    link_as_prev_next_pair(e2tp,e1);
    link_as_prev_next_pair(e1,e2tn);
    SFace_handle f = face(e2t);
    if ( is_boundary_object(e2t) )
    { undo_boundary_object(e2t,f); store_boundary_object(e1,f); }
    set_face(e1,f);
    if ( e2 == first_out_edge(source(e2)) )
      set_first_out_edge(source(e2),e1t);
    mark(e1,0) = mark(e1,0) || mark(e2,0);
    mark(e1,1) = mark(e1,1) || mark(e2,1);
    D.discard_info(e2);
    delete_edge_pair_only(e2);
  }
}

template <typename Refs_>
template <typename Selection>
void SNC_SM_overlayer<Refs_>::
select(const Selection& SP) const
{ 
  SVertex_iterator v;
  CGAL_nef3_forall_svertices(v,*this) {
    mark(v) = SP(mark(v,0),mark(v,1));
    discard_info(v); 
  }
  SHalfedge_iterator e;
  CGAL_nef3_forall_sedges(e,*this) {
    mark(e) = SP(mark(e,0),mark(e,1));
    discard_info(e);
  }
  SFace_iterator f;
  CGAL_nef3_forall_sfaces(f,*this) {
    mark(f) = SP(mark(f,0),mark(f,1));
    discard_info(f);
  }
  mark_of_halfsphere(-1) = SP(PI[0].mark_of_halfsphere(-1),
                              PI[1].mark_of_halfsphere(-1));
  mark_of_halfsphere(+1) = SP(PI[0].mark_of_halfsphere(+1),
                              PI[1].mark_of_halfsphere(+1));
}

template <typename Refs_>
void SNC_SM_overlayer<Refs_>::simplify() const
{

  TRACEN("simplifying"); 

  SVertex_handle vy;
  CGAL_nef3_forall_svertices(vy,*this)
    TRACEN(PH(vy)); // << " " << mark(vy,0) << " " << mark(vy,1));
  TRACEN(" ");

  SHalfedge_handle ey;
  CGAL_nef3_forall_shalfedges(ey,*this)
    TRACEN(PH(ey) << " " << mark(face(ey))); //  << " " << mark(ey,1));
  TRACEN(" ");

  /* typedef typename CGAL::Partition<SFace_handle>::item partition_item;
     CGAL::Unique_hash_map<SFace_handle,partition_item> Pitem;
     CGAL::Partition<SFace_handle> FP; */
  typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
  CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem;
  CGAL::Union_find< SFace_handle> UF;

  SFace_iterator f;
  CGAL_nef3_forall_sfaces(f,*this) {
     Pitem[f] = UF.make_set(f);
     clear_face_cycle_entries(f);
  }

  SHalfedge_iterator e, en;
  for(e = shalfedges_begin(); e != shalfedges_end(); e = en) { 
    en = e; ++en; if ( en==twin(e) ) ++en;
    SNC_decorator<Refs_> D;
    TRACEN("can simplify ? " << PH(e));
    if(!D.is_sedge_on_infibox(e)) {
      TRACEN(mark(e) << " " << mark(face(e)) << " " << mark(face(twin(e))));
      if (( mark(e) == mark(face(e)) && mark(e) == mark(face(twin(e))))){
	TRACEN("deleting "<<PH(e));
	if ( !UF.same_set(Pitem[face(e)],
			  Pitem[face(twin(e))]) ) {
	  
	  UF.unify_sets( Pitem[face(e)],
			 Pitem[face(twin(e))] );
	  TRACEN("unioning disjoint faces");
	}
	if ( is_closed_at_source(e) )
	  set_face(source(e),face(e));
	if ( is_closed_at_source(twin(e)))
	  set_face(target(e),face(e));
	delete_edge_pair(e);
      }
    }
  }

  CGAL_nef3_forall_svertices(vy,*this)
    TRACEN(PH(vy)); // << " " << mark(vy,0) << " " << mark(vy,1));
  TRACEN(" ");
    
  CGAL_nef3_forall_shalfedges(ey,*this)
    TRACEN(PH(ey)); //<< " " << mark(ey,0) << " " << mark(ey,1));
  TRACEN(" ");

  CGAL::Unique_hash_map<SHalfedge_handle,bool> linked(false);
  for (e = shalfedges_begin(); e != shalfedges_end(); ++e) {
    if ( linked[e] ) continue;
    SHalfedge_around_sface_circulator hfc(e),hend(hfc);
    SFace_handle f = *(UF.find( Pitem[face(e)]));
    CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
    store_boundary_object(e,f);
  }
  if ( has_loop() ) {
    SHalfloop_handle l = shalfloop();
    SFace_handle f = *(UF.find(Pitem[face(l)]));
    link_as_loop(l,f);
    f = *(UF.find(Pitem[face(twin(l))]));
    link_as_loop(twin(l),f);
  }

  SVertex_iterator v,vn;
  for(v = svertices_begin(); v != svertices_end(); v=vn) {
    vn=v; ++vn;
    if ( is_isolated(v) ) {
      if ( mark(v) == mark(face(v)) ) {
        TRACEN("removing isolated vertex"<<PH(v));
        delete_vertex_only(v);  
      } else
        store_boundary_object(v,face(v)); // isolated, but should stay
    } else { // v not isolated
      SHalfedge_handle e2 = first_out_edge(v), e1 = previous(e2);
      if ( has_outdeg_two(v) &&
           mark(v) == mark(e1) && mark(v) == mark(e2) &&
           circle(e1) == circle(e2) ) {
        TRACEN("collinear at "<<PH(v)<<PH(e1)<<PH(e2));
        if ( e1 == e2 ) convert_edge_to_loop(e1);
        else  merge_edge_pairs_at_target(e1); 
      }
    }
  }

  SFace_iterator fn;
  for (f = fn = sfaces_begin(); f != sfaces_end(); f=fn) { 
    ++fn;
    Union_find_handle pit = Pitem[f];
    if ( UF.find(pit) != pit ) 
      delete_face_only(f);
  }

  SVertex_handle vx;
  CGAL_nef3_forall_svertices(vx,*this)
    TRACEN(PH(vx)); // << " " << mark(vx,0) << " " << mark(vx,1));
  TRACEN(" ");

  SHalfedge_handle ex;
  CGAL_nef3_forall_shalfedges(ex,*this)
    TRACEN(PH(ex)); // << " " << mark(ex,0) << " " << mark(ex,1));
  TRACEN(" ");
}



CGAL_END_NAMESPACE
#endif //CGAL_SNC_SM_OVERLAYER_H


