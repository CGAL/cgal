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
// file          : include/CGAL/Nef_S2/SM_overlayer.h
// package       : Nef_S2 
// chapter       : Nef Polyhedra
//
// source        : nef_s2/SM_overlayer.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Overlay module for sphere maps
// ============================================================================

#ifndef CGAL_SM_OVERLAYER_H
#define CGAL_SM_OVERLAYER_H

#include <CGAL/basic.h>
#include <CGAL/Union_find.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#undef _DEBUG
#define _DEBUG 131
#include <CGAL/Nef_S2/debug.h>

#define USING(t) typedef typename Base::t t
#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif
CGAL_BEGIN_NAMESPACE

template <typename Decorator_, typename I>
struct SMO_from_segs {
  typedef Decorator_ SM_decorator;
  typedef typename Decorator_::Vertex_handle  Vertex_handle;
  typedef typename Decorator_::Halfedge_handle    Halfedge_handle;
  typedef typename Decorator_::Sphere_point    Point;
  typedef typename Decorator_::Sphere_segment  Segment;
  typedef CGAL::Unique_hash_map<I,bool>  Iterator_map;
  SM_decorator G;
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
    G.new_edge_pair_at_source(v,SM_decorator::BEFORE); 
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
  { CGAL_nef_assertion(G.mark(v1)==G.mark(v2)); }
  void discard_info(Vertex_handle v) const 
  { geninfo<Halfedge_handle>::clear(G.info(v)); }

  void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
  { CGAL_nef_assertion(G.mark(e1)==G.mark(e2)); }
  void transfer_marks(Halfedge_handle e) const 
  { G.unify_tmp_marks(e); }
  void discard_info(Halfedge_handle e) const {}

  void clear_temporary_vertex_info() const
  { Vertex_handle v;
    for(v = G.vertices_begin(); v != G.vertices_end(); ++v)
      geninfo<Halfedge_handle>::clear(G.info(v));
  }


}; // SMO_from_segs


template <typename Decorator_, typename IT, typename INFO>
struct SMO_from_sm {
  typedef Decorator_ SM_overlayer;
  typedef typename Decorator_::SM_decorator SM_decorator;
  typedef typename SM_decorator::Vertex_handle Vertex_handle;
  typedef typename SM_decorator::Halfedge_handle   Halfedge_handle;
  typedef typename SM_decorator::Sphere_point   Point;
  typedef typename SM_decorator::Sphere_segment Segment;
  SM_overlayer G;
  SM_decorator* pGI;
  CGAL::Unique_hash_map<IT,INFO>& M;
  SMO_from_sm(SM_overlayer Gi, 
              SM_decorator* pGIi, 
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
  G.new_edge_pair_at_source(v,SM_decorator::BEFORE); 
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
  TRACEN("   supporting "<<si._from<<" "<<*it);
}


void trivial_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  CGAL_nef_assertion( si._o != NULL );
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
  CGAL_nef_assertion(G.mark(v1,0)==G.mark(v2,0)&&
                 G.mark(v1,1)==G.mark(v2,1)); }
void discard_info(Vertex_handle v) const 
{ G.discard_info(v); }

void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
{ CGAL_nef_assertion(G.mark(e1,0)==G.mark(e2,0) && 
                 G.mark(e1,1)==G.mark(e2,1)); }
void transfer_marks(Halfedge_handle e) const 
{ Halfedge_handle et = G.twin(e);
  if (&*e < &*et) std::swap(e,et);
  for(int i=0; i<2; ++i) G.ginfo(e).m[i] = G.ginfo(et).m[i];
}

void discard_info(Halfedge_handle e) const 
{ G.discard_info(e); }




}; // SMO_from_sm

template <typename Decorator_, typename ITERATOR>
class SMO_decorator { 
public:
  typedef Decorator_ Graph;
  typedef typename Decorator_::Vertex_handle  Vertex_handle;
  typedef typename Decorator_::Halfedge_handle    Halfedge_handle;
  typedef typename Decorator_::Sphere_point    Point_2;
  typedef typename Decorator_::Sphere_segment  Segment_2;
  Decorator_ G;

SMO_decorator(Graph Gi) : G(Gi) {}

Vertex_handle new_vertex(const Point_2& p)
{ return G.new_vertex(p); }

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e)
{ G.link_as_target_and_append(v,e); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
{ return G.new_edge_pair_at_source(v,Graph::BEFORE); }

void supporting_segment(Halfedge_handle e, ITERATOR it) {}
void halfedge_below(Vertex_handle v, Halfedge_handle e) {}
void trivial_segment(Vertex_handle v, ITERATOR it) {}
void starting_segment(Vertex_handle v, ITERATOR it) {}
void passing_segment(Vertex_handle v, ITERATOR it) {}
void ending_segment(Vertex_handle v, ITERATOR it) {}


}; // SMO_decorator





/*{\Manpage {SM_overlayer}{Decorator_}{Overlay in the sphere}{O}}*/

template <typename Decorator_>
class SM_overlayer : public Decorator_ {
public:
  /*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
  decorator object offering sphere map overlay calculation. Overlay is
  either calculated from two sphere maps or from a set of halfspaces.
  The result is stored in a sphere map |M| that carries the geometry and
  the topology of the overlay.

  The template parameter provides the underlying topological interface
  to sphere maps. The template parameter |Decorator_| has to be a model
  conforming to our map decorator concept |SM_decorator|.  The concept
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
  of |SM_decorator|.} including isolated vertices. The overlay process
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

  typedef Decorator_                 Base;
  typedef typename Base::Sphere_map  Sphere_map;
  typedef SM_overlayer<Decorator_>   Self;
  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halfloop_handle);
  USING(Face_handle);
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Face_iterator);
  USING(Object_handle);
  USING(Halfedge_around_vertex_circulator);
  USING(Halfedge_around_face_circulator);
  typedef std::pair<Halfedge_handle,Halfedge_handle> Halfedge_pair;

  /*{\Mtypes 3}*/
  typedef Base SM_decorator;
  /*{\Mtypemember the sphere map decorator.}*/

  typedef typename Base::Kernel Kernel;
  /*{\Mtypemember the geometry kernel.}*/

  typedef typename Kernel::Sphere_point   Sphere_point;
  /*{\Mtypemember the point type of the sphere geometry.}*/
  typedef typename Kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember the segment type of the sphere geometry.}*/
  typedef typename Kernel::Sphere_circle  Sphere_circle;
  /*{\Mtypemember the circle type of the sphere geometry.}*/
  typedef typename SM_decorator::Mark       Mark;
  /*{\Mtypemember the mark of sphere map objects.}*/

  /*{\Mgeneralization SM_decorator}*/

protected:
  SM_decorator PI[2];
  const Kernel& K;

public:
  struct Seg_info { // to transport information from input to output
    Object_handle _o; int _from;

    Seg_info() : _o(), _from(-1) {}
    Seg_info(Vertex_handle v, int i) 
    { _o=Object_handle(v); _from=i; }
    Seg_info(Halfedge_handle e, int i) 
    { _o=Object_handle(e); _from=i; }
    Seg_info(Halfloop_handle l, int i) 
    { _o=Object_handle(l); _from=i; }
    Seg_info(const Seg_info& si) 
    { _o=si._o; _from=si._from; }
    Seg_info& operator=(const Seg_info& si) 
    { _o=si._o; _from=si._from; return *this; }
    LEDA_MEMORY(Seg_info)
  };

  typedef std::list<Sphere_segment>            Seg_list;
  typedef typename Seg_list::iterator          Seg_iterator;
  typedef std::pair<Seg_iterator,Seg_iterator> Seg_it_pair;
  typedef std::pair<Sphere_segment,Sphere_segment> Seg_pair;
  typedef CGAL::Unique_hash_map<Seg_iterator,Seg_info> Seg_map;

  struct vertex_info {
    Mark m[2];
    Object_handle o_supp[2];
    Halfedge_handle e_below;
    vertex_info() 
    { o_supp[0]=o_supp[1]=Object_handle(); }
    LEDA_MEMORY(vertex_info)
  };

  void assoc_info(Vertex_handle v) const
  { geninfo<vertex_info>::create(info(v)); }

  void discard_info(Vertex_handle v) const
  { geninfo<vertex_info>::clear(info(v)); }

  vertex_info& ginfo(Vertex_handle v) const
  { return geninfo<vertex_info>::access(info(v)); }

  Mark& mark(Vertex_handle v, int i) const
  { return ginfo(v).m[i]; }

  Object_handle& supp_object(Vertex_handle v, int i) const
  { return ginfo(v).o_supp[i]; }

  Halfedge_handle& halfedge_below(Vertex_handle v) const
  { return ginfo(v).e_below; }

  struct edge_info {
    Mark m[2];
    Mark mf[2];
    Object_handle o_supp[2];
    bool forw;
    edge_info()
    { m[0]=m[1]=mf[0]=mf[1]=Mark(); 
      o_supp[0]=o_supp[1]=Object_handle(); 
      forw=false; }
    LEDA_MEMORY(edge_info)
  };

  void assoc_info(Halfedge_handle e)  const
  { geninfo<edge_info>::create(info(e)); 
    geninfo<edge_info>::create(info(twin(e))); }

  void discard_info(Halfedge_handle e)  const
  { geninfo<edge_info>::clear(info(e)); 
    geninfo<edge_info>::clear(info(twin(e))); }

  edge_info& ginfo(Halfedge_handle e)  const
  { return geninfo<edge_info>::access(info(e)); }

  Mark& mark(Halfedge_handle e, int i)  const
  // uedge information we store in the smaller one 
  { if (&*e < &*(twin(e))) return ginfo(e).m[i]; 
    else                   return ginfo(twin(e)).m[i]; }

  Object_handle& supp_object(Halfedge_handle e, int i) const
  // uedge information we store in the smaller one 
  { if (&*e < &*(twin(e))) return ginfo(e).o_supp[i]; 
    else                   return ginfo(twin(e)).o_supp[i]; }

  Mark& incident_mark(Halfedge_handle e, int i)  const
  // biedge information we store in the edge
  { return ginfo(e).mf[i]; }

  bool& is_forward(Halfedge_handle e) const
  // biedge information we store in the edge
  { return ginfo(e).forw; }

  struct face_info {
    Mark m[2];
    face_info() { m[0]=m[1]=Mark(); }
    LEDA_MEMORY(face_info)
  };

  void assoc_info(Face_handle f)  const
  { geninfo<face_info>::create(info(f)); }

  void discard_info(Face_handle f)  const
  { geninfo<face_info>::clear(info(f)); }

  face_info& ginfo(Face_handle f)  const
  { return geninfo<face_info>::access(info(f)); }

  Mark& mark(Face_handle f, int i)  const
  { return ginfo(f).m[i]; }


  template <typename Below_accessor>
  Face_handle determine_face(Halfedge_handle e, 
    const std::vector<Halfedge_handle>& MinimalHalfedge,
    const CGAL::Unique_hash_map<Halfedge_handle,int>& FaceCycle,
    const Below_accessor& D) const
  { TRACEN("determine_face "<<PH(e));
    int fc = FaceCycle[e];
    Halfedge_handle e_min = MinimalHalfedge[fc];
    Halfedge_handle e_below = D.halfedge_below(target(e_min));
    CGAL_nef_assertion( e_below != Halfedge_handle() );
    Face_handle f = face(e_below);
    if ( f != Face_handle() ) return f; // has already a face 
    // e_below also has no face
    f = determine_face(e_below, MinimalHalfedge, FaceCycle,D);
    link_as_face_cycle(e_below,f);
    return f;
  }


  Sphere_segment segment(SM_decorator N, 
                         Halfedge_handle e) const
  { return K.construct_segment(
      N.point(N.source(e)),N.point(N.target(e)),N.circle(e)); }

  Sphere_segment trivial_segment(SM_decorator N, 
                                 Vertex_handle v) const
  { Sphere_point p = N.point(v); 
    return K.construct_segment(p,p); }

  Seg_pair two_segments(SM_decorator N, 
                        Halfedge_handle e) const
  // we know that source(e)==target(e)
  { return N.circle(e).split_at(N.point(N.source(e))); }

  Seg_pair two_segments(SM_decorator N, 
                        Halfloop_handle l) const
  { return N.circle(l).split_at_xy_plane(); }


  Mark& mark(Vertex_handle h) const
  { return Base::mark(h); }
  Mark& mark(Halfedge_handle h) const
  { return Base::mark(h); }
  Mark& mark(Halfloop_handle h) const
  { return Base::mark(h); }
  Mark& mark(Face_handle h) const
  { return Base::mark(h); }



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

  void subdivide(const Sphere_map& M0, const Sphere_map& M1);
  /*{\Mop constructs the overlay of the sphere maps |M0| and |M1| in
  |M|, where all objects (vertices, halfedges, faces) of |M| are
  \emph{enriched} by the marks of the supporting objects of the two
  input structures: e.g. let |v| be a vertex supported by a node |v0| in
  |M0| and by a face |f1| in |M1| and |D0|, |D1| be decorators of
  type |SM_decorator| on |M0|,|M1|. Then |\Mvar.mark(v,0) = D0.mark(v0)|
  and |\Mvar.mark(v,1) = D1.mark(f1)|.}*/

  template <typename Selection> 
  void select(const Selection& SP) const;
  /*{\Mop sets the marks of all objects according to the selection
  predicate |SP|. |Selection| has to be a function object type with a
  function operator\\
  [[Mark operator()(Mark m0, Mark m1) const]]\\
  For each object |u| of |M| enriched by the marks of the supporting
  objects according to the previous procedure |subdivide|, after this
  operation |\Mvar.mark(u) = SP ( \Mvar.mark(u,0),\Mvar.mark(u,1)
  )|. The additional marks are invalidated afterwards. }*/

  void simplify() const;
  /*{\Mop simplifies the structure of |M| according to the marks of
  its objects. An edge |e| separating two faces |f1| and |f2| and equal
  marks |mark(e) == mark(f1) == mark(f2)| is removed and the faces are
  unified.  An isolated vertex |v| in a face |f| with |mark(v)==mark(f)|
  is removed.  A vertex |v| with outdegree two, two collinear out-edges
  |e1|,|e2| and equal marks |mark(v) == mark(e1) == mark(e2)| is removed
  and the edges are unified.}*/

  /*{\Mcreation 6}*/
  SM_overlayer(SM_decorator D, 
    const Kernel& G = Kernel()) : Base(D), K(G) {}
  /*{\Mcreate |\Mvar| is a decorator object manipulating the map
  decorated by |D|.}*/

  /*{\Moperations 1.1 1}*/
  template <typename Iterator>
  void subdivide_segments(Iterator start, Iterator end) const;
  template <typename Iterator, typename T>
  void partition_to_halfsphere(Iterator start, Iterator end,
    Seg_list& L, CGAL::Unique_hash_map<Iterator,T>& M, int pos) const;
  template <typename Mark_accessor>
  void merge_halfsphere_maps(Vertex_handle v1, Vertex_handle v2,
    const Mark_accessor& D) const;
  template <typename Mark_accessor>
  void merge_nodes(Halfedge_handle e1, Halfedge_handle e2,
    const Mark_accessor& D) const;

  template <typename Below_accessor, typename Halfsphere_geometry>
  void create_face_objects(Halfedge_iterator e_start, Halfedge_iterator e_end,
    Vertex_iterator v_start, Vertex_iterator v_end,
    const Below_accessor& D, 
    const Halfsphere_geometry& SG) const;
  template <typename Below_accessor>
  void complete_face_support(Vertex_iterator v_start, Vertex_iterator v_end,
    const Below_accessor& D, int pos) const;

  void dump(std::ostream& os = std::cerr) const
  { SM_io_parser<Base>::dump(*this,os); }

}; // SM_overlayer<Decorator_>

template <typename Decorator_>
template <typename Forward_iterator>
void SM_overlayer<Decorator_>::
create_from_segments(Forward_iterator start, Forward_iterator end) const
{
  TRACEN("creating from segment iterator range");
  Seg_list L(start,end);
  Unique_hash_map<Seg_iterator,bool> From_input(false);
  Seg_iterator it;
  CGAL_forall_iterators(it,L) From_input[it]=true;
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input, -1);
  //TRACEN("L_pos="<<(MSDEBUG::print_elements(L_pos),""));
  //TRACEN("L_neg="<<(MSDEBUG::print_elements(L_neg),""));

  typedef SMO_from_segs<Self,Seg_iterator> SM_output;
  typedef typename Kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  Vertex_iterator v;
  Halfedge_iterator e;
  SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--vertices_end(); e=--halfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(halfedges_begin(), e, vertices_begin(), v, O,
                      K.get_positive_halfsphere_geometry());
  create_face_objects(e, halfedges_end(), v, vertices_end(), O,
                      K.get_negative_halfsphere_geometry());

  Halfedge_iterator u;
  CGAL_forall_edges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(vertices_begin(),v,O);
  check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Decorator_>
template <typename Forward_iterator>
void SM_overlayer<Decorator_>::
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
  CGAL_forall_iterators(it,L) From_input[it]=true;
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input, -1);


  typedef SMO_from_segs<Self,Seg_iterator> SM_output;
  typedef typename Kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  Vertex_iterator v;
  Halfedge_iterator e;
  SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--vertices_end(); e=--halfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(halfedges_begin(), e, vertices_begin(), v, O,
                      K.get_positive_halfsphere_geometry());
  create_face_objects(e, halfedges_end(), v, vertices_end(), O,
                      K.get_negative_halfsphere_geometry());

  Halfedge_iterator u;
  CGAL_forall_edges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(vertices_begin(),v,O);
  check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Decorator_>
void SM_overlayer<Decorator_>::
create(const Sphere_circle& c) const
{ Halfloop_handle l1 = new_loop_pair();
  Halfloop_handle l2 = twin(l1);
  circle(l1) = c; circle(l2) = c.opposite();
  Face_handle f1 = new_face();
  Face_handle f2 = new_face();
  link_as_loop(l1,f1);
  link_as_loop(l2,f2);
}

template <typename Decorator_>
void SM_overlayer<Decorator_>::
subdivide(const Sphere_map& M0, const Sphere_map& M1)
{
  PI[0] = SM_decorator(const_cast<Sphere_map&>(M0)); 
  PI[1] = SM_decorator(const_cast<Sphere_map&>(M1));
  Seg_list L;
  Seg_map  From;
  for (int i=0; i<2; ++i) {
    Vertex_iterator v;
    CGAL_forall_vertices(v,PI[i]) {
      if ( !PI[i].is_isolated(v) ) continue;
      L.push_back(trivial_segment(PI[i],v));
      From[--L.end()] = Seg_info(v,i);
    }
    Halfedge_iterator e;
    CGAL_forall_edges(e,PI[i]) {
      if ( source(e) == target(e) ) {
        Seg_pair p = two_segments(PI[i],e);
        L.push_back(p.first); 
        L.push_back(p.second);
        From[--L.end()] = From[--(--L.end())] = Seg_info(e,i);
      } else {
        L.push_back(segment(PI[i],e));
        From[--L.end()] = Seg_info(e,i);
      }
    }
    if ( PI[i].has_loop() ) {
      Seg_pair p = two_segments(PI[i],PI[i].halfloop());
      L.push_back(p.first); 
      L.push_back(p.second);
      From[--L.end()] = From[--(--L.end())] = 
        Seg_info(PI[i].halfloop(),i);
    }
  }

  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From, +1);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From, -1);
  //TRACEN("L_pos="<<(MSDEBUG::print_elements(L_pos),""));
  //TRACEN("L_neg="<<(MSDEBUG::print_elements(L_neg),""));

  typedef SMO_from_sm<Self,Seg_iterator,Seg_info> SM_output;
  typedef typename Kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  Vertex_handle v;
  Halfedge_handle e;
  SM_output O(*this,PI,From); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    K.get_positive_halfsphere_geometry());
  SP.sweep();
  v=--vertices_end(); e=--halfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    K.get_negative_halfsphere_geometry());
  SM.sweep();
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere
   
  create_face_objects(halfedges_begin(), e, vertices_begin(), v, O,
                      PH_geometry());
  create_face_objects(e, halfedges_end(), v, vertices_end(), O,
                      NH_geometry());

  Halfedge_iterator u;
  CGAL_forall_edges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle(); 
    circle(twin(u)) = s.sphere_circle().opposite();
  }



  complete_face_support(vertices_begin(), v, O, +1);
  complete_face_support(v, vertices_end(), O, -1);

  merge_halfsphere_maps(vertices_begin(),v,O);
  check_integrity_and_topological_planarity();

}


template <typename Decorator_>
template <typename Iterator, typename T>
void SM_overlayer<Decorator_>::
partition_to_halfsphere(Iterator start, Iterator beyond, Seg_list& L, 
  CGAL::Unique_hash_map<Iterator,T>& M, int pos) const
{ TRACEN("partition_to_halfsphere ");
  CGAL_nef_assertion(pos!=0);
  Sphere_segment s1,s2;
  Sphere_circle xycircle(0,0,pos);
  while ( start != beyond ) { 
    int i = start->intersection(xycircle,s1,s2);
    if (i>1) { L.push_back(s2); M[--L.end()] = M[start]; }
    if (i>0) { L.push_back(s1); M[--L.end()] = M[start]; }
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
  
  CGAL_forall_iterators(it,L) { TRACEN("  "<<*it);
    if ( equal_as_sets(it->sphere_circle(),xycircle) ) {
      TRACEN("  splitting xy seg "<<*it);
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
    }
  }
  CGAL_forall_iterators(it,L) {
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

template <typename Decorator_>
template <typename Below_accessor, typename Halfsphere_geometry>
void SM_overlayer<Decorator_>::
create_face_objects(Halfedge_iterator e_start, Halfedge_iterator e_end,
  Vertex_iterator v_start, Vertex_iterator v_end,
  const Below_accessor& D, 
  const Halfsphere_geometry& SG) const
{
  TRACEN("create_face_objects()");
  CGAL::Unique_hash_map<Halfedge_handle,int> FaceCycle(-1);
  std::vector<Halfedge_handle>  MinimalHalfedge;
  Halfedge_around_face_circulator hfc(last_out_edge(v_start)),hend(hfc);
  TRACEN("equator cycle "<<PH(hfc));
  CGAL_For_all(hfc,hend) FaceCycle[hfc]=0; // outer face cycle = 0
  MinimalHalfedge.push_back(twin(first_out_edge(v_start)));
  int i=1; 
  for (Halfedge_iterator e = e_start; e != e_end; ++e) {
    if ( FaceCycle[e] >= 0 ) continue; // already assigned
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Halfedge_handle e_min = e;
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
    Halfedge_handle e = MinimalHalfedge[j];
    TRACEN("  face cycle "<<j<<" minimal halfedge "<<PH(e));
    Sphere_point p1 = point(source(e)), 
                 p2 = point(target(e)), 
                 p3 = point(target(next(e)));
    if ( SG.orientation(p1,p2,p3) > 0 ) { // leftturn => outer face cycle
      Face_handle f = new_face();
      link_as_face_cycle(e,f);
      TRACEN("  creating new face object "<<&*f<<" bd "<<&*e);
    }
  }

  for (Halfedge_iterator e = e_start; e != e_end; ++e) {
    if ( face(e) != Face_handle() ) continue;
    if ( FaceCycle[e] == 0 ) continue;
    TRACEN("linking hole "<<PH(e));
    Face_handle f = determine_face(e,MinimalHalfedge,FaceCycle,D);
    link_as_face_cycle(e,f);
  }
  for (Vertex_iterator v = v_start; v != v_end; ++v) {
    if ( !is_isolated(v) ) continue;
    Halfedge_handle e_below = D.halfedge_below(v);
    CGAL_nef_assertion( e_below != Halfedge_handle() );
    link_as_isolated_vertex(v,face(e_below));
  }

}

template <typename Decorator_>
template <typename Below_accessor>
void SM_overlayer<Decorator_>::
complete_face_support(Vertex_iterator v_start, Vertex_iterator v_end,
  const Below_accessor& D, int pos) const
{ TRACEN("complete_face_support");
  for (Vertex_iterator v = v_start; v != v_end; ++v) { 
    TRACEN("VERTEX = "<<PH(v));
    Mark m_buffer[2];
    Halfedge_handle e_below = halfedge_below(v);
    if ( v == v_start ) {
      for (int i=0; i<2; ++i) 
        m_buffer[i] = PI[i].mark_of_halfsphere(-pos);
    } else if ( e_below != Halfedge_handle() ) {
      for (int i=0; i<2; ++i) 
        m_buffer[i] = incident_mark(e_below,i); 
    } else { // e_below does not exist
      CGAL_nef_assertion( point(v).hz() == 0 && 
                      ( pos > 0 ? (point(v).hx() >= 0) : (point(v).hx()<=0)) );
      for (int i=0; i<2; ++i) 
        m_buffer[i] = incident_mark(previous(first_out_edge(v)),i);
    } TRACEN(" faces right and below "<<m_buffer[0]<<" "<<m_buffer[1]);

    for (int i=0; i<2; ++i) {
      Object_handle o = supp_object(v,i);
      if ( o == NULL ) { mark(v,i) = m_buffer[i]; continue; }
      Vertex_handle vs;
      Halfedge_handle es;
      Halfloop_handle ls;
      if ( assign(vs,o) ) { mark(v,i) = PI[i].mark(vs); continue; }
      if ( assign(es,supp_object(v,i)) ) {
        if ( point(source(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(source(es)); continue; }
        if ( point(target(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(target(es)); continue; }
        mark(v,i) = PI[i].mark(es); continue;
      }
      if ( assign(ls,o) ) { mark(v,i) = PI[i].mark(ls); continue; }
      CGAL_nef_assertion_msg(0,"wrong handle");
    } TRACEN(" vertex marks "<<mark(v,0)<<" "<<mark(v,1));

    if ( is_isolated(v) ) continue;
    Halfedge_around_vertex_circulator e(first_out_edge(v)), hend(e);
    CGAL_For_all(e,hend) {
      if ( !is_forward(e) ) break;
      TRACEN("  forward edge "<<PH(e));
      for (int i=0; i<2; ++i) {
        if ( supp_object(e,i) != NULL ) {
          Halfedge_handle ei; 
          if ( assign(ei,supp_object(e,i)) ) { 
            if ( PI[i].circle(ei) != circle(e) ) { ei = PI[i].twin(ei); }
            CGAL_nef_assertion( PI[i].circle(ei) == circle(e) ); 
            TRACEN("  supporting edge "<<i<<" "<<PH(ei));
            incident_mark(twin(e),i) = 
              PI[i].mark(PI[i].face(PI[i].twin(ei)));
            mark(e,i) = PI[i].mark(ei);
            incident_mark(e,i) = m_buffer[i] =
              PI[i].mark(PI[i].face(ei)); 
          }
          Halfloop_handle li;
          if ( assign(li,supp_object(e,i)) ) { 
            if ( PI[i].circle(li) != circle(e) ) { li = PI[i].twin(li); }
            CGAL_nef_assertion( PI[i].circle(li) == circle(e) ); 
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
  Face_iterator f;
  for (f = faces_begin(); f != faces_end(); ++f) {
    assoc_info(f);
    Object_handle boundary_object = face_cycles_begin(f);
    CGAL_nef_assertion(boundary_object != NULL);
    Halfedge_handle e;
    if ( !assign(e,boundary_object) ) 
      CGAL_nef_assertion_msg(0,"Outer face cycle should be first.");
    for (int i=0; i<2; ++i) mark(f,i) = incident_mark(e,i);
  }


}

template <typename Decorator_>
template <typename Mark_accessor>
void SM_overlayer<Decorator_>::
merge_nodes(Halfedge_handle e1, Halfedge_handle e2,
  const Mark_accessor& D) const
{
  Vertex_handle v1 = source(e1), v2 = target(e2);
  TRACEN("merge_nodes "<<PH(v1)<<PH(v2));
  CGAL_nef_assertion(point(v1)==point(v2));
  Halfedge_handle ep1 = previous(e1), en2 = next(e2);
  Halfedge_around_vertex_circulator eav(out_edges(v2)),ee(eav);
  CGAL_For_all(eav,ee) { set_source(eav,v1); }
  link_as_prev_next_pair(e2,e1);  
  link_as_prev_next_pair(ep1,en2); 
  D.assert_equal_marks(v1,v2);
  D.discard_info(v2);
  delete_vertex_only(v2);
}

template <typename Decorator_>
template <typename Mark_accessor>
void SM_overlayer<Decorator_>::
merge_halfsphere_maps(Vertex_handle v1, Vertex_handle v2,
  const Mark_accessor& D) const
{ TRACEN("merging halfspheres "<<PH(v1)<<PH(v2));
  CGAL_nef_assertion(point(v1)==point(v2));
  std::list<Halfedge_pair> L_equator;
  Halfedge_around_face_circulator 
    ep(last_out_edge(v1)), en(twin(first_out_edge(v2)));
  do { 
   L_equator.push_back(Halfedge_pair(ep,en));
   merge_nodes(ep,en,D); ++ep; --en; 
  } while ( source(ep) != v1 );
  
  typename std::list<Halfedge_pair>::iterator it;
  CGAL_forall_iterators(it,L_equator) { 
    Halfedge_handle e1 = it->first, e2 = it->second;
    Halfedge_handle e1t = twin(e1), e2t = twin(e2);
    TRACEV(PH(e1));TRACEV(PH(e2));
    Halfedge_handle e2tp = previous(e2t);
    Halfedge_handle e2tn = next(e2t);
    link_as_prev_next_pair(e2tp,e1);
    link_as_prev_next_pair(e1,e2tn);
    Face_handle f = face(e2t);
    if ( is_boundary_object(e2t) )
    { undo_boundary_object(e2t,f); store_boundary_object(e1,f); }
    set_face(e1,f);
    if ( e2 == first_out_edge(source(e2)) )
      set_first_out_edge(source(e2),e1t);
    D.discard_info(e2);
    delete_edge_pair_only(e2);
  }
}

template <typename Decorator_>
template <typename Selection>
void SM_overlayer<Decorator_>::
select(const Selection& SP) const
{ 
  Vertex_iterator v;
  for(v = vertices_begin(); v != vertices_end(); ++v) {
    mark(v) = SP(mark(v,0),mark(v,1));
    discard_info(v); 
  }
  Halfedge_iterator e;
  for(e = halfedges_begin(); e != halfedges_end(); ++e) {
    if ( info(e) == 0 ) continue; // twin 
    mark(e) = SP(mark(e,0),mark(e,1));
    discard_info(e);
  }
  Face_iterator f;
  for(f = faces_begin(); f != faces_end(); ++f) {
    mark(f) = SP(mark(f,0),mark(f,1));
    discard_info(f);
  }
  mark_of_halfsphere(-1) = SP(PI[0].mark_of_halfsphere(-1),
                              PI[1].mark_of_halfsphere(-1));
  mark_of_halfsphere(+1) = SP(PI[0].mark_of_halfsphere(+1),
                              PI[1].mark_of_halfsphere(+1));
}

template <typename Decorator_>
void SM_overlayer<Decorator_>::simplify() const
{
  TRACEN("simplifying"); 
  typedef typename CGAL::Union_find<Face_handle>::handle Union_find_handle;
  CGAL::Unique_hash_map< Face_handle, Union_find_handle> Pitem;
  CGAL::Union_find<Face_handle> union_faces;
  Face_iterator f;
  for (f = faces_begin(); f != faces_end(); ++f) {
     Pitem[f] = union_faces.make_set(f);
     clear_face_cycle_entries(f);
  }

  Halfedge_iterator e, en;
  for(e = halfedges_begin(); e != halfedges_end(); e = en) { 
    en = e; ++en; if ( en==twin(e) ) ++en;
    if ( mark(e) == mark(face(e)) &&
         mark(e) == mark(face(twin(e))) ) {
      TRACEN("deleting "<<PH(e));
      if ( ! union_faces.same_set( Pitem[face(e)],
                                   Pitem[face(twin(e))]) ) {
        union_faces.unify_sets( Pitem[face(e)],
                                Pitem[face(twin(e))] );
        TRACEN("unioning disjoint faces");
      }
      if ( is_closed_at_source(e) ) 
        set_face(source(e),face(e));
      if ( is_closed_at_source(twin(e)) ) 
        set_face(target(e),face(e));
      delete_edge_pair(e);
    }
  }

  CGAL::Unique_hash_map<Halfedge_handle,bool> linked(false);
  for (e = halfedges_begin(); e != halfedges_end(); ++e) {
    if ( linked[e] ) continue;
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Face_handle f = *(union_faces.find( Pitem[face(e)]));
    CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
    store_boundary_object(e,f);
  }
  if ( has_loop() ) {
    Halfloop_handle l = halfloop();
    Face_handle f = *(union_faces.find( Pitem[face(l)]));
    link_as_loop(l,f);
    f = *(union_faces.find( Pitem[face(twin(l))]));
    link_as_loop(twin(l),f);
  }

  Vertex_iterator v,vn;
  for(v = vertices_begin(); v != vertices_end(); v=vn) {
    vn=v; ++vn;
    if ( is_isolated(v) ) {
      if ( mark(v) == mark(face(v)) ) {
        TRACEN("removing isolated vertex"<<PH(v));
        delete_vertex_only(v);  
      } else
        store_boundary_object(v,face(v)); // isolated, but should stay
    } else { // v not isolated
      Halfedge_handle e2 = first_out_edge(v), e1 = previous(e2);
      if ( has_outdeg_two(v) &&
           mark(v) == mark(e1) && mark(v) == mark(e2) &&
           circle(e1) == circle(e2) ) {
        TRACEN("collinear at "<<PH(v)<<PH(e1)<<PH(e2));
        if ( e1 == e2 ) convert_edge_to_loop(e1);
        else  merge_edge_pairs_at_target(e1); 
      }
    }
  }

  Face_iterator fn;
  for (f = fn = faces_begin(); f != faces_end(); f=fn) { 
    ++fn;
    Union_find_handle pit = Pitem[f];
    if ( union_faces.find(pit) != pit ) 
      delete_face_only(f);
  }


}

template <typename Decorator_>
template <typename Iterator>
void SM_overlayer<Decorator_>::
subdivide_segments(Iterator start, Iterator end) const
{
typedef SMO_decorator<SM_decorator,Iterator> SM_output;
typedef typename Kernel::Positive_halfsphere_geometry PH_geometry;
typedef CGAL::Segment_overlay_traits< 
          Iterator, SM_output, PH_geometry>  PHS_traits;
typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

typedef typename Kernel::Negative_halfsphere_geometry NH_geometry;
typedef CGAL::Segment_overlay_traits< 
          Iterator, SM_output, NH_geometry> NHS_traits;
typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  std::list<Sphere_segment> Lp,Lm;
  partition_xy( start, end, Lp , +1);
  partition_xy( start, end, Lm , -1);
  // both lists initialized with four quarter segments
  // supporting the xy-equator thereby separating the 
  // two halfspheres 
  // all other segments in the range are split into their
  // connected components with respect to the xy-plane.

  //TRACEN("Lp="<<(MSDEBUG::print_elements(Lp),""));
  //TRACEN("Lm="<<(MSDEBUG::print_elements(Lm),""));

  Vertex_handle v1,v2;
  SM_output O(*this); 
  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(Input_range(Lp.begin(),Lp.end()),O);
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v1= vertices_begin(); v2=--vertices_end();
  Negative_halfsphere_sweep SM(Input_range(Lm.begin(),Lm.end()),O);
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v2;
  // now two CCs of sphere graph calculated
  // v1 = first node of CC in positive xy-sphere
  // v2 = first node of CC in negative xy-sphere

  merge_halfsphere_maps(v1,v2,O);
  check_integrity_and_topological_planarity(false);
}


CGAL_END_NAMESPACE
#undef USING
#endif //CGAL_SM_OVERLAYER_H


