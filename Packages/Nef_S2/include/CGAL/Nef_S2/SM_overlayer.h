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

#ifndef CGAL_SM_OVERLAYER_H
#define CGAL_SM_OVERLAYER_H

#include <CGAL/basic.h>
#include <CGAL/Union_find.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_2/geninfo.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#undef _DEBUG
#define _DEBUG 131
#include <CGAL/Nef_S2/debug.h>

#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

CGAL_BEGIN_NAMESPACE

template <typename Decorator_, typename I>
struct SMO_from_segs {
  typedef Decorator_ SM_decorator;
  typedef typename SM_decorator::SVertex_handle     Vertex_handle;
  typedef typename SM_decorator::SHalfedge_handle   Halfedge_handle;
  typedef typename SM_decorator::Sphere_point       Point;
  typedef typename SM_decorator::Sphere_segment     Segment;
  typedef CGAL::Unique_hash_map<I,bool>             Iterator_map;
  SM_decorator G;
  const Iterator_map& M;
  SMO_from_segs(SM_decorator Gi, const Iterator_map& Mi) : G(Gi),M(Mi) {}

  Vertex_handle new_vertex(const Point& p)
  { Vertex_handle v = G.new_svertex(p); 
    geninfo<Halfedge_handle>::create(G.info(v));
    return v;
  }

  void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) 
  { G.link_as_target_and_append(v,e); }

  Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
  { Halfedge_handle e = 
    G.new_shalfedge_pair_at_source(v,SM_decorator::BEFORE); 
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
  { CGAL_assertion(G.mark(v1)==G.mark(v2)); }

  void discard_info(Vertex_handle v) const 
  { geninfo<Halfedge_handle>::clear(G.info(v)); }

  void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
  { CGAL_assertion(G.mark(e1)==G.mark(e2)); }

  void discard_info(Halfedge_handle e) const {}

  void clear_temporary_vertex_info() const
  { Vertex_handle v;
    CGAL_forall_svertices(v,G)
      geninfo<Halfedge_handle>::clear(G.info(v));
  }


}; // SMO_from_segs


template <typename SM_overlayer, typename IT, typename INFO>
struct SMO_from_sm {
  typedef typename SM_overlayer::SM_const_decorator      SM_const_decorator;
  typedef typename SM_overlayer::SVertex_handle          Vertex_handle;
  typedef typename SM_overlayer::SHalfedge_handle        Halfedge_handle;
  typedef typename SM_overlayer::Sphere_point            Point;
  typedef typename SM_overlayer::Sphere_segment          Segment;

  SM_overlayer G;
  SM_const_decorator* pGI;
  CGAL::Unique_hash_map<IT,INFO>& M;
  SMO_from_sm(SM_overlayer Gi, 
              SM_const_decorator* pGIi, 
              CGAL::Unique_hash_map<IT,INFO>& Mi) : 
    G(Gi), pGI(pGIi), M(Mi) {}

Vertex_handle new_vertex(const Point& p)
{ Vertex_handle v = G.new_svertex(p);
  G.assoc_info(v);
  return v;
}

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e)
{ G.link_as_target_and_append(v,e); }

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
{ Halfedge_handle e = 
  G.new_shalfedge_pair_at_source(v,SM_overlayer::BEFORE); 
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
  CGAL_assertion( si._o != NULL );
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
  CGAL_assertion(G.mark(v1,0)==G.mark(v2,0)&&
		 G.mark(v1,1)==G.mark(v2,1)); }
void discard_info(Vertex_handle v) const 
{ G.discard_info(v); }

void assert_equal_marks(Halfedge_handle e1, Halfedge_handle e2) const
{ CGAL_assertion(G.mark(e1,0)==G.mark(e2,0) && 
		 G.mark(e1,1)==G.mark(e2,1)); }

void discard_info(Halfedge_handle e) const 
{ G.discard_info(e); }

}; // SMO_from_sm

template <typename SM_decorator, typename ITERATOR>
class SMO_decorator { 
public:
  typedef SM_decorator Graph;
  typedef typename SM_decorator::SVertex_handle  SVertex_handle;
  typedef typename SM_decorator::SHalfedge_handle    SHalfedge_handle;
  typedef typename SM_decorator::Sphere_point    Point_2;
  typedef typename SM_decorator::Sphere_segment  Segment_2;
  SM_decorator G;

SMO_decorator(Graph Gi) : G(Gi) {}

SVertex_handle new_vertex(const Point_2& p)
{ return G.snew_vertex(p); }

void link_as_target_and_append(SVertex_handle v, SHalfedge_handle e)
{ G.link_as_target_and_append(v,e); }

SHalfedge_handle new_halfedge_pair_at_source(SVertex_handle v)
{ return G.new_shalfedge_pair_at_source(v,Graph::BEFORE); }

void supporting_segment(SHalfedge_handle e, ITERATOR it) {}
void halfedge_below(SVertex_handle v, SHalfedge_handle e) {}
void trivial_segment(SVertex_handle v, ITERATOR it) {}
void starting_segment(SVertex_handle v, ITERATOR it) {}
void passing_segment(SVertex_handle v, ITERATOR it) {}
void ending_segment(SVertex_handle v, ITERATOR it) {}


}; // SMO_decorator

// ============================================================================
// ============================================================================

/*{\Manpage {SM_overlayer}{SM_decorator}{Overlay in the sphere}{O}}*/

template <typename SM_decorator_>
class SM_overlayer : public SM_decorator_ {
public:

  /*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
  decorator object offering sphere map overlay calculation. Overlay is
  either calculated from two sphere maps or from a set of halfspaces.
  The result is stored in a sphere map |M| that carries the geometry and
  the topology of the overlay.

  The template parameter provides the underlying topological interface
  to sphere maps. The template parameter |SM_decorator| has to be a model
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
  edges.  SFaces are bounded by possibly several face cycles\footnote{For
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

  typedef SM_decorator_                         SM_decorator;
  typedef typename SM_decorator::Map            Map;
  typedef SM_decorator                          Base;
  typedef SM_overlayer<SM_decorator_>           Self;
  typedef CGAL::SM_const_decorator<Map>         SM_const_decorator;
  typedef SM_point_locator<SM_const_decorator>  SM_point_locator;

  //  typedef typename SM_const_decorator::Constructor_parameter 
  //                                       Constructor_const_parameter;
  typedef typename SM_const_decorator::SVertex_const_handle SVertex_const_handle;
  typedef typename SM_const_decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SM_const_decorator::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename SM_const_decorator::SFace_const_handle SFace_const_handle;
  typedef typename SM_const_decorator::SVertex_const_iterator SVertex_const_iterator;
  typedef typename SM_const_decorator::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename SM_const_decorator::SFace_const_iterator SFace_const_iterator;

  //  typedef typename Base::Constructor_parameter Constructor_parameter;
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

  typedef typename Base::Sphere_kernel           Sphere_kernel;
  /*{\Mtypemember the geometry kernel.}*/
  typedef typename Sphere_kernel::Sphere_point   Sphere_point;
  /*{\Mtypemember the point type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember the segment type of the sphere geometry.}*/
  typedef typename Sphere_kernel::Sphere_circle  Sphere_circle;
  /*{\Mtypemember the circle type of the sphere geometry.}*/
  typedef typename Base::Mark Mark;
  /*{\Mtypemember the mark of sphere map objects.}*/

  typedef typename Base::GenPtr GenPtr;

  /*{\Mgeneralization SM_decorator}*/

protected:
  SM_const_decorator PI[2];
  const Sphere_kernel& K;

public:

  // ---------------------------------------------------------------

  struct Seg_info { // to transport information from input to output
    Object_handle _o; int _from;

    Seg_info() : _o(), _from(-1) {}
    Seg_info(SVertex_const_handle v, int i) 
    { _o=Object_handle(v); _from=i; }
    Seg_info(SHalfedge_const_handle e, int i) 
    { _o=Object_handle(e); _from=i; }
    Seg_info(SHalfloop_const_handle l, int i) 
    { _o=Object_handle(l); _from=i; }
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
    Object_handle o_supp[2];
    SHalfedge_handle e_below;
    vertex_info() 
    { o_supp[0]=o_supp[1]=Object_handle(); }
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

  Object_handle& supp_object(SVertex_handle v, int i) const
  { return ginfo(v).o_supp[i]; }

  SHalfedge_handle& halfedge_below(SVertex_handle v) const
  { return ginfo(v).e_below; }

  // ---------------------------------------------------------------

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

  void assoc_info(SHalfedge_handle e)  const
  { geninfo<edge_info>::create(info(e)); 
    geninfo<edge_info>::create(info(twin(e))); }

  void discard_info(SHalfedge_handle e)  const
  { geninfo<edge_info>::clear(info(e)); 
    geninfo<edge_info>::clear(info(twin(e))); }

  edge_info& ginfo(SHalfedge_handle e)  const
  { return geninfo<edge_info>::access(info(e)); }

  Mark& mark(SHalfedge_handle e, int i)  const
    { return ginfo(e).m[i]; }

  Object_handle& supp_object(SHalfedge_handle e, int i) const
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

  template <typename Below_accessor>
  SFace_handle determine_face(SHalfedge_handle e, 
    const std::vector<SHalfedge_handle>& MinimalSHalfedge,
    const CGAL::Unique_hash_map<SHalfedge_handle,int>& SFaceCycle,
    const Below_accessor& D)
  { TRACEN("determine_face "<<PH(e));
    int fc = SFaceCycle[e];
    SHalfedge_handle e_min = MinimalSHalfedge[fc];
    SHalfedge_handle e_below = D.halfedge_below(target(e_min));
    if(e_below == SHalfedge_handle())
      return SFace_handle();
    SFace_handle f = face(e_below);
    if ( f != SFace_handle() ) return f; // has already a face 
    // e_below also has no face
    f = determine_face(e_below, MinimalSHalfedge, SFaceCycle,D);
    if(f != SFace_handle())
      link_as_face_cycle(e_below,f);
    return f;
  }

  Sphere_segment segment(SM_const_decorator N, 
                         SHalfedge_const_handle e) const
  { return K.construct_segment(
      N.point(N.source(e)),N.point(N.target(e)),N.circle(e)); }

  Sphere_segment trivial_segment(SM_const_decorator N, 
                                 SVertex_const_handle v) const
  { Sphere_point p = N.point(v); 
    return K.construct_segment(p,p); }

  Seg_pair two_segments(SM_const_decorator N, 
                        SHalfedge_const_handle e) const
  // we know that source(e)==target(e)
  { return N.circle(e).split_at(N.point(N.source(e))); }

  Seg_pair two_segments(SM_const_decorator N, 
                        SHalfloop_const_handle l) const
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

  SM_overlayer(Map* M, 
    const Sphere_kernel& G = Sphere_kernel()) : Base(M), K(G) {}
  /*{\Mcreate |\Mvar| is a decorator object manipulating the map
  of |v|.}*/

  /*{\Moperations 1.1 1}*/

  template <typename Forward_iterator>
  void create_from_segments(
    Forward_iterator start, Forward_iterator end); 
  /*{\Mop produces the sphere map which is the overlay of the
  segments from the iterator range |[start,end)|.  \precond
  |Forward_iterator| has value type |Sphere_segment|.}*/

  template <typename Forward_iterator>
  void create_from_circles(Forward_iterator start, Forward_iterator end);
  /*{\Mop produces the sphere map which is the overlay of the
  circles from the iterator range |[start,end)|.  \precond
  |Forward_iterator| has value type |Sphere_circle|.}*/

  void create(const Sphere_circle& c);
  /*{\Mop produces the sphere map which consists of one loop
  and the two halfspheres incident to it.}*/

  void subdivide(const Map* M0, const Map* M1);
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
  )|. The additional marks are invalidated afterwards.
  \precond subdivide() was called before.}*/

  void simplify();
  /*{\Mop simplifies the structure of |M| according to the marks of
  its objects. An edge |e| separating two faces |f1| and |f2| and equal
  marks |mark(e) == mark(f1) == mark(f2)| is removed and the faces are
  unified.  An isolated vertex |v| in a face |f| with |mark(v)==mark(f)|
  is removed.  A vertex |v| with outdegree two, two collinear out-edges
  |e1|,|e2| and equal marks |mark(v) == mark(e1) == mark(e2)| is removed
  and the edges are unified.}*/

  int check_sphere(const Seg_list& L, bool compute_halfsphere[3][2]) const;

  template <typename Iterator>
  void subdivide_segments(Iterator start, Iterator end) const;
  template <typename Iterator, typename T>
  void partition_to_halfsphere(Iterator start, Iterator end,
    Seg_list& L, CGAL::Unique_hash_map<Iterator,T>& M, 
    Sphere_circle xycircle, Sphere_circle yzcircle, bool include_equator) const;

  template <typename Mark_accessor>
  void merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2,
    const Mark_accessor& D);
  template <typename Mark_accessor>
  void merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2,
    const Mark_accessor& D);

  template <typename Below_accessor, typename Halfsphere_geometry>
  void create_face_objects(SHalfedge_iterator e_start, SHalfedge_iterator e_end,
			   SVertex_iterator v_start, SVertex_iterator v_end,
			   const Below_accessor& D, 
			   const Halfsphere_geometry& SG);

  template <typename Below_accessor>
  void complete_face_support(SVertex_iterator v_start, SVertex_iterator v_end,
    const Below_accessor& D, std::vector<Mark>& mohs, int offset) const;

  void set_outer_face_mark(int offset, const std::vector<Mark>& mohs);

  void dump(std::ostream& os = std::cerr) const
  { SM_io_parser<Base>::dump(*this,os); }

}; // SM_overlayer<SM_decorator>

template <typename Map>
template <typename Forward_iterator>
void SM_overlayer<Map>::
create_from_segments(Forward_iterator start, Forward_iterator end)
{
  TRACEN("creating from segment iterator range");
  Seg_list L(start,end);
  Unique_hash_map<Seg_iterator,bool> From_input(false);
  Seg_iterator it;
  CGAL_forall_iterators(it,L) From_input[it]=true;
  Seg_list L_pos,L_neg;
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input, 
			  Sphere_circle(0,0,1), Sphere_circle(1,0,0), true);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input, 
			  Sphere_circle(0,0,-1), Sphere_circle(1,0,0), true);
  //TRACEN("L_pos="<<(MSDEBUG::print_elements(L_pos),""));
  //TRACEN("L_neg="<<(MSDEBUG::print_elements(L_neg),""));

  typedef SMO_from_segs<Self,Seg_iterator> SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_iterator v;
  SHalfedge_iterator e;
  SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    PH_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--this->svertices_end(); e=--this->shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()),O,
    NH_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(this->shalfedges_begin(), e, this->svertices_begin(), v, O,
                      PH_geometry());
  create_face_objects(e, this->shalfedges_end(), v, this->svertices_end(), O,
                      NH_geometry());

  SHalfedge_iterator u;
  CGAL_forall_sedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(this->svertices_begin(),v,O);
  this->check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Map>
template <typename Forward_iterator>
void SM_overlayer<Map>::
create_from_circles(Forward_iterator start, Forward_iterator end)
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
  partition_to_halfsphere(L.begin(), L.end(), L_pos, From_input,
			  Sphere_circle(0,0,1), Sphere_circle(1,0,0), true);
  partition_to_halfsphere(L.begin(), L.end(), L_neg, From_input,
			  Sphere_circle(0,0,-1), Sphere_circle(1,0,0), true);

  typedef SMO_from_segs<Self,Seg_iterator> SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  typedef CGAL::Segment_overlay_traits< 
            Seg_iterator, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  SVertex_iterator v;
  SHalfedge_iterator e;
  SM_output O(*this,From_input); 

  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(
    Input_range(L_pos.begin(),L_pos.end()),O,
    PH_geometry());
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v=--this->svertices_end(); e=--this->shalfedges_end();

  Negative_halfsphere_sweep SM(
    Input_range(L_neg.begin(),L_neg.end()), O,
    NH_geometry());
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v; ++e;
  // now two CCs of sphere graph are calculated
  // v = first vertex of CC in negative x-sphere
  // e = first edge of CC in negative x-sphere

  create_face_objects(this->shalfedges_begin(), e, this->svertices_begin(), v, O,
                      PH_geometry());
  create_face_objects(e, this->shalfedges_end(), v, this->svertices_end(), O,
                      NH_geometry());

  SHalfedge_iterator u;
  CGAL_forall_sedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle();
    circle(twin(u)) = s.sphere_circle().opposite();
  }

  merge_halfsphere_maps(this->svertices_begin(),v,O);
  this->check_integrity_and_topological_planarity();

  O.clear_temporary_vertex_info();
}

template <typename Map>
int SM_overlayer<Map>::
check_sphere(const Seg_list& L, bool compute_halfsphere[3][2]) const {

  for(int i=0; i<6; i++)
    compute_halfsphere[i/2][i%2] = false;

  typename Seg_list::const_iterator it;
  CGAL_forall_iterators(it,L) {
    if(!compute_halfsphere[0][0])
      if(it->source().hx()>0 || it->target().hx()>0)
	compute_halfsphere[0][0] = true;
    if(!compute_halfsphere[0][1])
      if(it->source().hx()<0 || it->target().hx()<0)
	compute_halfsphere[0][1] = true;    
    if(!compute_halfsphere[1][0])
      if(it->source().hy()>0 || it->target().hy()>0)
	compute_halfsphere[1][0] = true;
    if(!compute_halfsphere[1][1])
      if(it->source().hy()<0 || it->target().hy()<0)
	compute_halfsphere[1][1] = true;    
    if(!compute_halfsphere[2][0])
      if(it->source().hz()>0 || it->target().hz()>0)
	compute_halfsphere[2][0] = true;
    if(!compute_halfsphere[2][1])
      if(it->source().hz()<0 || it->target().hz()<0)
	compute_halfsphere[2][1] = true;    
  }

  if(!compute_halfsphere[2][0]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(0,0,1),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hz()==0 && it->target().hz()==0) || it->is_long()) { 
	  compute_halfsphere[2][0] = true;
	  break;
	}
    }
  }
  
  if(!compute_halfsphere[2][0]) {
    compute_halfsphere[2][1] = true;
    return 4;
  }
  
  if(!compute_halfsphere[2][1]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(0,0,-1),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hz()==0 && it->target().hz()==0) || it->is_long()) { 
	  compute_halfsphere[2][1] = true;
	  break;
	}
    }
  }

  if(!compute_halfsphere[2][1])
    return 5;

  if(!compute_halfsphere[0][0]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(1,0,0),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hx()==0 && it->target().hx()==0) || it->is_long()) { 
	  compute_halfsphere[0][0] = true;
	  break;
	}
    }
  }
  
  if(!compute_halfsphere[0][0]) {
    compute_halfsphere[0][1] = true;
    return 0;
  }
  
  if(!compute_halfsphere[0][1]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(-1,0,0),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hx()==0 && it->target().hx()==0) || it->is_long()) { 
	  compute_halfsphere[0][1] = true;
	  break;
	}
    }
  }

  if(!compute_halfsphere[0][1])
    return 1;
  

  if(!compute_halfsphere[1][0]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(0,1,0),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hy()==0 && it->target().hy()==0) || it->is_long()) { 
	  compute_halfsphere[1][0] = true;
	  break;
	}
    }
  }
  
  if(!compute_halfsphere[1][0]) {
    compute_halfsphere[1][1] = true;
    return 2;
  }
  
  if(!compute_halfsphere[1][1]) {
    CGAL_forall_iterators(it,L) {
      if(CGAL::orientation(it->source(),it->target(), Sphere_point(0,-1,0),
			   CGAL::ORIGIN+it->sphere_circle().orthogonal_vector())
	 == CGAL::NEGATIVE) 
	if((it->source().hy()==0 && it->target().hy()==0) || it->is_long()) { 
	  compute_halfsphere[1][1] = true;
	  break;
	}
    }
  }

  if(!compute_halfsphere[1][1])
    return 3;
  return -1;
}

template <typename Map>
void SM_overlayer<Map>::
create(const Sphere_circle& c)
{ SHalfloop_handle l1 = this->new_shalfloop_pair();
  SHalfloop_handle l2 = this->twin(l1);
  circle(l1) = c; circle(l2) = c.opposite();
  SFace_handle f1 = this->new_sface();
  SFace_handle f2 = this->new_sface();
  link_as_loop(l1,f1);
  link_as_loop(l2,f2);
}

template <typename Map>
void SM_overlayer<Map>::
subdivide(const Map* M0, const Map* M1)
{
  PI[0] = SM_const_decorator(M0); 
  PI[1] = SM_const_decorator(M1);
  Seg_list L;
  Seg_map  From;
  for (int i=0; i<2; ++i) {
    SVertex_const_iterator v;
    CGAL_forall_svertices(v,PI[i]) {
      TRACEN(v->point() << " from " << i << " mark " << v->mark());
      if ( !PI[i].is_isolated(v) ) continue;
      L.push_back(trivial_segment(PI[i],v));
      From[--L.end()] = Seg_info(v,i);
    }
    SHalfedge_const_iterator e;
    CGAL_forall_sedges(e,PI[i]) {
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
    if ( PI[i].has_shalfloop() ) {
      SHalfloop_const_handle shl = PI[i].shalfloop();
      Seg_pair p = two_segments(PI[i],shl);
      L.push_back(p.first); 
      L.push_back(p.second);
      From[--L.end()] = From[--(--L.end())] = 
        Seg_info(shl,i);
    }
  }

  CGAL_assertion_code(typename Seg_list::iterator it);
  CGAL_assertion_code(CGAL_forall_iterators(it,L) TRACEN("  "<<*it));

  bool compute_halfsphere[3][2];
  int cs = check_sphere(L, compute_halfsphere);

  Seg_list L_pos,L_neg;

  switch(cs) {
  case 1: 
    partition_to_halfsphere(L.begin(), L.end(), L_pos, From, 
			    Sphere_circle(1,0,0), Sphere_circle(0,0,-1),
			    compute_halfsphere[0][1]);
    break;
  case 0:
    partition_to_halfsphere(L.begin(), L.end(), L_neg, From, 
			    Sphere_circle(-1,0,0), Sphere_circle(0,0,-1),
			    compute_halfsphere[0][0]);
    break;
  case 3: 
    partition_to_halfsphere(L.begin(), L.end(), L_pos, From, 
			    Sphere_circle(0,1,0), Sphere_circle(1,0,0),
			    compute_halfsphere[1][1]);
    break;
  case 2:
    partition_to_halfsphere(L.begin(), L.end(), L_neg, From, 
			    Sphere_circle(0,-1,0), Sphere_circle(1,0,0), 
			    compute_halfsphere[1][0]);
    break;
  case 5: 
    partition_to_halfsphere(L.begin(), L.end(), L_pos, From, 
			    Sphere_circle(0,0,1), Sphere_circle(1,0,0),
			    compute_halfsphere[2][1]);
    break;
  case 4:
    partition_to_halfsphere(L.begin(), L.end(), L_neg, From, 
			    Sphere_circle(0,0,-1), Sphere_circle(1,0,0), 
			    compute_halfsphere[2][0]);
    break;
  case -1:
    partition_to_halfsphere(L.begin(), L.end(), L_pos, From, 
			    Sphere_circle(0,0,1), Sphere_circle(1,0,0), true);
    partition_to_halfsphere(L.begin(), L.end(), L_neg, From, 
			    Sphere_circle(0,0,-1), Sphere_circle(1,0,0), true);
    break;
  default: CGAL_assertion_msg(0, "wrong value");
  }

  cs = cs==-1 ? 2 : cs/2;

  typedef SMO_from_sm<Self,Seg_iterator,Seg_info> SM_output;
  typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
  typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
  
  typedef CGAL::Segment_overlay_traits< 
    Seg_list, SM_output, PH_geometry>  PHS_traits;
  typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;
    
  typedef CGAL::Segment_overlay_traits< 
    Seg_list, SM_output, NH_geometry> NHS_traits;
  typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  typedef typename PHS_traits::INPUT Input_range;

  SVertex_handle v;
  SHalfedge_handle e;
  SM_output O(*this,PI,From); 

  // DEBUG CODE: to do: have all svertices a halfedge below associated?
  TRACEN("Vertex info before swep");
  CGAL_assertion_code(SVertex_iterator svi);
  //  CGAL_assertion_code(
  //    for( svi=svertices_begin(); svi!=svertices_end(); svi++) {
  //      GenPtr i = info(svi);
  //      TRACEN("vertex "<<point(svi)<<" info "<<i<<
  //	     " marks "<<mark(svi,0)<<" "<<mark(svi,1));
  //    }
  //  )
  
  if(compute_halfsphere[cs][0]) {
    Positive_halfsphere_sweep SP(
	Input_range(L_pos.begin(),L_pos.end()),O,
	PH_geometry(cs));
    SP.sweep();
    v=--this->svertices_end(); e=--this->shalfedges_end();
  }
  
  if(compute_halfsphere[cs][1]) {
    Negative_halfsphere_sweep SM(
        Input_range(L_neg.begin(),L_neg.end()),O,
        NH_geometry(cs));
    SM.sweep();
  }

  if(compute_halfsphere[cs][0]) {
    ++v; 
    ++e;
  }
  else {
    v = this->svertices_begin(); 
    e = this->shalfedges_begin();
  }
  
  if(compute_halfsphere[cs][0])
    create_face_objects(this->shalfedges_begin(), e, this->svertices_begin(), v, O,
                        PH_geometry(cs));
  if(compute_halfsphere[cs][1])
    create_face_objects(e, this->shalfedges_end(), v, this->svertices_end(), O,
			NH_geometry(cs));

  SHalfedge_iterator u;
  CGAL_forall_shalfedges(u,*this) {
    Sphere_segment s(point(source(u)),point(target(u)));
    circle(u) = s.sphere_circle(); 
    circle(twin(u)) = s.sphere_circle().opposite();
    TRACEN(PH(u) << " with circle " << circle(u));
  }

  std::vector<Mark> mohs(4);
  SM_point_locator L0(M0);
  SM_point_locator L1(M1);
  
  L0.marks_of_halfspheres(mohs, 0, cs);
  L1.marks_of_halfspheres(mohs, 2, cs);

  if(compute_halfsphere[cs][0])
    complete_face_support(this->svertices_begin(), v, O, mohs, 0);
  if(compute_halfsphere[cs][1])
    complete_face_support(v, this->svertices_end(), O, mohs, 1);

  // DEBUG CODE: to do: have all svertices a halfedge below associated?
  TRACEN("Vertex info after swep");
  CGAL_assertion_code(
    for( svi=this->svertices_begin(); svi!=this->svertices_end(); svi++) {
      GenPtr i = info(svi);
      TRACEN("vertex "<<point(svi)<<" info "<<i<<
	     " marks "<<mark(svi,0)<<" "<<mark(svi,1));
    }
  )

  if(compute_halfsphere[cs][0] && compute_halfsphere[cs][1])
    merge_halfsphere_maps(this->svertices_begin(),v,O);
  else
    set_outer_face_mark(compute_halfsphere[cs][1], mohs);
  
  this->check_integrity_and_topological_planarity();

  TRACEN("subdivided");
  CGAL_assertion_code(CGAL_forall_svertices(v,*this) TRACEN(PH(v)));
}

template <typename Map>
void SM_overlayer<Map>::
set_outer_face_mark(int offset, const std::vector<Mark>& mohs) {

  SFace_handle sf = this->new_sface();
  mark(sf, 0) = mohs[offset];
  mark(sf, 1) = mohs[offset+2];

  SHalfedge_iterator e;
  CGAL_forall_shalfedges(e, *this) {
    if ( face(e) != SFace_handle() ) continue;
    link_as_face_cycle(e,sf); 
  }
  
  SVertex_handle v;
  CGAL_forall_svertices(v, *this) {
    if(!is_isolated(v) || face(v) != SFace_handle()) continue;
    link_as_isolated_vertex(v,sf);
  }
}

template <typename Map>
template <typename Iterator, typename T>
void SM_overlayer<Map>::
partition_to_halfsphere(Iterator start, Iterator beyond, Seg_list& L, 
			CGAL::Unique_hash_map<Iterator,T>& M, 
			Sphere_circle xycircle, Sphere_circle yzcircle, 
			bool include_equator) const
{ TRACEN("partition_to_halfsphere ");
//  CGAL_assertion(pos!=0);
  Sphere_segment s1,s2;
  //  Sphere_circle xycircle(0,0,pos);
  if(include_equator) {
    while ( start != beyond) {
      int i = start->intersection(xycircle,s1,s2);
      TRACEN("segment " << start->source() << " " << start->target());
      if (i>1) { 
	L.push_back(s2); M[--L.end()] = M[start];    
	TRACEN(">1 " << s2.source() << " " << s2.target()); 
      }
      if (i>0) { 
	L.push_back(s1); M[--L.end()] = M[start];    
	TRACEN(">0 " << s1.source() << " " << s1.target()); 
      }
      ++start;
    }
  }
  else {
    while(start != beyond) {
      L.push_back(*start);
      M[--L.end()] = M[start];
      ++start;
    }
  }

  // now all segments are split into hemispheres
  // we still have to:
  // - split segments containing our special poles y^-, y^+
  // - split halfcircles
  // - add four equator segments 

  //  Sphere_circle yzcircle(1,0,0);
  typename Seg_list::iterator it, itl;
  
  CGAL_forall_iterators(it,L) { TRACEN("  "<<*it);
    if ( equal_as_sets(it->sphere_circle(),xycircle) ) {
      TRACEN("  splitting xy seg "<<*it);
      bool added=false;
      int n1 =  it->intersection(yzcircle,s1,s2);
      if (n1 > 1 && !s2.is_degenerate()) { 
	M[ L.insert(it,s2) ] = M[it]; 
	added=true; 
	TRACEN(">1 " << s2.source() << " " << s2.target()); 
      }
      if (n1 > 0 && !s1.is_degenerate()) { 
	M[ L.insert(it,s1) ] = M[it]; 
	added = true; 
	TRACEN(">1 " << s1.source() << " " << s1.target()); 
      }
      int n2 =  it->intersection(yzcircle.opposite(),s1,s2);
      if (n2 > 1 && !s2.is_degenerate()) { 
	M[ L.insert(it,s2) ] = M[it]; 
	added=true; 
	TRACEN(">1 " << s2.source() << " " << s2.target()); 
      }
      if (n2 > 0 && !s1.is_degenerate()) { 
	M[ L.insert(it,s1) ] = M[it]; 
	added=true; 
	TRACEN(">1 " << s1.source() << " " << s1.target()); 
      }
      if(added) {
	itl = it; --it; L.erase(itl); 
	//	M[itl] = T();
      }
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

  if(include_equator) {
    // append 4 xy-equator segments:
    Sphere_point S(0,-1,0),N(0,1,0);
    Sphere_segment sp(S,N,xycircle);
    Sphere_segment sm(S,N,xycircle.opposite());
    Sphere_segment s[4];
    sp.split_halfcircle(s[0],s[1]);
    sm.split_halfcircle(s[2],s[3]);
    L.insert(L.end(),s,s+4);
  }
}

template <typename Map>
template <typename Below_accessor, typename Halfsphere_geometry>
void SM_overlayer<Map>::
create_face_objects(SHalfedge_iterator e_start, SHalfedge_iterator e_end,
  SVertex_iterator v_start, SVertex_iterator v_end,
  const Below_accessor& D, 
  const Halfsphere_geometry& SG)
{
  TRACEN("create_face_objects()");
  if(e_start != e_end) {
    CGAL::Unique_hash_map<SHalfedge_handle,int> SFaceCycle(-1);
    std::vector<SHalfedge_handle>  MinimalSHalfedge;
    SHalfedge_around_sface_circulator hfc(last_out_edge(v_start)),hend(hfc);
    TRACEN("equator cycle "<<PH(hfc));
    CGAL_For_all(hfc,hend) SFaceCycle[hfc]=0; // outer face cycle = 0
    MinimalSHalfedge.push_back(twin(first_out_edge(v_start)));
    int i=1; 
    for (SHalfedge_iterator e = e_start; e != e_end; ++e) {
      if ( SFaceCycle[e] >= 0 ) continue; // already assigned
      SHalfedge_around_sface_circulator hfc(e),hend(hfc);
      SHalfedge_handle e_min = e;
      TRACEN(""); 
      TRACEN("  face cycle numbering "<<i);
      CGAL_For_all(hfc,hend) {
	SFaceCycle[hfc]=i; // assign face cycle number
	if ( SG.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
	  e_min = hfc;      
	TRACEN(PH(hfc));
      } TRACEN("");
      MinimalSHalfedge.push_back(e_min);
      ++i;
    }
    
    for (int j=1; j<i; ++j) {
      SHalfedge_handle e = MinimalSHalfedge[j];
      TRACEN("  face cycle "<<j<<" minimal halfedge "<<PH(e));
      Sphere_point p1 = point(source(e)), 
	p2 = point(target(e)), 
	p3 = point(target(next(e)));
      if ( SG.orientation(p1,p2,p3) > 0 ) { // left_turn => outer face cycle
	SFace_handle f = this->new_sface();
	link_as_face_cycle(e,f);
	TRACEN("  creating new face object "<<&*f<<" bd "<<&*e);
      }
    }
    
    for (SHalfedge_iterator e = e_start; e != e_end; ++e) {
      if ( face(e) != SFace_handle() ) continue;
      if ( SFaceCycle[e] == 0 ) continue;
      TRACEN("linking hole "<<PH(e));
      SFace_handle f = determine_face(e,MinimalSHalfedge,SFaceCycle,D);
      if(f != SFace_handle())
	link_as_face_cycle(e,f);
    }
  }    

  for (SVertex_iterator v = v_start; v != v_end; ++v) {
    if ( !is_isolated(v) ) continue;
    SHalfedge_handle e_below = D.halfedge_below(v);
    CGAL_assertion( e_below != SHalfedge_handle() || e_start == e_end );
    if(e_below != SHalfedge_handle())
      link_as_isolated_vertex(v,face(e_below));
  }

}

template <typename Map>
template <typename Below_accessor>
void SM_overlayer<Map>::
complete_face_support(SVertex_iterator v_start, SVertex_iterator v_end,
  const Below_accessor& D, std::vector<Mark>& mohs, int offset) const
{ TRACEN("complete_face_support");
  for (SVertex_iterator v = v_start; v != v_end; ++v) { 
    TRACEN("VERTEX = "<<PH(v));
    Mark m_buffer[2];
    SHalfedge_handle e_below = halfedge_below(v);
    if ( v == v_start ) {     
      for (int i=0; i<2; ++i){
	m_buffer[i] = mohs[offset+2*i];
      } 
    } else if ( e_below != SHalfedge_handle() ) {
      for (int i=0; i<2; ++i) {
	TRACEN("edge below "<< PH(e_below) << " " << mark(e_below,i));
        m_buffer[i] = incident_mark(e_below,i); 
      }
    } else { // e_below does not exist
      //      CGAL_assertion( point(v).hz() == 0 && 
      //		   ( offset == 0 ? (point(v).hx() >= 0) : (point(v).hx()<=0)) );
      if(!is_isolated(v))
	for (int i=0; i<2; ++i) 
	  m_buffer[i] = incident_mark(previous(first_out_edge(v)),i);
    } TRACEN(" faces right and below "<<m_buffer[0]<<" "<<m_buffer[1]);

    for (int i=0; i<2; ++i) {
      Object_handle o = supp_object(v,i);
      if ( o == NULL ) { 
	TRACEN("no vertex support"); 
	mark(v,i) = m_buffer[i]; continue; 
      }
      SVertex_const_handle vs;
      SHalfedge_const_handle es;
      SHalfloop_const_handle ls;
      if ( assign(vs,o) ) { mark(v,i) = PI[i].mark(vs); continue; }
      if ( assign(es,supp_object(v,i)) ) {
        if ( point(source(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(source(es)); continue; }
        if ( point(target(es)) == point(v) ) 
        { mark(v,i) = PI[i].mark(target(es)); continue; }
        mark(v,i) = PI[i].mark(es); continue;
      }
      if ( assign(ls,o) ) { 
	mark(v,i) = PI[i].mark(ls); 
	TRACEN("loop " << PI[i].circle(ls)); continue; }
      CGAL_assertion_msg(0,"wrong handle");
    } TRACEN(" vertex marks "<<mark(v,0)<<" "<<mark(v,1));

    if ( is_isolated(v) ) continue;
    SHalfedge_around_svertex_circulator e(first_out_edge(v)), hend(e);
    CGAL_For_all(e,hend) {
      if ( !is_forward(e) ) break;
      TRACEN("  forward edge "<<PH(e));
      for (int i=0; i<2; ++i) {
        if ( supp_object(e,i) != NULL ) {
          SHalfedge_const_handle ei; 
          if ( assign(ei,supp_object(e,i)) ) { 
            if ( PI[i].circle(ei) != circle(e) ) { ei = PI[i].twin(ei); }
            CGAL_assertion( PI[i].circle(ei) == circle(e) ); 
            TRACEN("  supporting edge "<<i<<" "<<PH(ei));
            incident_mark(twin(e),i) =
              PI[i].mark(PI[i].face(PI[i].twin(ei)));
            mark(e,i) = mark(twin(e),i) = PI[i].mark(ei);
            incident_mark(e,i) = m_buffer[i] =
              PI[i].mark(PI[i].face(ei)); 
          }
          SHalfloop_const_handle li;
          if ( assign(li,supp_object(e,i)) ) { 
            if ( PI[i].circle(li) != circle(e) ) { li = PI[i].twin(li); }
            CGAL_assertion( PI[i].circle(li) == circle(e) ); 
            TRACEN("  supporting loop "<<i<<" "<<PH(li));
            incident_mark(twin(e),i) =
              PI[i].mark(PI[i].face(PI[i].twin(li)));
            mark(e,i) = mark(twin(e),i) = PI[i].mark(li);
            incident_mark(e,i) = m_buffer[i] =
              PI[i].mark(PI[i].face(li)); 
          }
        } else { TRACEN("  support from face below "<<i);
	  incident_mark(twin(e),i) = mark(e,i) = mark(twin(e),i) =
          incident_mark(e,i) = m_buffer[i];
        }
      } TRACEN("  face marks "<<m_buffer[0]<<" "<<m_buffer[1]);
    }

    TRACEN(" mark of "<<PH(v)<<" "<<mark(v,0)<<" "<<mark(v,1));
  }
 
  SFace_iterator f;
  for (f = this->sfaces_begin(); f != this->sfaces_end(); ++f) {
    assoc_info(f);
    Object_handle boundary_object = sface_cycles_begin(f);
    CGAL_assertion(boundary_object != NULL);
    SHalfedge_handle e;
    if ( !CGAL::assign(e,boundary_object) ) 
      CGAL_assertion_msg(0,"Outer face cycle should be first.");
    for (int i=0; i<2; ++i) mark(f,i) = incident_mark(e,i);
  }
}

template <typename Map>
template <typename Mark_accessor>
void SM_overlayer<Map>::
merge_nodes(SHalfedge_handle e1, SHalfedge_handle e2,
  const Mark_accessor& D)
{
  SVertex_handle v1 = source(e1), v2 = target(e2);
  TRACEN("merge_nodes "<<PH(v1)<<PH(v2));
  CGAL_assertion(point(v1)==point(v2));
  SHalfedge_handle ep1 = previous(e1), en2 = next(e2);
  SHalfedge_around_svertex_circulator eav(out_edges(v2)),ee(eav);
  CGAL_For_all(eav,ee) { set_source(eav,v1); }
  link_as_prev_next_pair(e2,e1);  
  link_as_prev_next_pair(ep1,en2); 
  D.assert_equal_marks(v1,v2);
  D.discard_info(v2);
  delete_vertex_only(v2);
}

template <typename Map>
template <typename Mark_accessor>
void SM_overlayer<Map>::
merge_halfsphere_maps(SVertex_handle v1, SVertex_handle v2,
  const Mark_accessor& D)
{ TRACEN("merging halfspheres "<<PH(v1)<<PH(v2));
  CGAL_assertion(point(v1)==point(v2));
  std::list<SHalfedge_pair> L_equator;
  SHalfedge_around_sface_circulator 
    ep(last_out_edge(v1)), en(twin(first_out_edge(v2)));
  do { 
   L_equator.push_back(SHalfedge_pair(ep,en));
   merge_nodes(ep,en,D); ++ep; --en; 
  } while ( source(ep) != v1 );
  
  typename std::list<SHalfedge_pair>::iterator it;
  CGAL_forall_iterators(it,L_equator) { 
    SHalfedge_handle e1 = it->first, e2 = it->second;
    SHalfedge_handle e1t = twin(e1), e2t = twin(e2);
    TRACEV(PH(e1));TRACEV(PH(e2));
    SHalfedge_handle e2tp = previous(e2t);
    SHalfedge_handle e2tn = next(e2t);
    link_as_prev_next_pair(e2tp,e1);
    link_as_prev_next_pair(e1,e2tn);
    SFace_handle f = face(e2t);
    if ( is_sm_boundary_object(e2t) )
    { undo_sm_boundary_object(e2t,f); store_sm_boundary_object(e1,f); }
    set_face(e1,f);
    if ( e2 == first_out_edge(source(e2)) )
      set_first_out_edge(source(e2),e1t);
    D.discard_info(e2);
    delete_edge_pair_only(e2);
  }
}

template <typename Map>
template <typename Selection>
void SM_overlayer<Map>::
select(const Selection& SP) const
{ 
  SVertex_iterator v;
  CGAL_forall_svertices(v,*this) {
    mark(v) = SP(mark(v,0),mark(v,1));
    discard_info(v); 
  }
  SHalfedge_iterator e;
  CGAL_forall_sedges(e,*this) {
    mark(e) = SP(mark(e,0),mark(e,1));
    mark(twin(e)) = SP(mark(twin(e),0),mark(twin(e),1));
    CGAL_assertion(mark(e) == mark(twin(e)));
    discard_info(e);
  }
  SFace_iterator f;
  CGAL_forall_sfaces(f,*this) {
    mark(f) = SP(mark(f,0),mark(f,1));
    discard_info(f);
  }
}

template <typename Map>
void SM_overlayer<Map>::simplify()
{
  TRACEN("simplifying"); 

  typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
  CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem(NULL);
  CGAL::Unique_hash_map< SVertex_handle, Union_find_handle> Vitem(NULL);
  CGAL::Union_find< SFace_handle> UF;
  
  SFace_iterator f;
  CGAL_forall_sfaces(f,*this) {
     Pitem[f] = UF.make_set(f);
     clear_face_cycle_entries(f);
  }

  if ( this->has_shalfloop() ) {
    SHalfloop_handle l = this->shalfloop();
    SFace_handle f = *(UF.find(Pitem[face(l)]));
    link_as_loop(l,f);
    f = *(UF.find(Pitem[face(twin(l))]));
    link_as_loop(twin(l),f);
  }

  SHalfedge_iterator e, en;
  for(e = this->shalfedges_begin(); e != this->shalfedges_end(); e = en) { 
    en = e; ++en; if ( en==twin(e) ) ++en;
    TRACEN("can simplify ? " << PH(e));
    TRACEN(mark(e) << " " << mark(face(e)) << " " << mark(face(twin(e))));
    if (( mark(e) == mark(face(e)) && mark(e) == mark(face(twin(e))))){
      TRACEN("deleting "<<PH(e));
      if ( !UF.same_set(Pitem[face(e)],
			Pitem[face(twin(e))]) ) {
	
	UF.unify_sets( Pitem[face(e)],
		       Pitem[face(twin(e))] );
	TRACEN("unioning disjoint faces");
      }
      
      TRACEN("is_closed_at_source " << is_closed_at_source(e) << 
	     " " << is_closed_at_source(twin(e)));
      
      if ( is_closed_at_source(e) )
	Vitem[source(e)] = Pitem[face(e)];
      
      if ( is_closed_at_source(twin(e)))
	Vitem[target(e)] = Pitem[face(e)];
      
      delete_edge_pair(e);
    }
  }

  CGAL::Unique_hash_map<SHalfedge_handle,bool> linked(false);
  for (e = this->shalfedges_begin(); e != this->shalfedges_end(); ++e) {
    if ( linked[e] ) continue;
    SHalfedge_around_sface_circulator hfc(e),hend(hfc);
    SFace_handle f = *(UF.find( Pitem[face(e)]));
    CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
    store_sm_boundary_object(e,f);
  }

  SVertex_iterator v,vn;
  for(v = this->svertices_begin(); v != this->svertices_end(); v=vn) {
    vn=v; ++vn;
    if ( is_isolated(v) ) {
    
      if(Vitem[v] != NULL) {
	set_face(v,*(UF.find(Vitem[v])));
	TRACEN("incident face of " << PH(v) << " set to " << &*(face(v)));
      }
      else {
	set_face(v, *(UF.find(Pitem[face(v)])));
	TRACEN("isolated svertex " << PH(v) << 
	       " already has incident face " << &*(face(v)));
      }

      if ( mark(v) == mark(face(v)) ) {
        TRACEN("removing isolated vertex"<<PH(v));
        delete_vertex_only(v);  
      } 
      else 
        store_sm_boundary_object(v,face(v)); // isolated, but should stay
    } else { // v not isolated
      SHalfedge_handle e2 = first_out_edge(v), e1 = previous(e2);
      if ( has_outdeg_two(v) &&
           mark(v) == mark(e1) && mark(v) == mark(e2) &&
           circle(e1) == circle(e2) ) {
        TRACEN("collinear at "<<PH(v)<<PH(e1)<<PH(e2));
        if ( e1 == e2 ){ TRACEN("edge_to_loop"); convert_edge_to_loop(e1);}
        else {TRACEN("merge_edge_pairs"); merge_edge_pairs_at_target(e1); } 
      }
    }
  }

  SFace_iterator fn;
  for (f = fn = this->sfaces_begin(); f != this->sfaces_end(); f=fn) { 
    ++fn;
    Union_find_handle pit = Pitem[f];
    if ( UF.find(pit) != pit ) {
      TRACEN("delete face " << &*f);
      delete_face_only(f);
    }
  }
}

template <typename Map>
template <typename Iterator>
void SM_overlayer<Map>::
subdivide_segments(Iterator start, Iterator end) const
{
typedef SMO_decorator<SM_decorator,Iterator> SM_output;
typedef typename Sphere_kernel::Positive_halfsphere_geometry PH_geometry;
typedef CGAL::Segment_overlay_traits< 
          Iterator, SM_output, PH_geometry>  PHS_traits;
typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

typedef typename Sphere_kernel::Negative_halfsphere_geometry NH_geometry;
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

  SVertex_handle v1,v2;
  SM_output O(*this); 
  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(Input_range(Lp.begin(),Lp.end()),O);
  SP.sweep();
  //TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v1= this->vertices_begin(); v2=--this->vertices_end();
  Negative_halfsphere_sweep SM(Input_range(Lm.begin(),Lm.end()),O);
  SM.sweep();
  //TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  ++v2;
  // now two CCs of sphere graph calculated
  // v1 = first node of CC in positive xy-sphere
  // v2 = first node of CC in negative xy-sphere

  merge_halfsphere_maps(v1,v2,O);
  this->check_integrity_and_topological_planarity(false);
}


CGAL_END_NAMESPACE
#endif //CGAL_SM_OVERLAYER_H


