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
// file          : include/CGAL/Nef_2/PM_overlayer.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/PM_overlayer.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Overlay module for plane maps
// ============================================================================

#ifndef CGAL_PM_OVERLAYER_H
#define CGAL_PM_OVERLAYER_H

#include <CGAL/basic.h>
#include <CGAL/Hash_map.h>
#include <CGAL/Partition.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#undef _DEBUG
#define _DEBUG 13
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

CGAL_BEGIN_NAMESPACE

/*{\Manpage {PM_overlayer}{PMDEC,GEOM}{A HDS overlay calculator}{O}}*/
template <typename PMDEC, typename GEOM>
class PM_overlayer : public PMDEC {
  typedef PM_overlayer<PMDEC,GEOM>  Self;
  typedef PMDEC                     Base;
  const GEOM& K; // geometry reference

/*{\Mdefinition An instance |\Mvar| of data type |\Mname| is a
decorator object offering plane map overlay calculation. Overlay is
either calculated from two plane maps or from a set of segments.  The
result is stored in a plane map |P| that carries the geometry and the
topology of the overlay.

The two template parameters allow to adapt the overlay calculation to
different scenarios.  The template parameter |PMDEC| has to be a model
conforming to our plane map decorator concept |PM_decorator|.  The
concept also describes the interface how the topological information
stored in |P| can be extracted.  The geometry |GEOM| has to be a model
conforming to our two dimensional geometry kernel concept
|Affine_geometry|.

The overlay of a set of segments $S$ is stored in a plane map $P =
(V,E,F)$. Vertices are either the endpoints of segments (trivial
segments are allowed) or the result of the internal intersection of
two segments. Between two vertices there's an edge if there's 
a segment that supports the straight line embedding of $e$ and
if there's no vertex in the relative interior of the embedding of 
$e$.

The faces refer to the maximal connected open point sets of the
planar subdivision implied by the embedding of the vertices and edges.
Faces are bounded by possibly several face cycles\footnote{For the
definition of plane maps and their concepts see the manual page of
|PM_decorator|.} including isolated vertices. The overlay process in
the method |create| creates the objects, the topology of the result
and allows to link the plane map objects to input segments by means of
a data accessor. The method starts from zero- and one-dimensional
geometric objects in $S$ and produces a planar structure where each
point of the plane can be assigned to an object (vertex, edge, or
face) of |P|.

The overlay of two plane maps $P_i = (V_i, E_i, F_i)$ has the
additional aspect that we already start from two planar subdivisions.
We use the index $i=0,1$ defining the reference to $P_i$, unindexed
variables refer to the resulting plane map $P$.  The $1$-skeleta of
the two maps subdivide the edges and faces of the complementary
structure into smaller units. This means vertices and edges of $P_i$
can split edges of $P_{1-i}$ and face cycles of $P_i$ subdivide faces
of $P_{1-i}$. The 1-skeleton $G$ of $P$ is defined by the overlay of
the embedding of the 1-skeleta of $P_0$ and $P_1$ (Take a trivial
segment for each vertex and a segment for each edge and use the
overlay definition of a set of segments above). The faces of $P$ refer
to the maximal connected open point sets of the planar subdivision
implied by the embedding of $G$. Each object from the output tuple
$(V,E,F)$ has a \emph{supporting} object $u_i$ in each of the two
input structures.  Imagine the two maps to be transparencies, which we
stack. Then each point of the plane is covered by an object from each
of the input structures.  This support relationship from the input
structures to the output structure defines an information flow. Each
supporting object $u_i$ of $u$ $(i=0,1)$ carries an associated
information $|mark|(u_i)$. After the subdivision operation this
information is attributed to the output object $u$ by
$|mark|(u,i)$.}*/

/*{\Mgeneralization PM_decorator}*/

public:
/*{\Mtypes 3}*/
  typedef PMDEC                     PM_decorator;
  /*{\Mtypemember the plane map decorator |PMDEC|.}*/
  typedef typename PMDEC::Plane_map Plane_map;
  /*{\Mtypemember the plane map type decorated by |PMDEC|.}*/
  typedef GEOM                      Geometry;
  /*{\Mtypemember the geometry kernel |GEOM|.}*/
  typedef typename GEOM::Point_2    Point;
  /*{\Mtypemember the point type of the geometric kernel, 
     \precond |Point_2| equals |Plane_map::Point_2|.}*/
  typedef typename GEOM::Segment_2  Segment;
  /*{\Mtypemember the segment type of the geometric kernel.}*/
  typedef typename PMDEC::Mark      Mark;
  /*{\Mtypemember the mark of plane map objects.}*/

  #define USING(t) typedef typename PMDEC::t t
  typedef typename Base::Base PM_const_decorator;
  USING(Halfedge_handle);
  USING(Vertex_handle);
  USING(Face_handle);
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Face_iterator);
  USING(Halfedge_const_handle);
  USING(Vertex_const_handle);
  USING(Face_const_handle);
  USING(Halfedge_const_iterator);
  USING(Vertex_const_iterator);
  USING(Face_const_iterator);
  USING(Halfedge_around_vertex_circulator);
  USING(Halfedge_around_face_circulator);
  USING(Hole_iterator);
  USING(Isolated_vertex_iterator);
  #undef USING

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

  struct Seg_info { // to transport information from input to output
    Halfedge_const_handle _e;
    Vertex_const_handle   _v;
    int                   _i;

    Seg_info() : _i(-1) {}
    Seg_info(Halfedge_const_handle e, int i) 
    { _e=e; _i=i; }
    Seg_info(Vertex_const_handle v, int i) 
    { _v=v; _i=i; }
    Seg_info(const Seg_info& si) 
    { _e=si._e; _v=si._v; _i=si._i; }
    Seg_info& operator=(const Seg_info& si) 
    { _e=si._e; _v=si._v; _i=si._i; return *this; }
    LEDA_MEMORY(Seg_info)
  };

  typedef std::list<Segment>                   Seg_list;
  typedef typename Seg_list::const_iterator    Seg_iterator;
  typedef std::pair<Seg_iterator,Seg_iterator> Seg_it_pair;


/*{\Mcreation 6}*/
PM_overlayer(Plane_map& P, const Geometry& g = Geometry()) : 
/*{\Mcreate |\Mvar| is a decorator object manipulating |P|.}*/
  Base(P), K(g) {}

/*{\Moperations 1.1 1}*/

template <typename Forward_iterator, typename VE_data_accessor>
void create(Forward_iterator start, Forward_iterator end, 
            VE_data_accessor& A) const;
/*{\Mop produces in |P| the plane map consistent with the overlay
of the segments from the iterator range |[start,end)|. The data accessor 
|A| allows to initialize created vertices and edges with respect to the
segments in the iterator range. |A| requires the following methods:\\
[[void supporting_segment(Halfedge_handle e, Forward_iterator it)]]\\
[[void trivial_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void starting_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void passing_segment(Vertex_handle v, Forward_iterator it)]]\\
[[void ending_segment(Vertex_handle v, Forward_iterator it)]]\\
Where |supporting_segment| is called for each segment |*it|
supporting a newly created edge |e|, |trivial_segment| is called for
each trivial segment |*it| supporting a newly created vertex |v|, and
the three last operations are called for each non-trivial segment
|*it| starting at/passing through/ending at the embedding of a newly
created vertex |v|.
\precond |Forward_iterator| has value type |Segment|.}*/

void subdivide(const Plane_map& P0, const Plane_map& P1) const;
/*{\Mop constructs the overlay of the plane maps |P0| and |P1| in
|P|, where all objects (vertices, halfedges, faces) of |P| are
\emph{enriched} by the marks of the supporting objects of the two
input structures: e.g. let |v| be a vertex supported by a node |v0| in
|P0| and by a face |f1| in |P1| and |D0|, |D1| be decorators of
type |PM_decorator| on |P0|,|P1|. Then |\Mvar.mark(v,0) = D0.mark(v0)|
and |\Mvar.mark(v,1) = D1.mark(f1)|.}*/

template <typename SELECTION> 
void select(SELECTION& SP) const;
/*{\Mop sets the marks of all objects according to the selection
predicate |SP|. |SELECTION| has to be a function object type with a
function operator\\
[[Mark operator()(Mark m0, Mark m1)]]\\
For each object |u| of |P| enriched by the marks of the supporting
objects according to the previous procedure |subdivide|, after this
operation |\Mvar.mark(u) = SP ( \Mvar.mark(u,0),\Mvar.mark(u,1)
)|. The additional marks are invalidated afterwards. }*/

void simplify() const;
/*{\Mop simplifies the structure of |P| according to the marks of
its objects. An edge |e| separating two faces |f1| and |f2| and equal
marks |mark(e) == mark(f1) == mark(f2)| is removed and the faces are
unified.  An isolated vertex |v| in a face |f| with |mark(v)==mark(f)|
is removed.  A vertex |v| with outdegree two, two collinear out-edges
|e1|,|e2| and equal marks |mark(v) == mark(e1) == mark(e2)| is removed
and the edges are unified.}*/

bool is_outer_face_cycle_edge(Halfedge_handle h) const
/*{\Xop returns |true| when |h| or |twin(h)| is part of the outer 
face cycle.}*/
{ return ( face(h) == faces_begin() ||
           face(twin(h)) == faces_begin() ); }


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
{ geninfo<vertex_info>::create(info(v)); }

void discard_info(Vertex_handle v) const
{ geninfo<vertex_info>::clear(info(v)); }

vertex_info& ginfo(Vertex_handle v) const
{ return geninfo<vertex_info>::access(info(v)); }

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
{ geninfo<halfedge_info>::create(info(e)); 
  geninfo<halfedge_info>::create(info(twin(e))); }

void discard_info(Halfedge_handle e)  const
{ geninfo<halfedge_info>::clear(info(e)); 
  geninfo<halfedge_info>::clear(info(twin(e))); }

halfedge_info& ginfo(Halfedge_handle e)  const
{ return geninfo<halfedge_info>::access(info(e)); }

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
{ geninfo<face_info>::create(info(f)); }

void discard_info(Face_handle f)  const
{ geninfo<face_info>::clear(info(f)); }

face_info& ginfo(Face_handle f)  const
{ return geninfo<face_info>::access(info(f)); }

Mark& mark(Face_handle f, int i)  const
{ return ginfo(f).m[i]; }

void clear_associated_info_of_all_objects() const 
{
  Vertex_iterator vit;
  for (vit = vertices_begin(); vit != vertices_end(); ++vit)
    discard_info(vit);
  Halfedge_iterator hit;
  for (hit = halfedges_begin(); hit != halfedges_end(); ++hit) 
    discard_info(hit);
  Face_iterator fit;
  for (fit = faces_begin(); fit != faces_end(); ++fit) 
    discard_info(fit);
}

template <class Below_info>
void create_face_objects(const Below_info& D) const
{
  TRACEN("create_face_objects()");
  CGAL::Hash_map<Halfedge_handle,int> FaceCycle(-1);
  std::vector<Halfedge_handle>  MinimalHalfedge;
  int i=0;
  Halfedge_iterator e, eend = halfedges_end();
  for (e=halfedges_begin(); e != eend; ++e) {
    if ( FaceCycle[e] >= 0 ) continue; // already assigned
    Halfedge_around_face_circulator hfc(e),hend(hfc);
    Halfedge_handle e_min = e;
    TRACE("face cycle "<<i<<"\n");
    CGAL_For_all(hfc,hend) {
      FaceCycle[hfc]=i; // assign face cycle number
      if ( K.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      TRACE(PE(hfc));
    } 
    TRACEN("");
    MinimalHalfedge.push_back(e_min);
    ++i;
  }

  Face_handle f_outer = new_face(); assoc_info(f_outer);
  for (int j=0; j<i; ++j) {
    Halfedge_handle e = MinimalHalfedge[j];
      TRACEN("  face cycle "<<j);TRACEN("  minimal halfedge "<<PE(e));
    Point p1 = point(source(e)), 
          p2 = point(target(e)), 
          p3 = point(target(next(e)));
    if ( K.leftturn(p1,p2,p3) ) { // leftturn => outer face cycle
        TRACEN("  creating new face object");
      Face_handle f = new_face();
      link_as_outer_face_cycle(f,e);
    }
  }

  for (e = halfedges_begin(); e != eend; ++e) {
    if ( face(e) != Face_handle() ) continue;
    TRACEN("linking hole "<<PE(e));
    Face_handle f = determine_face(e,MinimalHalfedge,FaceCycle,D);
    link_as_hole(f,e);
  }
  Vertex_iterator v, v_end = vertices_end();
  for (v = vertices_begin(); v != v_end; ++v) {
    if ( !is_isolated(v) ) continue;
    Halfedge_handle e_below = D.halfedge_below(v);
    if ( e_below == Halfedge_handle() ) 
      link_as_isolated_vertex(f_outer,v);
    else
      link_as_isolated_vertex(face(e_below),v);    
  }

}

template <class Below_info>
Face_handle determine_face(Halfedge_handle e, 
  const std::vector<Halfedge_handle>& MinimalHalfedge,
  const CGAL::Hash_map<Halfedge_handle,int>& FaceCycle,
  const Below_info& D) const
{ TRACEN("determine_face "<<PE(e));
  int fc = FaceCycle[e];
  Halfedge_handle e_min = MinimalHalfedge[fc];
  Halfedge_handle e_below = D.halfedge_below(target(e_min));
  if ( e_below == Halfedge_handle() ) { // below is nirwana
    return faces_begin();
  }
  Face_handle f = face(e_below);
  if (f != Face_handle()) return f; // has face already
  f = determine_face(e_below, MinimalHalfedge, FaceCycle,D);
  link_as_hole(f,e_below);
  return f;
}

Segment segment(const PM_const_decorator& N, 
                Halfedge_const_handle e) const
{ return K.construct_segment(
    N.point(N.source(e)),N.point(N.target(e))); }

Segment segment(const PM_const_decorator& N, 
                Vertex_const_handle v) const
{ Point p = N.point(v); 
  return K.construct_segment(p,p); }

bool is_forward_edge(const PM_const_decorator& N, 
                     Halfedge_const_iterator hit) const
{ Point p1 = N.point(N.source(hit));
  Point p2 = N.point(N.target(hit));
  return (K.compare_xy(p1,p2) < 0); }

void assert_type_precondition() const
{ typename PMDEC::Point p1;
  Point p2;
  assert_equal_types(p1,p2); }



}; // PM_overlayer<PMDEC,GEOM>


template <typename PMDEC, typename I, typename DA>
struct PMO_from_segs {
  typedef PMDEC PM_decorator;
  typedef typename PMDEC::Vertex_handle   Vertex_handle;
  typedef typename PMDEC::Halfedge_handle Halfedge_handle;
  typedef typename PMDEC::Point           Point;
  const PMDEC& G;
  DA& D; 
  PMO_from_segs(const PMDEC& Gi, DA& Di) : 
    G(Gi),D(Di) {}

  Vertex_handle new_vertex(const Point& p)
  { Vertex_handle v = G.new_vertex(p); 
    geninfo<Halfedge_handle>::create(G.info(v));
    return v;
  }

  void link_as_target_and_append(Vertex_handle v, Halfedge_handle e) 
  { G.link_as_target_and_append(v,e); }

  Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
  { Halfedge_handle e = 
    G.new_halfedge_pair_at_source(v,PM_decorator::BEFORE); 
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
  { geninfo<Halfedge_handle>::access(G.info(v)) = e; }

  Halfedge_handle halfedge_below(Vertex_handle v) const
  { return geninfo<Halfedge_handle>::access(G.info(v)); }

  void clear_temporary_vertex_info() const
  { Vertex_handle v;
    for(v = G.vertices_begin(); v!= G.vertices_end(); ++v)
      geninfo<Halfedge_handle>::clear(G.info(v));
  }


}; // PMO_from_segs


template <typename PMDEC, typename IT, typename INFO>
struct PMO_from_pm {
  typedef PMDEC PM_decorator;
  typedef typename PMDEC::PM_const_decorator PM_const_decorator;
  typedef typename PM_decorator::Vertex_handle Vertex_handle;
  typedef typename PM_decorator::Halfedge_handle Halfedge_handle;
  typedef typename PM_decorator::Vertex_const_handle Vertex_const_handle;
  typedef typename PM_decorator::Halfedge_const_handle Halfedge_const_handle;
  typedef typename PM_decorator::Point Point;

  const PM_decorator& G;
  const PM_const_decorator* pGI[2];
  CGAL::Hash_map<IT,INFO>& M;
  PMO_from_pm(const PM_decorator& Gi, 
              const PM_const_decorator* pG0, 
              const PM_const_decorator* pG1,
              CGAL::Hash_map<IT,INFO>& Mi) : G(Gi),M(Mi) 
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
  G.new_halfedge_pair_at_source(v,PM_decorator::BEFORE); 
  G.assoc_info(e);
  return e;
}

void halfedge_below(Vertex_handle v, Halfedge_handle e) const
{ G.halfedge_below(v) = e; }

void supporting_segment(Halfedge_handle e, IT it) const
{ INFO& si = M[it];
  assert( si._e != Halfedge_const_handle() );
  G.supp_halfedge(e,si._i) = si._e;
  G.is_forward(e) = true;
  TRACEN("   supporting "<<si._i<<" "<<*it);
}


void trivial_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  assert( si._v != Vertex_const_handle() );
  G.supp_vertex(v,si._i) = si._v; 
}

void starting_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  G.supp_vertex(v,si._i) = pGI[si._i]->source(si._e);
}

void ending_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  G.supp_vertex(v,si._i) = pGI[si._i]->target(si._e);
}

void passing_segment(Vertex_handle v, IT it) const
{ INFO& si = M[it];
  G.supp_halfedge(v,si._i) = si._e; 
}

Halfedge_handle halfedge_below(Vertex_handle v) const
{ return G.halfedge_below(v); }


}; // PMO_from_pm

template <typename PMDEC, typename GEOM>
template <typename Forward_iterator, typename Object_data_accessor>
void PM_overlayer<PMDEC,GEOM>::
create(Forward_iterator start, Forward_iterator end, 
       Object_data_accessor& Obj) const
{
  TRACEN("creating from iterator range");
  typedef PMO_from_segs<Self,Forward_iterator,Object_data_accessor> PMO;
  typedef Segment_overlay_traits<Forward_iterator,PMO,GEOM> seg_overlay;
  typedef gen_plane_sweep< seg_overlay > seg_overlay_sweep;
  PMO Out(*this, Obj);
  seg_overlay_sweep SOS( seg_overlay::INPUT(start, end), Out, K);
  SOS.sweep();
  create_face_objects(Out);
  Out.clear_temporary_vertex_info();
}


template <typename PMDEC, typename GEOM>
void PM_overlayer<PMDEC,GEOM>::
subdivide(const Plane_map& P0, const Plane_map& P1) const
{
  PM_const_decorator PI[2];
  PI[0] = PM_const_decorator(P0); PI[1] = PM_const_decorator(P1);
  Seg_list Segments;
  CGAL::Hash_map<Seg_iterator,Seg_info> From;
  for (int i=0; i<2; ++i) {
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


  typedef PMO_from_pm<Self,Seg_iterator,Seg_info> PMO;
  typedef Segment_overlay_traits<Seg_iterator,PMO,GEOM> pm_overlay;
  typedef gen_plane_sweep< pm_overlay > pm_overlay_sweep;
  PMO Out(*this,&PI[0],&PI[1],From);
  pm_overlay_sweep SOS(Seg_it_pair(Segments.begin(),Segments.end()),Out,K);
  SOS.sweep();
  create_face_objects(Out);


  TRACEN("creating marks");
  Vertex_iterator v, vend = vertices_end();
  for (v = vertices_begin(); v != vend; ++v) {
    TRACEN("mark at "<<PV(v));
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

    for (int i=0; i<2; ++i) 
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
        TRACEN("   halfedge "<<PE(e));
        Halfedge_const_handle ei;
        bool supported;
        for (int i=0; i<2; ++i) {
          supported = ( supp_halfedge(e,i) != Halfedge_const_handle() );
          if ( supported ) {
            ei = supp_halfedge(e,i);
            TRACEN("   supp halfedge "<<i<<" "<<PE(ei));
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
  Face_iterator f = faces_begin();
  for (int i=0; i<2; ++i) mark(f,i) = PI[i].mark(PI[i].faces_begin());
  for (++f; f != faces_end(); ++f) { // skip first face
    assoc_info(f);
    for (int i=0; i<2; ++i) mark(f,i) = incident_mark(halfedge(f),i);
  }


}


template <typename PMDEC, typename GEOM>
template <typename SELECTION>
void PM_overlayer<PMDEC,GEOM>::
select(SELECTION& predicate) const
{ 
  Vertex_iterator vit = vertices_begin(),
                  vend = vertices_end();
  for( ; vit != vend; ++vit) {
    mark(vit) = predicate(mark(vit,0),mark(vit,1));
    discard_info(vit); 
  }
  Halfedge_iterator hit = halfedges_begin(),
                    hend = halfedges_end();
  for(; hit != hend; ++(++hit)) {
    mark(hit) = predicate(mark(hit,0),mark(hit,1));
    discard_info(hit);
  }
  Face_iterator fit = faces_begin(),
                fend = faces_end();
  for(; fit != fend; ++fit) {
    mark(fit) = predicate(mark(fit,0),mark(fit,1));
    discard_info(fit);
  }
}

template <typename PMDEC, typename GEOM>
void PM_overlayer<PMDEC,GEOM>::simplify() const
{
  TRACEN("simplifying"); 
  typedef CGAL::Partition<Face_handle>::item partition_item;
  CGAL::Hash_map<Face_iterator,partition_item> Pitem;
  CGAL::Partition<Face_handle> FP;

  Face_iterator fit, fend = faces_end();
  for (fit = faces_begin(); fit!= fend; ++fit) {
     Pitem[fit] = FP.make_block(fit);
     clear_face_cycle_entries(fit);
  }


  Halfedge_iterator hit = halfedges_begin(), hitn,
                    hend = halfedges_end();
  for(; hitn=hit, ++(++hitn), hit != hend; hit=hitn) { 
    if ( is_outer_face_cycle_edge(hit) ) continue;
    if ( mark(hit) == mark(face(hit)) &&
         mark(hit) == mark(face(twin(hit))) ) {
        TRACEN("deleting "<<PE(hit));
      if ( !FP.same_block(Pitem[face(hit)],
                          Pitem[face(twin(hit))]) ) {
        FP.union_blocks( Pitem[face(hit)],
                         Pitem[face(twin(hit))] );
        TRACEN("unioning disjoint faces");
      }
      if ( is_closed_at_source(hit) ) 
        set_face(source(hit),face(hit));
      if ( is_closed_at_source(twin(hit)) ) 
        set_face(target(hit),face(hit));
      delete_halfedge_pair(hit);
    }
  }

  CGAL::Hash_map<Halfedge_handle,bool> linked(false);
  Halfedge_iterator eit, eend = halfedges_end();
  for (eit = halfedges_begin(); eit != eend; ++eit) {
    if ( linked[eit] ) continue;
    Halfedge_around_face_circulator hfc(eit),hend(hfc);
    Halfedge_handle e_min = eit;
    Face_handle f = FP.inf(FP.find(Pitem[face(eit)]));
    CGAL_For_all(hfc,hend) {
      set_face(hfc,f);
      if ( K.compare_xy(point(target(hfc)), point(target(e_min))) < 0 )
        e_min = hfc;
      linked[hfc]=true;
    }
    Point p1 = point(source(e_min)),
          p2 = point(target(e_min)),
          p3 = point(target(next(e_min)));
    if ( K.orientation(p1,p2,p3) > 0 ) 
      // leftturn => outer face cycle 
      set_halfedge(f,e_min); // store as outer
    else // inner hole cycle
      set_hole(f,e_min); // store as inner
  }


  Vertex_iterator vit,vitn, vend = vertices_end();
  for(vit = vertices_begin(); vit != vend; vit=vitn) {
    vitn=vit; ++vitn;
    if ( is_isolated(vit) ) {
      if ( mark(vit) == mark(face(vit)) )
        delete_vertex_only(vit);  // isolated, to be removed
      else
        set_isolated_vertex(face(vit),vit); // isolated, but should stay
    } else { // vit not isolated
      Halfedge_handle e2 = first_out_edge(vit), e1 = previous(e2);
      Point p1 = point(source(e1)),
            p2 = point(vit),
            p3 = point(target(e2));
      if ( has_outdeg_two(vit) &&
           mark(vit) == mark(e1) && mark(vit) == mark(e2) &&
           (K.orientation(p1,p2,p3) == 0) ) 
        merge_halfedge_pairs_at_target(e1); 
    }
  }

  Face_iterator fitn;
  for (fit = faces_begin(); fit != fend; fit=fitn) {
    fitn=fit; ++fitn;
    partition_item pit = Pitem[fit];
    if ( FP.find(pit) != pit ) delete_face(fit);
  }


}


CGAL_END_NAMESPACE
#endif // CGAL_PM_OVERLAYER_H

