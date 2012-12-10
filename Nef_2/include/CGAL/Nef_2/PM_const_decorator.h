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

#ifndef CGAL_PM_CONST_DECORATOR_H
#define CGAL_PM_CONST_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <string>
#include <list>
#include <sstream>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_2/iterator_tools.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 7
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

template <typename Iter, typename Move>
inline CGAL::Circulator_tag  
query_circulator_or_iterator(const CircFromIt<Iter,Move>& )
{ return CGAL::Circulator_tag(); }


template <typename HE>
class move_halfedge_around_vertex {
public:
  void forward(HE& e) const  { e = (e->prev()->opposite()); }
  void backward(HE& e) const { e = (e->opposite()->next()); }
};

template <typename HE>
struct move_halfedge_around_face {
  void forward(HE& e)  const { e = (e->next()); }
  void backward(HE& e) const { e = (e->prev()); }
};


/*{\Moptions print_title=yes}*/ 
/*{\Moptions constref=yes}*/ 
/*{\Msubst 
PM_const_decorator#PMConstDecorator
PM_decorator#PMDecorator
}*/ 
/*{\Moptions outfile=PMConstDecorator.man }*/
/*{\Manpage {PMConstDecorator}{}{Topological plane map exploration}{D}}*/

template <typename HDS> class PM_decorator;

template <typename HDS>
class PM_const_decorator 
{ typedef PM_const_decorator<HDS> Self;
/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a 
decorator for interfacing the topological structure of a plane map |P| 
(read-only).

A plane map |P| consists of a triple $(V, E, F)$ of vertices, edges, and
faces. We collectively call them objects. An edge |e| is a pair of
vertices |(v,w)| with incidence operations |v = source(e)|, |w =
target(e)|. The list of all edges with source |v| is called the
adjacency list |A(v)|.

Edges are paired into twins. For each edge |e = (v,w)| there's an edge
|twin(e) = (w,v)| and |twin(twin(e)) == e|\cgalFootnote{The existence of 
the edge pairs makes |P| a bidirected graph, the |twin| links make
|P| a map.}.

An edge |e = (v,w)| knows two adjacent edges |en = next(e)| and |ep
= previous(e)| where |source(en) = w|, |previous(en) = e| and
|target(ep) = v| and |next(ep) = e|. By this symmetric |previous-next|
relationship all edges are partitioned into face cycles.  Two edges
$e$ and $e'$ are in the same face cycle if $e = |next|^*(e')$.  All
edges |e| in the same face cycle have the same incident face $f =
|face|(e)$. The cyclic order on the adjacency list of a vertex |v =
source(e)| is given by |cyclic_adj_succ(e) = twin(previous(e))| and
|cyclic_adj_pred(e) = next(twin(e))|.

A vertex |v| is embedded via coordinates |point(v)|. By the
embedding of its source and target an edge corresponds to a
segment. |P| has the property that the embedding is always
\emph{order-preserving}.  This means a ray fixed in |point(v)| of a
vertex |v| and swept around counterclockwise meets the embeddings of
|target(e)| ($e \in A(v)$) in the cyclic order defined by the list
order of |A|.

The embedded face cycles partition the plane into maximal connected
subsets of points. Each such set corresponds to a face. A face is
bounded by its incident face cycles. For all the edges in the
non-trivial face cycles it holds that the face is left of the edges.
There can also be trivial face cycles in form of isolated vertices in
the interior of a face. Each such vertex |v| knows its surrounding
face |f = face(v)|.

We call the embedded map |(V,E)| also the $1$-skeleton of |P|.

Plane maps are attributed. To each object $u \in V \cup E \cup F$
we attribute a value |mark(u)| of type |Mark|. |Mark| fits the
concepts assignable, default-constructible, and equal-comparable.}*/

protected: 
HDS* phds; 
friend class PM_decorator<HDS>;

public:
/*{\Mtypes 6}*/
typedef HDS Plane_map;
/*{\Mtypemember The underlying plane map type}*/

typedef typename HDS::Traits Traits;

typedef typename Traits::Point  Point;
/*{\Mtypemember The point type of vertices.}*/
typedef typename Traits::Mark   Mark;
/*{\Mtypemember All objects (vertices, edges, faces) are attributed by a 
|Mark| object.}*/
typedef size_t Size_type;
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
/*{\Mtypemember The size type.}*/
typedef void*  GenPtr;
#else
typedef boost::any GenPtr;
#endif



typedef typename HDS::Vertex                  Vertex; 
typedef typename HDS::Vertex_base             Vertex_base;
typedef typename HDS::Vertex_handle           Vertex_handle;
typedef typename HDS::Vertex_const_handle     Vertex_const_handle;
typedef typename HDS::Vertex_const_iterator   Vertex_const_iterator;
typedef typename HDS::Halfedge                Halfedge ; 
typedef typename HDS::Halfedge_base           Halfedge_base ;
typedef typename HDS::Halfedge_handle         Halfedge_handle; 
typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;
typedef typename HDS::Face                    Face;
typedef typename HDS::Face_base               Face_base;
typedef typename HDS::Face_handle             Face_handle;
typedef typename HDS::Face_const_handle       Face_const_handle;
typedef typename HDS::Face_const_iterator     Face_const_iterator;


/*{\Mtext Local types are handles, iterators and circulators of the
following kind: |Vertex_const_handle|,
|Vertex_const_iterator|, |Halfedge_const_handle|,
|Halfedge_const_iterator|, |Face_const_handle|, |Face_\-const_\-ite\-rator|.
Additionally the following circulators are defined.}*/

typedef CircFromIt<
        Halfedge_const_iterator, 
        move_halfedge_around_vertex<Halfedge_const_iterator> > 
        Halfedge_around_vertex_const_circulator;
/*{\Mtypemember circulating the outgoing halfedges in $A(v)$.}*/

typedef CircFromIt<
        Halfedge_const_iterator, 
        move_halfedge_around_face<Halfedge_const_iterator> > 
        Halfedge_around_face_const_circulator;
/*{\Mtypemember circulating the halfedges in the face cycle of a 
face |f|.}*/

typedef typename Face::Hole_const_iterator Hole_const_iterator;
/*{\Mtypemember iterating all holes of a face |f|. The type is 
convertible to |Halfedge_const_handle|.}*/

typedef typename Face::Isolated_vertex_const_iterator 
        Isolated_vertex_const_iterator;
/*{\Mtypemember iterating all isolated vertices of a face |f|. 
The type generalizes |Vertex_const_handle|.}*/

typedef PntItFromVertIt<Vertex_const_iterator,Point> 
  Point_const_iterator;

/*{\Mcreation 3}*/

PM_const_decorator() : phds(0) {}
PM_const_decorator(const PM_const_decorator<HDS>& D) :
    phds(D.phds) {}
PM_const_decorator& operator=(
  const PM_const_decorator<HDS>& D)
{ phds=D.phds; return *this; }

PM_const_decorator(const Plane_map& P) 
/*{\Mcreate constructs a plane map decorator exploring |P|.}*/
 : phds(const_cast<HDS*>(&P)) {}

/*{\Moperations 4 4}*/

Vertex_const_handle source(Halfedge_const_handle e) const
/*{\Mop returns the source of |e|.}*/
{ return e->opposite()->vertex(); }

Vertex_const_handle target(Halfedge_const_handle e) const
/*{\Mop returns the target of |e|.}*/
{ return e->vertex(); }

Halfedge_const_handle twin(Halfedge_const_handle e) const 
/*{\Mop returns the twin of |e|.}*/
{ return e->opposite(); }

bool is_isolated(Vertex_const_handle v) const
/*{\Mop returns |true| iff $A(v) = \emptyset$.}*/
{ return v->is_isolated(); }

Halfedge_const_handle first_out_edge(Vertex_const_handle v) const
/*{\Mop returns one halfedge with source |v|. It's the starting point for
  the circular iteration over the halfedges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->halfedge()->opposite(); }

Halfedge_const_handle last_out_edge(Vertex_const_handle v) const
/*{\Mop returns the halfedge with source |v| that is the last
  in the circular iteration before encountering |first_out_edge(v)| 
  again. \precond |!is_isolated(v)|.}*/
{ return v->halfedge()->next(); }

Halfedge_const_handle cyclic_adj_succ(
                        Halfedge_const_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->prev()->opposite(); }

Halfedge_const_handle cyclic_adj_pred(
                        Halfedge_const_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->opposite()->next(); }


Halfedge_const_handle next(Halfedge_const_handle e) const
/*{\Mop returns the next edge in the face cycle containing |e|.}*/
{ return e->next(); }

Halfedge_const_handle previous(Halfedge_const_handle e) const
/*{\Mop returns the previous edge in the face cycle containing |e|.}*/
{ return e->prev(); }

Face_const_handle face(Halfedge_const_handle e) const
/*{\Mop returns the face incident to |e|.}*/
{ return e->face(); }

Face_const_handle face(Vertex_const_handle v) const
/*{\Mop returns the face incident to |v|.
   \precond |is_isolated(v)|.}*/
{ return v->face(); }

Halfedge_const_handle halfedge(Face_const_handle f) const
/*{\Mop returns a halfedge in the bounding face cycle of |f| 
(|Halfedge_const_handle()| if there is no bounding face cycle).}*/
{ return f->halfedge(); }


/*{\Mtext \headerline{Iteration} \setopdims{7.5cm}{0cm}}*/
  
Vertex_const_iterator   vertices_begin() const
{ return phds->vertices_begin(); }
Halfedge_const_iterator halfedges_begin() const
{ return phds->halfedges_begin(); }
Face_const_iterator     faces_begin() const
{ return phds->faces_begin(); }
Vertex_const_iterator   vertices_end() const
{ return phds->vertices_end(); }
Halfedge_const_iterator halfedges_end() const
{ return phds->halfedges_end(); }
Face_const_iterator     faces_end() const
{ return phds->faces_end(); }
Point_const_iterator points_begin() const
{ return Point_const_iterator(vertices_begin()); }
Point_const_iterator points_end() const
{ return Point_const_iterator(vertices_end()); }

Halfedge_around_vertex_const_circulator 
  out_edges(Vertex_const_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.}*/
{ return Halfedge_around_vertex_const_circulator(first_out_edge(v)); }

Halfedge_around_face_const_circulator 
  face_cycle(Face_const_handle f) const
/*{\Mop returns a circulator for the outer face cycle of |f|.}*/
{ return Halfedge_around_face_const_circulator(f->halfedge()); }

Hole_const_iterator  holes_begin(Face_const_handle f) const
/*{\Mop returns an iterator for all holes in the interior of |f|.
   A |Hole_iterator| can be assigned to a 
   |Halfedge_around_face_const_circulator|.}*/
{ return f->fc_begin(); }

Hole_const_iterator  holes_end(Face_const_handle f) const
/*{\Mop returns the past-the-end iterator of |f|.}*/
{ return f->fc_end(); }

Isolated_vertex_const_iterator 
  isolated_vertices_begin(Face_const_handle f) const
/*{\Mop returns an iterator for all isolated vertices in the interior 
   of |f|.}*/
{ return f->iv_begin(); }

Isolated_vertex_const_iterator 
  isolated_vertices_end(Face_const_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->iv_end(); }


/*{\Mtext \headerline{Associated Information}\setopdims{2.5cm}{4cm}
The type |Mark| is the general attribute of an object. The type
|GenPtr| is equal to type |void*|.}*/

const Point& point(Vertex_const_handle v) const
/*{\Mop returns the embedding of |v|.}*/
{ return v->point(); }

const Mark& mark(Vertex_const_handle v) const  
/*{\Mop returns the mark of |v|.}*/
{ return v->mark(); }

const Mark& mark(Halfedge_const_handle e) const 
/*{\Mop returns the mark of |e|.}*/
{ if (&*e < &*(e->opposite())) return e->mark();
  else return e->opposite()->mark();
} // we store the mark in the container with smaller memory address !

const Mark& mark(Face_const_handle f) const
/*{\Mop returns the mark of |f|.}*/
{ return f->mark(); }

const GenPtr& info(Vertex_const_handle v) const
/*{\Mop returns a generic information slot.}*/
{ return v->info(); }

const GenPtr& info(Halfedge_const_handle e) const
/*{\Mop returns a generic information slot.}*/
{ return e->info(); }

const GenPtr& info(Face_const_handle f) const
/*{\Mop returns a generic information slot.}*/
{ return f->info(); }


/*{\Mtext \headerline{Statistics and Integrity}}*/

Size_type number_of_vertices() const 
/*{\Mop returns the number of vertices.}*/
{ return phds->size_of_vertices(); }

Size_type number_of_halfedges() const 
/*{\Mop returns the number of halfedges.}*/
{ return phds->size_of_halfedges(); }

Size_type number_of_edges() const    
/*{\Mop returns the number of halfedge pairs.}*/
{ return phds->size_of_halfedges()/2; }

Size_type number_of_faces() const    
/*{\Mop returns the number of faces.}*/
{ return phds->size_of_faces(); }

Size_type number_of_face_cycles() const;
/*{\Mop returns the number of face cycles.}*/

Size_type number_of_connected_components() const;
/*{\Mop calculates the number of connected components of |P|.}*/

void print_statistics(std::ostream& os = std::cout) const
/*{\Mop print the statistics of |P|: the number of vertices, edges, and 
   faces.}*/
{
  os << "Plane Map - Statistics\n";
  os << "|V| = " << number_of_vertices() << std::endl;
  os << "|E| = " << number_of_edges() 
  << " (" << 2*number_of_edges() << ")" << std::endl;
  os << "|F| = " << number_of_faces() << std::endl;
  os << "|Fcs| = " << number_of_face_cycles() << std::endl << std::endl;
}
 
void check_integrity_and_topological_planarity(bool faces=true) const;
/*{\Mop checks the link structure and the genus of |P|.}*/

}; // PM_const_decorator

template <class VH>
std::string PV(VH v)
{ std::ostringstream os; CGAL::set_pretty_mode(os);
  if (v != VH()) os << v->point();
  else           os << "nil";
  return os.str();
}

template <class HH>
std::string PE(HH e)
{ std::ostringstream os;
  if (e==HH()) return "nil";
  os << "[" << PV(e->opposite()->vertex()) << ","
            << PV(e->vertex()) << " " 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  << e->info()
  #endif
  << "]";
  return os.str();
}

template <typename HDS>
void PM_const_decorator<HDS>::
check_integrity_and_topological_planarity(bool faces) const
{
  CGAL_NEF_TRACEN("check_integrity_and_topological_planarity:");
  using CGAL::Object_index;
  Object_index<Vertex_const_iterator>   
    VI(vertices_begin(),vertices_end(),'v');
  Object_index<Halfedge_const_iterator> 
    EI(halfedges_begin(),halfedges_end(),'e');
  Object_index<Face_const_iterator> 
    FI(faces_begin(),faces_end(),'f');
  Vertex_const_handle vit, vend = phds->vertices_end();
  int iso_vert_num=0;
  /* check the source links of out edges and count isolated vertices */
  for (vit = vertices_begin() ; vit != vend; ++vit) {
    if ( is_isolated(vit) ) {
      if ( faces )
        CGAL_assertion_msg( vit->face() != Face_const_handle(),
                            VI(vit).c_str());
      ++iso_vert_num;
    } else {
      CGAL_assertion_msg( vit->halfedge() != Halfedge_const_handle(),
      VI(vit).c_str());
      CGAL_assertion_msg( vit->halfedge()->vertex() == vit ,VI(vit).c_str());
    }
  }

  /* check the bidirected links and the face pointer init */
  Halfedge_const_iterator eit, eend = phds->halfedges_end();
  for (eit = phds->halfedges_begin() ; eit != eend; ++eit) {
    CGAL_assertion( twin(twin(eit)) == eit );
    CGAL_assertion( eit->vertex() != Vertex_const_handle() );
    CGAL_assertion( next(eit) != Halfedge_const_handle() );
    CGAL_assertion( previous(next(eit)) == eit );
    CGAL_assertion( target(eit) == source(next(eit)) );
    CGAL_assertion( previous(eit) != Halfedge_const_handle() );
    CGAL_assertion( next(previous(eit)) == eit );
    CGAL_assertion( target(previous(eit)) == source(eit) );
    if ( !faces ) continue;
    CGAL_assertion( face(eit) != Face_const_handle() );
    CGAL_assertion( face(next(eit)) == face(eit) );
    CGAL_assertion( face(previous(eit)) == face(eit) );
  }

  bool first=true;
  int fc_num(0),iv_num(0);
  Face_const_iterator fit;
  for (fit = faces_begin(); fit != faces_end(); ++fit) {
    if (!first) {
      CGAL_assertion( face(halfedge(fit))==fit ); ++fc_num;
    }
    Hole_const_iterator fcit;
    for( fcit = holes_begin(fit); fcit != holes_end(fit); ++fcit) {
      CGAL_assertion( face(fcit)==fit ); ++fc_num;
    }
    Isolated_vertex_const_iterator ivit;
    for(ivit = isolated_vertices_begin(fit); 
        ivit != isolated_vertices_end(fit); ++ivit) {
      CGAL_assertion( face(ivit)==fit ); ++iv_num;
    }
    first=false;
  }

  std::size_t v_num = number_of_vertices() - iso_vert_num;
  std::size_t e_num = number_of_edges();
  std::size_t c_num = number_of_connected_components() - iso_vert_num;
  std::size_t f_num = number_of_face_cycles() - c_num + 1;
  CGAL_NEF_TRACEV(fc_num);CGAL_NEF_TRACEV(iv_num);CGAL_NEF_TRACEV(iso_vert_num);
  CGAL_NEF_TRACEV(v_num);CGAL_NEF_TRACEV(e_num);CGAL_NEF_TRACEV(c_num);CGAL_NEF_TRACEV(f_num);
  // CGAL_assertion(fc_num == f_num && iv_num == iso_vert_num);
  /* this means all face cycles and all isolated vertices are 
     indeed referenced from a face */
  /* every isolated vertex increases the component count
       one face cycle per component is redundent except one
       finally check the Euler formula: */
  CGAL_assertion( v_num - e_num + f_num == 1 + c_num );
}

template <typename HDS>
typename PM_const_decorator<HDS>::Size_type 
PM_const_decorator<HDS>::
number_of_face_cycles() const
{
  Size_type fc_num=0;
  CGAL::Unique_hash_map<Halfedge_const_handle,bool> visited; 
    // init with bool() == false
  Halfedge_const_iterator eit =  phds->halfedges_begin();
  Halfedge_const_iterator eend = phds->halfedges_end();
  for ( ; eit != eend; ++eit) {
    if (visited[eit]) continue;
    Halfedge_around_face_const_circulator hfc(eit), hend(hfc);
    CGAL_For_all(hfc,hend) visited[hfc]=true;
    ++fc_num;
  }
  return fc_num;
}

template <typename HDS>
size_t PM_const_decorator<HDS>::
number_of_connected_components() const
{
  typedef Vertex_const_iterator vc_handle;
  typedef Halfedge_around_vertex_const_circulator hvc_circulator;
  int comp_num=0;
  CGAL::Unique_hash_map< vc_handle, bool> handled(false); 
  vc_handle vit = vertices_begin(), vend = vertices_end();
  for ( ; vit != vend; ++vit) {
    if (handled[vit]) continue;
    std::list<vc_handle> L;
    L.push_back(vit); handled[vit]=true; 
    /* we keep the invariant that all nodes which have been stacked
       are marked handled */
    while (!L.empty()) {
      vc_handle v=L.front(); L.pop_front();
      if ( is_isolated(v) ) continue;
      hvc_circulator havc(first_out_edge(v)), hend(havc);
      CGAL_For_all(havc,hend) {
        if (!handled[target(havc)]) {
          L.push_back(target(havc)); handled[target(havc)]=true; 
        }
      }
    }
    ++comp_num;
  }
  return comp_num;   
}

struct KERNELPNT {
  template <typename PNT>
  std::string operator() (const PNT& p) const
  { std::ostringstream os;
    os << "(" << CGAL::to_double(p.x()) << ","
              << CGAL::to_double(p.y()) << ")";
    return os.str();
  }
};

template <typename PMCDEC, typename POINTDA>
void print_as_leda_graph(std::ostream& os, const PMCDEC& D, 
  const POINTDA& P)
{
  typedef typename PMCDEC::Vertex_const_iterator   
  Vertex_const_iterator;
  typedef typename PMCDEC::Halfedge_const_iterator 
  Halfedge_const_iterator;
  int vn(1), en(1);  
  CGAL::Unique_hash_map<Vertex_const_iterator,int>   v_num;
  CGAL::Unique_hash_map<Halfedge_const_iterator,int> e_num;
  os << "LEDA.GRAPH\n" << "point\n" << "int\n";
  os << D.number_of_vertices() << std::endl;
  Vertex_const_iterator vit;
  for (vit = D.vertices_begin(); vit != D.vertices_end(); ++vit) {
    v_num[vit] = vn++;
    os << "|{(" << P(D.point(vit)) << ")}|\n";
     typename PMCDEC::Halfedge_around_vertex_const_circulator
       ecirc(D.first_out_edge(vit)),ecend(ecirc);
    int l=0;
    CGAL_For_all(ecirc,ecend) e_num[ecirc]=l++;
  }
  os << 2* D.number_of_edges() << std::endl;
  Halfedge_const_iterator eit;
  for (eit = D.halfedges_begin(); eit != D.halfedges_end(); ++eit) {
    e_num[eit] = en++;
  }
  for (eit = D.halfedges_begin(); eit != D.halfedges_end(); ++eit) {
    os << v_num[D.source(eit)] << " "
       << v_num[D.target(eit)] << " "
       << e_num[D.twin(eit)]   << " ";
    os << "|{" << e_num[eit] << "}|\n";
  }
  os << std::flush;
}



} //namespace CGAL
#endif // CGAL_PM_CONST_DECORATOR_H
