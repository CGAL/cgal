#ifndef CGAL_SM_CONST_DECORATOR_H 
#define CGAL_SM_CONST_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <string>
#include <list>
#include <strstream>
#undef _DEBUG
#define _DEBUG 127
#include <CGAL/Nef_S2/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename HE>
class move_edge_around_vertex {
public:
  void forward(HE& e) const  { e = (e->prev_->twin_); }
  void backward(HE& e) const { e = (e->twin_->next_); }
};

template <typename HE>
struct move_edge_around_face {
  void forward(HE& e)  const { e = (e->next_); }
  void backward(HE& e) const { e = (e->prev_); }
};


template <typename Sphere_map_, typename Kernel_>
class SM_decorator;

/*{\Moptions print_title=yes }*/ 
/*{\Moptions outfile=SM_const_decorator.man }*/
/*{\Manpage {SM_const_decorator}{Sphere_map,Kernel}
  {Topological sphere map decorator}{D}}*/

template <typename Sphere_map_, typename Kernel_>
class SM_const_decorator 
{ typedef SM_const_decorator<Sphere_map_,Kernel_> Self;
public:

/*{\Mdefinition ...}*/

/*{\Mtypes 5}*/

typedef Kernel_ Kernel;

typedef typename Kernel_::Sphere_point Sphere_point;
/*{\Mtypemember embedding vertices.}*/

typedef typename Kernel_::Sphere_segment Sphere_segment;
/*{\Mtypemember embedding edges.}*/

typedef typename Kernel_::Sphere_circle Sphere_circle;
/*{\Mtypemember embedding loops.}*/

typedef typename Sphere_map_::Mark   Mark;
/*{\Mtypemember attributes of objects (vertices, edges, faces).}*/

typedef size_t Size_type;
/*{\Mtypemember size type.}*/

typedef void*  GenPtr;

#define USING(t) typedef typename Sphere_map_::t t
USING(Vertex_const_handle);
USING(Vertex_const_iterator);
USING(Vertex); 
USING(Halfedge); 
USING(Halfedge_const_handle); 
USING(Halfedge_const_iterator);
USING(Face);
USING(Face_const_handle);
USING(Face_const_iterator);
USING(Halfloop);
USING(Halfloop_const_handle);
USING(Halfloop_const_iterator);
#undef USING

/*{\Mtext Local types are handles, iterators and circulators of the
following kind: |Vertex_handle|, |Vertex_iterator|, |Halfedge_handle|,
|Halfedge_iterator|, |Halfloop_handle|, |Halfloop_iterator|,
|Face_handle|, |Face_iterator|.  Additionally the following
circulators are defined.}*/

typedef CircFromIt<
        Halfedge_const_iterator, 
        move_edge_around_vertex<Halfedge_const_iterator> > 
        Halfedge_around_vertex_const_circulator;
/*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/

typedef CircFromIt<
        Halfedge_const_iterator, 
        move_edge_around_face<Halfedge_const_iterator> > 
        Halfedge_around_face_const_circulator;
/*{\Mtypemember circulating the face cycle of an face |f|.}*/

typedef typename Sphere_map_::Face_cycle_const_iterator 
  Face_cycle_const_iterator;
/*{\Mtypemember iterating all face cycles of an face |f|.
The iterator has method |bool is_vertex()|, |bool is_halfedge()|,
|bool is_halfloop()|, and can be converted to the corresponding
handles |Vertex_const_handle|, |Halfedge_const_handle|, or 
|Halfloop_const_handle|.}*/

protected: 
  Sphere_map_* psm_;
  friend class SM_decorator<Sphere_map_,Kernel_>;
public:

/*{\Mcreation 3}*/
SM_const_decorator() : psm_(0) {}
SM_const_decorator(const Self& D) : psm_(D.psm_) {}
Self& operator=(const Self& D) { psm_=D.psm_; return *this; }

SM_const_decorator(const Sphere_map_& M) : 
  psm_(const_cast<Sphere_map_*>(&M)) {}
/*{\Mcreate constructs a plane map decorator exploring |M|.}*/

/*{\Moperations 4 4}*/

Vertex_const_handle source(Halfedge_const_handle e) const
/*{\Mop returns the source of |e|.}*/
{ return e->source_; }

Vertex_const_handle target(Halfedge_const_handle e) const
/*{\Mop returns the target of |e|.}*/
{ return e->twin_->source_; }

Halfedge_const_handle twin(Halfedge_const_handle e) const 
/*{\Mop returns the twin of |e|.}*/
{ return e->twin_; }

Halfloop_const_handle twin(Halfloop_const_handle l) const 
/*{\Mop returns the twin of |l|.}*/
{ return l->twin_; }

bool is_isolated(Vertex_const_handle v) const
/*{\Mop returns |true| when |v| is linked to the interior of a face.}*/
{ return (Halfedge_const_handle(v->edge_) == Halfedge_const_handle()); }

Halfedge_const_handle first_out_edge(Vertex_const_handle v) const
/*{\Mop returns one edge with source |v|. It's the starting point for
  the circular iteration over the edges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->edge_; }

Halfedge_const_handle last_out_edge(Vertex_const_handle v) const
/*{\Mop returns one edge with source |v|. \precond |!is_isolated(v)|.}*/
{ return cap(v->edge_); }

Halfedge_const_handle cyclic_adj_succ(Halfedge_const_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->prev_->twin_; }

Halfedge_const_handle cyclic_adj_pred(Halfedge_const_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->twin_->next_; }


Halfedge_const_handle next(Halfedge_const_handle e) const
/*{\Mop returns the next edge in the face cycle containing |e|.}*/
{ return e->next_; }

Halfedge_const_handle previous(Halfedge_const_handle e) const
/*{\Mop returns the previous edge in the face cycle containing |e|.}*/
{ return e->prev_; }

Face_const_handle face(Halfedge_const_handle e) const
/*{\Mop returns the face incident to |e|.}*/
{ return e->face_; }

Face_const_handle face(Halfloop_const_handle l) const
/*{\Mop returns the face incident to |l|.}*/
{ return l->face_; }

Face_const_handle face(Vertex_const_handle v) const
/*{\Mop returns the face incident to |v|.
   \precond |is_isolated(v)|.}*/
{ return v->face_; }

/*{\Mtext \headerline{Iteration} \setopdims{3.3cm}{0cm}}*/
  
Vertex_const_iterator vertices_begin() const
{ return psm_->vertices_begin(); }
Vertex_const_iterator vertices_end() const
{ return psm_->vertices_end(); }
Halfedge_const_iterator halfedges_begin() const
{ return psm_->halfedges_begin(); }
Halfedge_const_iterator halfedges_end() const
{ return psm_->halfedges_end(); }
Face_const_iterator faces_begin() const
{ return psm_->faces_begin(); }
Face_const_iterator faces_end() const
{ return psm_->faces_end(); }
Halfloop_const_iterator halfloops_begin() const
{ return psm_->halfloops_begin(); }
Halfloop_const_iterator halfloops_end() const
{ return psm_->halfloops_end(); }

Halfloop_const_handle halfloop() const
/*{\Mop returns access to the loop.}*/
{ return psm_->loops_; }

bool has_loop() const
/*{\Mop returns true iff there is a loop.}*/
{ return halfloop() != 0; }

Halfedge_around_vertex_const_circulator 
  out_edges(Vertex_const_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.
\precond the adjacency list is not empty.}*/
{ return Halfedge_around_vertex_const_circulator(first_out_edge(v)); }

Face_cycle_const_iterator face_cycles_begin(Face_const_handle f) const
/*{\Mop returns an iterator for all bounding face cycles of |f|.
The iterator is is convertable to |Vertex_const_handle|, 
|Halfloop_const_handle|, or |Halfedge_const_handle|.}*/
{ return f->bounday_.begin(); }

Face_cycle_const_iterator face_cycles_end(Face_const_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->boundary_.end(); }

/*{\Mtext \headerline{Statistics and Integrity}}*/

Size_type number_of_vertices() const 
/*{\Mop returns the number of vertices.}*/
{ return psm_->number_of_vertices(); }

Size_type number_of_halfedges() const 
/*{\Mop returns the number of halfedges.}*/
{ return psm_->number_of_halfedges(); }

Size_type number_of_edges() const 
/*{\Mop returns the number of edges.}*/
{ return number_of_halfedges()/2; }

Size_type number_of_halfloops() const 
/*{\Mop returns the number of halfloops.}*/
{ return psm_->number_of_halfloops(); }

Size_type number_of_loops() const 
/*{\Mop returns the number of loops.}*/
{ return psm_->number_of_halfloops()/2; }

Size_type number_of_faces() const    
/*{\Mop returns the number of faces.}*/
{ return psm_->number_of_faces(); }

Size_type number_of_face_cycles() const;
/*{\Mop returns the number of face cycles.}*/

Size_type number_of_connected_components() const;
/*{\Mop calculates the number of connected components of |P|.}*/

void print_statistics(std::ostream& os = std::cout) const
/*{\Mop print the statistics of |P|: the number of vertices, edges, 
and faces.}*/
{
  os << "Sphere Map - Statistics\n";
  os << "|V| = " << number_of_vertices() << std::endl;
  os << "|E| = " << number_of_edges() << std::endl;
  os << "|L| = " << number_of_halfloops()/2 << std::endl;
  os << "|F| = " << number_of_faces() << std::endl;
  os << "|Fcs| = " << number_of_face_cycles() << std::endl << std::endl;
}
 
void check_integrity_and_topological_planarity(bool faces=true) const;
/*{\Mop checks the link structure and the genus of |P|.}*/

Halfedge_const_handle cas(Halfedge_const_handle e) const 
{ return cyclic_adj_succ(e); } 

Halfedge_const_handle cap(Halfedge_const_handle e) const
{ return cyclic_adj_pred(e); }

template <typename H>
bool is_boundary_object(H h) const
{ return psm_->is_boundary_object(h); }

/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

const Sphere_point& point(Vertex_const_handle v) const
/*{\Mop returns the embedding of |v|.}*/
{ return v->point_; }

const Sphere_circle& circle(Halfedge_const_handle e) const
/*{\Mop returns the circle supporting |e|.}*/
{ return e->circle_; }

const Sphere_circle& circle(Halfloop_const_handle l) const
/*{\Mop returns the circle supporting |l|.}*/
{ return l->circle_; }

Mark mark(Vertex_const_handle v) const
/*{\Mop returns the mark of |v|.}*/
{ return v->mark_; }

Mark mark(Halfedge_const_handle e) const
/*{\Mop returns the mark of |e|.}*/
{ return ( &*e < &*twin(e) ) ? e->mark_ : twin(e)->mark_; }

Mark mark(Halfloop_const_handle l) const
/*{\Mop returns the mark of |l|.}*/
{ return ( &*l < &*twin(l) ) ? l->mark_ : twin(l)->mark_; }

Mark mark(Face_const_handle f) const
/*{\Mop returns the mark of |f|.}*/
{ return f->mark_; }

Mark mark_of_halfsphere(int i) const
{ CGAL_nef_assertion(i);
  if (i<0) return psm_->m_neg_; 
  return psm_->m_pos_; }

/*{\Mtext \headerline{Iteration}}*/
/*{\Mtext The list of all objects can be accessed via iterator ranges.
For comfortable iteration we also provide iterations macros. 
The iterator range access operations are of the following kind:\\
|Vertex_iterator   vertices_begin()/vertices_end()|\\
|Halfedge_iterator halfedges_begin()/halfedges_end()|\\
|Halfloop_iterator halfloops_begin()/halfloops_end()|\\
|Face_iterator     faces_begin()/faces_end()|

The macros are then |CGAL_forall_vertices(v,M)|,
|CGAL_forall_halfedges(e,M)|, |CGAL_forall_edges(e,M)|,
|CGAL_forall_faces(f,M)|, |CGAL_forall_face_cycles_of(fc,F)| where |M|
is a sphere map and |F| is a face.}*/

}; // SM_const_decorator



template <typename SM_, typename K_>
void SM_const_decorator<SM_,K_>::
check_integrity_and_topological_planarity(bool faces) const
{
  TRACEN("check_integrity_and_topological_planarity:");
  using CGAL::Object_index;
  Object_index<Vertex_const_iterator>   
    VI(vertices_begin(),vertices_end(),'v');
  Object_index<Halfedge_const_iterator> 
    EI(halfedges_begin(),halfedges_end(),'e');
  Object_index<Face_const_iterator> 
    FI(faces_begin(),faces_end(),'f');
  typedef Halfedge_around_vertex_const_circulator hvc_circulator;
  typedef Halfedge_around_face_const_circulator   hfc_circulator;
  Vertex_const_handle v;
  int iso_vert_num=0;
  /* check the source links of out edges and count isolated vertices */
  CGAL_forall_vertices(v,*this) {
    if ( is_isolated(v) ) {
      if ( faces )
        CGAL_nef_assertion_msg( face(v) != Face_const_handle(), VI(v).c_str());
      ++iso_vert_num;
    } else {
      CGAL_nef_assertion_msg( first_out_edge(v) != Halfedge_const_handle(),
      VI(v).c_str());
      TRACEN(point(v)<<" "<<EI(first_out_edge(v)));
      CGAL_nef_assertion_msg( source(first_out_edge(v)) == v ,
                          VI(v).c_str());
    }
  }

  /* check the bidirected links and the face pointer init */
  Halfedge_const_iterator e;
  CGAL_forall_halfedges(e,*this) {
    CGAL_nef_assertion( twin(twin(e)) == e );
    CGAL_nef_assertion( source(e) != Vertex_const_handle() );
    CGAL_nef_assertion( next(e) != Halfedge_const_handle() );
    CGAL_nef_assertion( previous(next(e)) == e );
    CGAL_nef_assertion( target(e) == source(next(e)) );
    CGAL_nef_assertion( previous(e) != Halfedge_const_handle() );
    CGAL_nef_assertion( next(previous(e)) == e );
    CGAL_nef_assertion( target(previous(e)) == source(e) );
    if ( !faces ) continue;
    CGAL_nef_assertion( face(e) != Face_const_handle() );
    CGAL_nef_assertion( face(next(e)) == face(e) );
    CGAL_nef_assertion( face(previous(e)) == face(e) );
  }

  int fc_num(0),iv_num(0);
  Face_const_iterator f;
  Face_cycle_const_iterator fci;
  CGAL_forall_faces(f,*this) {
    CGAL_forall_face_cycles_of(fci,f) {
      if ( fci.is_halfedge() ) {
        CGAL_nef_assertion( face(Halfedge_const_handle(fci)) == f ); ++fc_num;
      } else if ( fci.is_vertex() ) {
        CGAL_nef_assertion( face(Vertex_const_handle(fci)) == f ); ++iv_num;
      } else CGAL_nef_assertion_msg(0,"damn generic handle.");
    }
  }

  int v_num = number_of_vertices() - iso_vert_num;
  int e_num = number_of_edges();
  int c_num = number_of_connected_components() - iso_vert_num;
  int f_num = number_of_face_cycles() - c_num + 1;
  TRACEV(fc_num);TRACEV(iv_num);TRACEV(iso_vert_num);
  TRACEV(v_num);TRACEV(e_num);TRACEV(c_num);TRACEV(f_num);
  /* this means all face cycles and all isolated vertices are 
     indeed referenced from a face */
  /* every isolated vertex increases the component count
       one face cycle per component is redundent except one
       finally check the Euler formula: */
  CGAL_nef_assertion( v_num - e_num + f_num == 1 + c_num );
}

template <typename SM_, typename K_>
typename SM_const_decorator<SM_,K_>::Size_type
SM_const_decorator<SM_,K_>::
number_of_face_cycles() const
{
  unsigned int fc_num=0;
  CGAL::Unique_hash_map<Halfedge_const_handle,bool> visited; 
    // init with bool() == false
  Halfedge_const_iterator e;
  CGAL_forall_halfedges(e,*this) {
    if (visited[e]) continue;
    Halfedge_around_face_const_circulator hfc(e), hend(hfc);
    CGAL_For_all(hfc,hend) visited[hfc]=true;
    ++fc_num;
  }
  if ( has_loop() ) fc_num += 2;
  return fc_num;
}

template <typename SM_, typename K_>
typename SM_const_decorator<SM_,K_>::Size_type
SM_const_decorator<SM_,K_>::
number_of_connected_components() const
{
  int comp_num=0;
  CGAL::Unique_hash_map<Vertex_const_iterator,bool> visited(false); 
  Vertex_const_iterator v;
  CGAL_forall_vertices(v,*this) {
    if (visited[v]) continue;
    std::list<Vertex_const_iterator> L;
    L.push_back(v); visited[v]=true; 
    /* we keep the invariant that all nodes which have been stacked
       are marked visited */
    while (!L.empty()) {
      Vertex_const_iterator vc = L.front(); L.pop_front();
      if ( is_isolated(vc) ) continue;
      Halfedge_around_vertex_const_circulator 
	havc(first_out_edge(vc)), hend(havc);
      CGAL_For_all(havc,hend) {
        if (!visited[target(havc)]) {
          L.push_back(target(havc)); visited[target(havc)]=true; 
        }
      }
    }
    ++comp_num;
  }
  return comp_num;   
}

CGAL_END_NAMESPACE
#endif // CGAL_SM_CONST_DECORATOR_H

