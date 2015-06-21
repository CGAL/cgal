#ifndef ONE_RING_ITERATORS_H
#define ONE_RING_ITERATORS_H
#include "Polyhedron_type.h"
/***************************************************************/
/* For iterating on one ring neighbor, sample usage:
   Halfedge_handle h; //or Vertex_handle or Facet_handle
   for(One_ring_iterator<Halfedge_handle> circ(h); circ; ++circ) {
     Halfedge_handle h_neighbor = circ;
   } 
*/
template<class T>
struct One_ring_iterator;

template<>
struct One_ring_iterator<Polyhedron::Vertex_handle> {
  One_ring_iterator(Polyhedron::Vertex_handle v) 
    : circ(v->vertex_begin()), end(circ), first(true) { }

  operator bool() const { return first || circ != end; }
  operator Polyhedron::Vertex_handle() const { return circ->opposite()->vertex(); }
  One_ring_iterator& operator++() { 
    first = false;
    ++circ; 
    return *this;
  }

  Polyhedron::Halfedge_around_vertex_circulator circ;
  Polyhedron::Halfedge_around_vertex_circulator end;
  bool first;
  // to be used in One_ring_iterator<Halfedge_handle>
  operator Polyhedron::Halfedge_handle() const { 
    if ( &*circ < &*circ->opposite() ) return circ;
    return circ->opposite();
  }
};

template<>
struct One_ring_iterator<Polyhedron::Facet_handle> {
  One_ring_iterator(Polyhedron::Facet_handle f)
    : circ(f->facet_begin()), end(circ), first(true)
  {
    iterate_to_non_border(); // move it to valid location
  }

  operator bool() const { return first || circ != end; }
  operator Polyhedron::Facet_handle() const {
    CGAL_assertion(!circ->opposite()->is_border());
    return circ->opposite()->facet(); 
  }
  One_ring_iterator& operator++() {
    first = false;
    ++circ;
    if(circ != end) { iterate_to_non_border(); }
    return *this;
  }
  
  void iterate_to_non_border() {
    while(circ->opposite()->is_border()) {
      first = false;
      ++circ;
      if(circ == end) { break; }
    }
  }
  Polyhedron::Halfedge_around_facet_circulator circ;
  Polyhedron::Halfedge_around_facet_circulator end;
  bool first;
};

template<>
struct One_ring_iterator<boost::graph_traits<Polyhedron>::edge_descriptor> {
  One_ring_iterator(boost::graph_traits<Polyhedron>::edge_descriptor h)
    : it_1(h.halfedge()->vertex()), it_2(h.halfedge()->opposite()->vertex())
  { }

  operator bool() const { return it_1 || it_2; }
  operator boost::graph_traits<Polyhedron>::edge_descriptor() const {
     Polyhedron tmp;
     if(it_1) { return edge(it_1, tmp); }
     return edge(it_2, tmp);
  }
  One_ring_iterator& operator++() {
    it_1 ? ++it_1 : ++it_2;
    return *this;
  }

  One_ring_iterator<Polyhedron::Vertex_handle> it_1;
  One_ring_iterator<Polyhedron::Vertex_handle> it_2;
};

#endif
