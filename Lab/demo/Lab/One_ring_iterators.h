#ifndef ONE_RING_ITERATORS_H
#define ONE_RING_ITERATORS_H

/***************************************************************/
/* For iterating on one ring neighbor, sample usage:
   Halfedge_handle h; //or Vertex_handle or Facet_handle
   for(One_ring_iterator<Halfedge_handle> circ(h); circ; ++circ) {
     Halfedge_handle h_neighbor = circ;
   }
*/
template<typename Mesh, class T>
struct One_ring_iterator;

template<typename Mesh>
struct One_ring_iterator<Mesh, typename boost::graph_traits<Mesh>::vertex_descriptor> {
  One_ring_iterator(typename boost::graph_traits<Mesh>::vertex_descriptor v, const Mesh& mesh)
    : circ(v, mesh), end(circ), first(true), mesh(mesh) { }

  operator bool() const { return first || circ != end; }
  operator typename boost::graph_traits<Mesh>::vertex_descriptor() const { return target(opposite(*circ, mesh), mesh); }
  One_ring_iterator& operator++() {
    first = false;
    ++circ;
    return *this;
  }

  CGAL::Halfedge_around_target_circulator<Mesh> circ;
  CGAL::Halfedge_around_target_circulator<Mesh> end;
  bool first;
  const Mesh& mesh;
  // to be used in One_ring_iterator<Halfedge_handle>
  operator typename boost::graph_traits<Mesh>::halfedge_descriptor() const {
    typename boost::graph_traits<Mesh>::halfedge_descriptor op =
        opposite(*circ, mesh);
    if ( &*circ < &op ) return *circ;
    return op;
  }
};

template<typename Mesh>
struct One_ring_iterator<Mesh, typename boost::graph_traits<Mesh>::face_descriptor> {
  One_ring_iterator(typename boost::graph_traits<Mesh>::face_descriptor f, const Mesh& mesh)
    : circ(halfedge(f, mesh), mesh), end(circ), first(true), mesh(mesh)
  {
    iterate_to_non_border(); // move it to valid location
  }

  operator bool() const { return first || circ != end; }
  operator typename boost::graph_traits<Mesh>::face_descriptor() const {
    CGAL_assertion(!is_border_edge(opposite(*circ, mesh), mesh));
    return face(opposite(*circ, mesh), mesh);
  }
  One_ring_iterator& operator++() {
    first = false;
    ++circ;
    if(circ != end) { iterate_to_non_border(); }
    return *this;
  }

  void iterate_to_non_border() {
    while(is_border_edge(opposite(*circ, mesh), mesh)) {
      first = false;
      ++circ;
      if(circ == end) { break; }
    }
  }
  CGAL::Halfedge_around_face_circulator<Mesh> circ;
  CGAL::Halfedge_around_face_circulator<Mesh> end;
  bool first;
  const Mesh& mesh;
};

template<typename Mesh>
struct One_ring_iterator<Mesh, typename boost::graph_traits<Mesh>::edge_descriptor> {
  One_ring_iterator(typename boost::graph_traits<Mesh>::edge_descriptor e, const Mesh& mesh)
    : it_1(target(halfedge(e, mesh), mesh), mesh), it_2(target(opposite(halfedge(e, mesh), mesh), mesh), mesh),
           mesh(mesh)
  { }

  operator bool() const { return it_1 || it_2; }
  operator typename boost::graph_traits<Mesh>::edge_descriptor() const {
     Mesh tmp;
     if(it_1) { return edge(it_1, tmp); }
     return edge(it_2, tmp);
  }
  One_ring_iterator& operator++() {
    it_1 ? ++it_1 : ++it_2;
    return *this;
  }

  One_ring_iterator<Mesh, typename boost::graph_traits<Mesh>::vertex_descriptor> it_1;
  One_ring_iterator<Mesh, typename boost::graph_traits<Mesh>::vertex_descriptor> it_2;
  const Mesh& mesh;
};

#endif
