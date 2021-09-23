#ifndef _PRINT_ARR_H_
#define _PRINT_ARR_H_

//-----------------------------------------------------------------------------
// Print all neighboring vertices to a given arrangement vertex.
//
template<class Arrangement>
void print_neighboring_vertices(typename Arrangement::Vertex_const_handle v) {
  if (v->is_isolated()) {
    std::cout << "The vertex (" << v->point() << ") is isolated\n";
    return;
  }

  std::cout << "The neighbors of the vertex (" << v->point() << ") are:";
  typename Arrangement::Halfedge_around_vertex_const_circulator first, curr;
  first = curr = v->incident_halfedges();
  do std::cout << " (" << curr->source()->point() << ")";
  while (++curr != first);
  std::cout << std::endl;
}

//-----------------------------------------------------------------------------
// Print all vertices (points) and edges (curves) along a connected component
// boundary.
//
template <typename Arrangement>
void print_ccb(typename Arrangement::Ccb_halfedge_const_circulator circ) {
  std::cout << "(" << circ->source()->point() << ")";
  typename Arrangement::Ccb_halfedge_const_circulator curr = circ;
  do {
    typename Arrangement::Halfedge_const_handle e = curr;
    std::cout << "   [" << e->curve() << "]   "
              << "(" << e->target()->point() << ")";
  } while (++curr != circ);
  std::cout << std::endl;
}

//-----------------------------------------------------------------------------
// Print the boundary description of an arrangement face.
//
template <typename Arrangement>
void print_face(typename Arrangement::Face_const_handle f) {
  // Print the outer boundary.
  if (f->is_unbounded()) std::cout << "Unbounded face.\n";
  else {
    std::cout << "Outer boundary: ";
    print_ccb<Arrangement>(f->outer_ccb());
  }

  // Print the boundary of each of the holes.
  size_t index = 1;
  for (auto hole = f->holes_begin(); hole != f->holes_end(); ++hole, ++index) {
    std::cout << "    Hole #" << index << ": ";
    // The following statement pacifies msvc.
    typename Arrangement::Ccb_halfedge_const_circulator circ = *hole;
    print_ccb<Arrangement>(circ);
  }

  // Print the isolated vertices.
  index = 1;
  for (auto iv = f->isolated_vertices_begin();
       iv != f->isolated_vertices_end(); ++iv, ++index)
  {
    std::cout << "    Isolated vertex #" << index << ": "
              << "(" << iv->point() << ")" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Print the given arrangement.
//
template <typename Arrangement>
void print_arrangement(const Arrangement& arr) {
  CGAL_precondition(arr.is_valid());

  // Print the arrangement vertices.
  std::cout << arr.number_of_vertices() << " vertices:\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "(" << vit->point() << ")";
    if (vit->is_isolated()) std::cout << " - Isolated.\n";
    else std::cout << " - degree " << vit->degree() << std::endl;
  }

  // Print the arrangement edges.
  std::cout << arr.number_of_edges() << " edges:\n";
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
    std::cout << "[" << eit->curve() << "]\n";

  // Print the arrangement faces.
  std::cout << arr.number_of_faces() << " faces:\n";
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    print_face<Arrangement>(fit);
}

//-----------------------------------------------------------------------------
// Print the size of the given arrangement.
//
template <typename Arrangement>
void print_arrangement_size(const Arrangement& arr) {
  std::cout << "The arrangement size:\n"
            << "   |V| = " << arr.number_of_vertices()
            << ",  |E| = " << arr.number_of_edges()
            << ",  |F| = " << arr.number_of_faces() << std::endl;
}

//-----------------------------------------------------------------------------
// Print the size of the given unbounded arrangement.
//
template <typename Arrangement>
void print_unbounded_arrangement_size(const Arrangement& arr) {
  std::cout << "The arrangement size:\n"
            << "   |V| = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  |E| = " << arr.number_of_edges()
            << ",  |F| = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces() << " unbounded)\n\n";
}

#endif
