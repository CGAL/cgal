#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>          Triangulation;
typedef Triangulation::Vertex_handle      Vertex_handle;
typedef Triangulation::Face_handle        Face_handle;
typedef Triangulation::All_faces_iterator All_faces_iterator;
typedef Triangulation::All_edges_iterator All_edges_iterator;
typedef Triangulation::Point              Point;

int main() {

  Point p(0,0), q(1,0);

  Triangulation t;

  Vertex_handle inf = t.infinite_vertex();
  Face_handle fh = inf->face();
  assert(fh->vertex(0) == inf);
  assert(fh->vertex(1) == Vertex_handle());
  assert(fh->vertex(2) == Vertex_handle());

  assert(t.all_faces_begin() == t.all_faces_end());
  assert(t.all_edges_begin() == t.all_edges_end());

  t.insert(p);
  Vertex_handle pvh = t.finite_vertices_begin();
  Face_handle pfh = pvh->face();
  assert(pfh->neighbor(0) == fh);

  t.insert(q);

  assert(t.infinite_vertex()->face() == fh);

  assert( (fh->vertex(0) == inf) || (fh->vertex(1) == inf) );

  std::cout << "After the insertion of the second point" <<std::endl;
  std::cout << "|V| = " << t.number_of_vertices() << std::endl;
  std::cout << "|F| = " << t.number_of_faces() << std::endl;


  // Even now we have not really faces
  assert(t.all_faces_begin() == t.all_faces_end());

  for (Triangulation::All_edges_iterator it = t.all_edges_begin();
      it != t.all_edges_end();  ++it) {
    Face_handle fh = it->first;
    Vertex_handle v0 = fh->vertex(0);
    Vertex_handle v1 = fh->vertex(1);
    std::cout << "Edge: ";
    if (v0 == inf) {std::cout << "inf -- ";}else{ std::cout<< v0->point() << " -- ";}
    if (v1 == inf) {std::cout << "inf\n";}else{ std::cout<< v1->point() << std::endl;}
  }

  std::cout << "Edge traversal by hand" << std::endl;
  Face_handle done = fh;
  do {
      assert(fh->vertex(2) == Vertex_handle());
      assert(fh->neighbor(2) == Face_handle());
      Vertex_handle v0 = fh->vertex(0);
      Vertex_handle v1 = fh->vertex(1);
      std::cout << "Edge: ";
      if (v0 == inf) {std::cout << "inf -- ";}else{ std::cout<< v0->point() << " -- ";}
      if (v1 == inf) {std::cout << "inf\n";}else{ std::cout<< v1->point() << std::endl;}
      fh = fh->neighbor(0);
  } while (fh != done);

  std::cout << std::endl;

  return 0;
}
