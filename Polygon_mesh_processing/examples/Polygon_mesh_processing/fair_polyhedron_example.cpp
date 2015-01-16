#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/fair.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef Polyhedron::Vertex_handle   Vertex_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

// extract vertices which are at most k (inclusive) far from vertex v
std::vector<Vertex_handle> extract_k_ring(Vertex_handle v, int k)
{
  std::map<Vertex_handle, int>  D;
  std::vector<Vertex_handle>    Q;
  Q.push_back(v); D[v] = 0;
  std::size_t current_index = 0;

  int dist_v;
  while( current_index < Q.size() && (dist_v = D[ Q[current_index] ]) < k ) {
    v = Q[current_index++];

    Halfedge_around_vertex_circulator e(v->vertex_begin()), e_end(e);
    do {
      Vertex_handle new_v = e->opposite()->vertex();
      if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        Q.push_back(new_v);
      }
    } while(++e != e_end);
  }
  return Q;
}

int main() {
  Polyhedron poly;
  std::ifstream input("data/max.off");
  if ( !input || !(input >> poly) || poly.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  Vertex_iterator v = poly.vertices_begin();
  std::advance(v, 8286);
  const std::vector<Vertex_handle>& region = extract_k_ring(v, 45);

  bool success = CGAL::Polygon_mesh_processing::fair(poly,
    region.begin(), region.end());
  std::cout << "Is fairing successful: " << success << std::endl;

  std::ofstream faired_off("data/faired.off");
  faired_off << poly;
  faired_off.close();
}
