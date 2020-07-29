#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>

#include <fstream>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef Polyhedron::Vertex_handle   Vertex_handle;

// extract vertices which are at most k (inclusive)
// far from vertex v in the graph of edges
void extract_k_ring(Vertex_handle v,
                    int k,
                    std::vector<Vertex_handle>& qv)
{
  std::map<Vertex_handle, int>  D;
  qv.push_back(v);
  D[v] = 0;
  std::size_t current_index = 0;

  int dist_v;
  while (current_index < qv.size() && (dist_v = D[qv[current_index]]) < k)
  {
    v = qv[current_index++];

    Polyhedron::Halfedge_around_vertex_circulator e(v->vertex_begin()), e_end(e);
    do {
      Vertex_handle new_v = e->opposite()->vertex();
      if (D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        qv.push_back(new_v);
      }
    } while (++e != e_end);
  }
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/blobby.off";
  std::ifstream input(filename);

  Polyhedron poly;
  if ( !input || !(input >> poly) || poly.empty()
              || !CGAL::is_triangle_mesh(poly)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  std::vector<Polyhedron::Facet_handle>  new_facets;
  std::vector<Vertex_handle> new_vertices;

  CGAL::Polygon_mesh_processing::refine(poly,
                  faces(poly),
                  std::back_inserter(new_facets),
                  std::back_inserter(new_vertices),
                  CGAL::Polygon_mesh_processing::parameters::density_control_factor(2.));

  std::ofstream refined_off("refined.off");
  refined_off.precision(17);
  refined_off << poly;
  refined_off.close();
  std::cout << "Refinement added " << new_vertices.size() << " vertices." << std::endl;

  Polyhedron::Vertex_iterator v = poly.vertices_begin();
  std::advance(v, 82/*e.g.*/);
  std::vector<Vertex_handle> region;
  extract_k_ring(v, 12/*e.g.*/, region);

  bool success = CGAL::Polygon_mesh_processing::fair(poly, region);
  std::cout << "Fairing : " << (success ? "succeeded" : "failed") << std::endl;

  std::ofstream faired_off("faired.off");
  faired_off.precision(17);
  faired_off << poly;
  faired_off.close();

  return 0;
}
